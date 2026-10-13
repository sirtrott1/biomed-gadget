"""Independent exact endpoint, sector, metric and complete native map contracts."""
from pathlib import Path
import copy,hashlib,json,shutil,subprocess,sys
import sympy as S

def main(binary,consumer,fixture,root):
    shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples/bioinformatics/factored-temporal'))
    from maps import Expansion
    B=[S.Matrix([[-1,0],[1,0],[0,0]]),S.Matrix([[-1,0,0],[1,-1,0],[0,1,0]]),S.Matrix([[-1,0],[1,0],[0,0]])]
    C=[S.zeros(2,0),S.Matrix([[0],[0],[1]]),S.zeros(2,0)]
    coordinates=[[0,2],[0,1,2],[0,2]];names=['baseline','followup','final']
    xs=[S.Matrix([1,2]),S.Matrix([3,1,6]),S.Matrix([13,-7])];checks=0;cases=[]
    def projector(frame,metric):
        independent=frame.columnspace()
        if not independent:return S.zeros(metric.rows)
        f=S.Matrix.hstack(*independent);return f*(f.T*metric*f).inv()*f.T*metric
    for profile in ['identity','weighted']:
        case=root/profile;study=case/'study'
        subprocess.run([str(fixture),str(study)],capture_output=True,text=True,check=True)
        d=json.loads((study/'study.json').read_text())
        if profile=='weighted':
            for t,context in enumerate(d['contexts']):
                context['metrics']=[dict(kind='diagonal',values=[str(1+t+i) for i in range(3)]),
                    dict(kind='diagonal',values=[str(2+2*t+i) for i in range(3)]),dict(kind='diagonal',values=[str(7+t)])]
            (study/'study.json').write_text(json.dumps(d))
        declared={context['id']:context for context in d['contexts']}
        metrics=[];operators=[];sectors=[]
        for t in range(3):
            declarations=declared[names[t]]['metrics']
            m0=S.diag(*[S.Rational(x) for x in declarations[0].get('values',['1']*3)])
            m1=S.diag(*[S.Rational(declarations[1].get('values',['1']*3)[i]) for i in coordinates[t]])
            m2=S.Rational(declarations[2].get('values',['1'])[0]);metrics.append(m1)
            gradient=m1.inv()*B[t].T*m0;curl=C[t];g=projector(gradient,m1);c=projector(curl,m1)
            sectors.append([g,c,S.eye(m1.rows)-g-c]);operators.append(gradient*B[t]+C[t]*(S.eye(C[t].cols)/m2)*C[t].T*m1)
        requests={}
        for t in range(2):
            j=S.zeros(len(coordinates[t+1]),len(coordinates[t]))
            for column,original in enumerate(coordinates[t]):
                if original in coordinates[t+1]:j[coordinates[t+1].index(original),column]=1
            shift=S.Rational(3,2);old_g=(S.eye(j.cols)+shift*operators[t]).inv();new_g=(S.eye(j.rows)+shift*operators[t+1]).inv()
            innovation=xs[t+1]-j*xs[t]
            expected=dict(correspondence=j,metric_pullback_defect=j.T*metrics[t+1]*j-metrics[t],
                boundary_defect=B[t+1]*j-B[t],old_transported=j*xs[t],field_innovation=innovation,
                realization_term=S.zeros(j.rows,1),innovation_term=new_g*innovation,
                operator_term=(new_g*j-j*old_g)*xs[t],total=new_g*xs[t+1]-j*old_g*xs[t],
                green_operator_defect=operators[t+1]*j-j*operators[t],green_response_defect=new_g*j-j*old_g)
            for mode in ['enumerated','factored']:
                query=dict(operation='structure_transition',context=names[t],other_context=names[t+1],
                    field=names[t]+'.x',new_field=names[t+1]+'.x',**{'lambda':'3/2'},map_output=mode)
                name=str(t)+'-'+mode;requests[name]=query;path=case/(name+'.json');path.write_text(json.dumps(query))
                process=subprocess.run([str(binary),'query','--study',str(study/'study.json'),'--query',str(path)],capture_output=True,text=True)
                assert process.returncode==0,process.stderr;actual=json.loads(process.stdout)['result']
                expansion=Expansion(actual['factored_maps']) if mode=='factored' else None
                def matrix(value):
                    if 'expression' in value:return S.Matrix(expansion.map(value))
                    if 'values' in value:return S.Matrix([[S.Rational(x) for x in row] for row in value['values']])
                    out=S.zeros(value['rows'],value['columns'])
                    for item in value['entries']:out[item['row'],item['column']]=S.Rational(item['value'])
                    return out
                for key,value in expected.items():
                    selected=matrix(actual[key]) if key in ['correspondence','metric_pullback_defect','boundary_defect','green_operator_defect','green_response_defect'] else S.Matrix([S.Rational(x) for x in actual[key]])
                    assert selected==value,(profile,name,key,selected,value);checks+=value.rows*value.cols
                for label,index,green in [('old',t,old_g),('new',t+1,new_g)]:
                    for sector,p in zip(['gradient','curl','harmonic'],sectors[index]):
                        field=p*xs[index];assert S.Matrix([S.Rational(x) for x in actual[label][sector]])==field
                        assert S.Rational(actual[label][sector+'_quadrance'])==(field.T*metrics[index]*field)[0];checks+=field.rows+1
                    assert S.Matrix([S.Rational(x) for x in actual[label]['green']])==green*xs[index];checks+=green.rows
                for a in range(3):
                    for b in range(3):
                        block=sectors[t+1][b]*j*sectors[t][a]
                        assert matrix(actual['sector_maps'][a][b])==block
                        assert S.Matrix([S.Rational(x) for x in actual['sector_transport'][a][b]])==block*xs[t];checks+=block.rows*(block.cols+1)
        path=case/'queries.json';path.write_text(json.dumps(requests))
        result=subprocess.run([str(consumer),str(study/'study.json'),str(path),str(case/'native')],capture_output=True,text=True)
        assert result.returncode==0,result.stderr;cases.append(json.loads(result.stdout))
        for edit in [{'lambda':'0'},{'field':'baseline.partial'},{'other_context':'global'}]:
            bad={**requests['0-enumerated'],**edit};path=case/'refused.json';path.write_text(json.dumps(bad))
            result=subprocess.run([str(binary),'query','--study',str(study/'study.json'),'--query',str(path)],capture_output=True,text=True)
            assert result.returncode!=0 and result.stderr;checks+=1
    case=root/'129-coordinates';inputs=case/'inputs';inputs.mkdir(parents=True)
    table='id,before,after\n'+''.join(f'r{i},{i%7},{(i+3)%11}\n' for i in range(129))
    dictionary='Every authored row is retained. Two fields on 129 source coordinates; declaration stages and identity metrics.\n'
    (inputs/'values.csv').write_text(table);(inputs/'dictionary.txt').write_text(dictionary)
    def source_file(name):return dict(path=name,sha256='sha256:'+hashlib.sha256((inputs/name).read_bytes()).hexdigest())
    columns=[dict(id=name,channel='observation',quantity='authored_value',units='units',value_type='quantitative',dictionary_start=0,dictionary_end=len(dictionary)) for name in ['before','after']]
    columns.insert(0,dict(id='id',channel='observation',quantity='identity',units='dimensionless',value_type='categorical',dictionary_start=0,dictionary_end=len(dictionary)))
    specification=dict(schema='gadget.table.import.v1',id='complete-129-coordinates',source_namespace='913013',
        source={**source_file('values.csv'), 'delimiter':',','header':True},dictionary=source_file('dictionary.txt'),
        columns=columns,missing_tokens=['NA'],relation_columns=['id'],identity_column='id',
        question='What does the declared full comparison retain?',population='All 129 authored rows',selection='All rows',
        limitations=['Authored algebra contract'],source_url='authored://complete-coordinates',attribution='GADGET method fixture',license='Authored fixture')
    path=inputs/'import.json';path.write_text(json.dumps(specification));study=case/'study'
    imported=subprocess.run([str(binary),'import-table','--spec',str(path),'--out',str(study)],capture_output=True,text=True)
    assert imported.returncode==0,imported.stderr
    metadata=json.loads((study/'study.json').read_text());context=metadata['contexts'][0]
    for name,position in [('before','0'),('after','1')]:
        declared=copy.deepcopy(context);declared['id']=name
        declared['time']=dict(axis='declared_comparison_stages',unit='ordinal_stage',position=position,population='All 129 authored rows')
        metadata['contexts'].append(declared)
    for field in metadata['fields']:
        if field['id'] in ['all.before','all.after']:field['context']=field['id'].removeprefix('all.')
    (study/'study.json').write_text(json.dumps(metadata));queries={}
    for mode in ['enumerated','factored']:
        query=dict(operation='structure_transition',context='before',other_context='after',field='all.before',new_field='all.after',**{'lambda':'2/3'},map_output=mode)
        queries[mode]=query;path=case/(mode+'.json');path.write_text(json.dumps(query))
        result=subprocess.run([str(binary),'query','--study',str(study/'study.json'),'--query',str(path)],capture_output=True,text=True)
        assert result.returncode==0,result.stderr;actual=json.loads(result.stdout)['result'];expansion=Expansion(actual['factored_maps']) if mode=='factored' else None
        for a in range(3):
            for b in range(3):
                descriptor=actual['sector_maps'][a][b]
                block=expansion.map(descriptor) if expansion else S.Matrix([[S.Rational(x) for x in row] for row in descriptor['values']])
                assert block==(S.eye(129) if a==b==1 else S.zeros(129));checks+=129*129
        assert list(map(S.Rational,actual['total']))==[S.Rational(3,5)*((i+3)%11-i%7) for i in range(129)];checks+=129
    path=case/'queries.json';path.write_text(json.dumps(queries));result=subprocess.run([str(consumer),str(study/'study.json'),str(path),str(case/'native')],capture_output=True,text=True)
    assert result.returncode==0,result.stderr;cases.append(json.loads(result.stdout))
    receipt=dict(status='PASS',independent_checks=checks,native_cases=cases)
    (root/'contracts.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))

if __name__=='__main__':main(*(Path(x).resolve() for x in sys.argv[1:]))
