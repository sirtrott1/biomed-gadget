"""Native comparison owners with independent exact count and control oracles."""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json,subprocess,sys,shutil
from research_cli import projection,dot


def run(binary,consumer,root):
    if root.exists():shutil.rmtree(root)
    root.mkdir(parents=True)
    checks=0
    def invoke(*args):
        p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
        assert p.returncode==0,p.stderr
        return json.loads(p.stdout)
    def import_rows(name,header,rows,kinds):
        case=root/name;case.mkdir();dictionary='Authored original rows, complete declared units and missing support.\n'
        raw='\n'.join(['\t'.join(header),*['\t'.join(row) for row in rows]])+'\n'
        (case/'table.tsv').write_text(raw);(case/'dictionary.txt').write_text(dictionary)
        columns=[]
        for i,column in enumerate(header):
            channel,quantity,units,value_type=kinds[i]
            columns.append(dict(id=column,channel=channel,quantity=quantity,units=units,value_type=value_type,dictionary_start=0,dictionary_end=len(dictionary)))
        spec=dict(schema='gadget.table.import.v1',id=name,source_namespace=str(915011+len(header)),source=dict(path='table.tsv',sha256='sha256:'+hashlib.sha256(raw.encode()).hexdigest(),delimiter='\t',header=True),dictionary=dict(path='dictionary.txt',sha256='sha256:'+hashlib.sha256(dictionary.encode()).hexdigest()),columns=columns,identity_column=header[0],missing_tokens=['NA'],relation_columns=[header[1]],question='Authored exact native comparison contracts',population='all original authored rows',selection='all rows',limitations=['authored nonclinical fixture'],source_url='synthetic://native-comparison',attribution='GADGET contracts',license='synthetic')
        (case/'import.json').write_text(json.dumps(spec));invoke('import-table','--spec',case/'import.json','--out',case/'study')
        return case,dictionary
    categorical=('observation','identity','dimensionless','categorical')
    count=('E','read_count','counts','quantitative')
    pair_rows=[['gA','gene','4','12','3','7','-1','0.01'],['gB','gene','0','0','0','0','0','NA'],['gC','gene','5','1','1','9','1','0.1'],['gD','gene','2','2','1','NA','0','NA']]
    case,dictionary=import_rows('paired',['gene','kind','a0','a1','b0','b1','effect','padj'],pair_rows,[categorical,categorical,*[count]*4,('observation','model_log2_fold_change','dimensionless','quantitative'),('observation','adjusted_p_value','dimensionless','quantitative')])
    paired=dict(operation='paired_field_contrasts',pairs=[dict(id='A',left='all.a1',right='all.a0',left_scale='3/2',right_scale='1'),dict(id='B',left='all.b1',right='all.b0',left_scale='2',right_scale='1/2')],selection=[],support='joint_reported',normalization=dict(method='declared_scales',evidence=[dict(document='dictionary',start=0,end=len(dictionary),text=dictionary)]),output_units='scaled_counts',workers=1,benchmark=dict(effect='all.effect',adjusted_p='all.padj',alpha='0.05'))
    query=case/'query.json';query.write_text(json.dumps(paired));r=invoke('query','--study',case/'study/study.json','--query',query)['result'];assert r['reported']==3 and r['missing']==1
    for row,actual in zip(pair_rows,r['features']):
        left=[F(row[3])*F(3,2),F(row[5])*2];right=[F(row[2]),F(row[4])*F(1,2)];delta=[a-b for a,b in zip(left,right)];mean=sum(delta)/2;residual=[x-mean for x in delta]
        for key,expected in [('left',left),('right',right),('delta',delta),('consensus',[mean]*2),('residual',residual)]:assert list(map(F,actual[key]))==expected;checks+=len(expected)
        assert F(actual['heterogeneity_quadrance'])==dot(residual,residual);checks+=1
    assert r['features'][1]['mean_ratio'] is None and r['features'][1]['benchmark']['adjusted_p'] is None
    paired_queries={'serial':paired,'parallel':{**paired,'workers':2}}
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    huge=copy.deepcopy(paired);huge['pairs'][0]['left_scale']=str(1<<65537);paired_queries['large_exact_scale']=huge
    (case/'queries.json').write_text(json.dumps(paired_queries));p=subprocess.run([str(consumer),str(case/'study/study.rfcxs'),str(case/'queries.json'),str(case/'native'),'study'],capture_output=True,text=True);assert p.returncode==0,p.stderr;print(p.stdout.strip())
    rows=[]
    for i in range(16):
        size=i+1;y=i%2;x=17+3*size+2*y;z=17+3*size
        rows.append([f'person-{i}','batch',str(size),str(size),str(y),str(x),str(z),'NA' if i==1 else str(x),'NA' if i==0 else str(y)])
    quantitative=('observation','measurement','synthetic_measurement','quantitative');indicator=('observation','membership','dimensionless','indicator')
    case,_=import_rows('controlled',['person','group','size','copy','y','x','z','missing_x','partial_y'],rows,[categorical,categorical,quantitative,quantitative,indicator,quantitative,quantitative,quantitative,indicator])
    controlled=dict(operation='controlled_pairings',fields=['all.x','all.z'],endpoint='all.y',controls=['all.size','all.copy'],selection=[],support='joint_reported',metric='identity',intercept=True)
    cases={'dependent_controls':controlled,'parallel':{**controlled,'workers':2},'empty_gram':{**controlled,'gram_fields':[]},'reverse_gram':{**controlled,'gram_fields':['all.z','all.x']},'selected':{**controlled,'selection':['all.y']},'missing':{**controlled,'fields':['all.missing_x'],'endpoint':'all.partial_y'},'intercept':{**controlled,'controls':[]}}
    for name,q in cases.items():
        path=case/(name+'.json');path.write_text(json.dumps(q));r=invoke('query','--study',case/'study/study.json','--query',path)['result'];selected=[i for i,row in enumerate(rows) if (not q['selection'] or row[4]=='1') and not (name=='missing' and i in [0,1])]
        columns=[[F(1)]*len(selected),*[list(map(F,[rows[i][2] for i in selected])) for _ in q['controls']]];endpoint=[F(rows[i][4]) for i in selected];py=projection(columns,endpoint);ry=[a-b for a,b in zip(endpoint,py)]
        assert list(map(F,r['endpoint_residual']))==ry and r['coordinates']==selected;checks+=len(selected)*2
        for actual in r['fields']:
            col=5 if actual['field'] in ['all.x','all.missing_x'] else 6;x=[F(rows[i][col]) for i in selected];p=projection(columns,x);residual=[a-b for a,b in zip(x,p)]
            assert list(map(F,actual['residual']))==residual and F(actual['pairing'])==dot(residual,ry);checks+=len(selected)+1
        assert r['control_kernel']['columns']==len(q['controls'])+1-r['control_rank'];checks+=1
    (case/'queries.json').write_text(json.dumps(cases));p=subprocess.run([str(consumer),str(case/'study/study.rfcxs'),str(case/'queries.json'),str(case/'native'),'study'],capture_output=True,text=True);assert p.returncode==0,p.stderr;print(p.stdout.strip())
    print(json.dumps(dict(status='PASS',independent_checks=checks,paired_queries=len(paired_queries),controlled_queries=len(cases),complete_source_rows=20,large_scale_bits=65538)))

if __name__=='__main__':run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),Path(sys.argv[3]).resolve())
