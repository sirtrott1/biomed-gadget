"""Complete native history publication, source ownership and exact map contracts."""
from pathlib import Path
import copy,importlib.util,json,shutil,subprocess,sys

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def main(binary,fixture,consumer,root):
    repository=Path(__file__).resolve().parents[1];replay=load(repository/'examples/bioinformatics/versioned-history/replay.py','history_replay')
    verify=load(repository/'examples/bioinformatics/versioned-history/verify.py','history_verify')
    sys.path.insert(0,str(repository/'examples/bioinformatics/factored-temporal'))
    from verify import verify_pair
    shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    independent=0;native=[]
    for series in [False,True]:
        case=root/('series' if series else 'publications');query,result=replay.authored(binary,fixture,case,series)
        independent+=verify.authored(case)['checks'];requests={}
        for name in ['certificate','relations','completion']:
            declared=json.loads((case/(name+'-query.json')).read_text());full=json.loads((case/(name+'.json')).read_text());requests[name+'-enumerated']=declared
            factored=copy.deepcopy(declared);factored['map_output']='factored';other=replay.evaluate(binary,case,factored,name+'-factored');requests[name+'-factored']=factored
            matched=verify_pair(full,other);independent+=matched['complete_map_coordinates']+matched['source_binding_coordinates']
        # Every boundary and adjoint position has its own actual grade axes.
        ordered=copy.deepcopy(query);ordered['endpoints']=[dict(publication=id,context='all',field='r') for id in ['p0','p1','p2']]
        ordered['word']=[dict(kind='boundary'),dict(kind='adjoint'),dict(kind='green',**{'lambda':'2/3'})]
        full=replay.evaluate(binary,case,ordered,'ordered');requests['ordered-enumerated']=ordered
        factored=copy.deepcopy(ordered);factored['map_output']='factored';other=replay.evaluate(binary,case,factored,'ordered-factored');requests['ordered-factored']=factored
        matched=verify_pair(full,other);independent+=matched['complete_map_coordinates']+matched['source_binding_coordinates']
        diagonal=copy.deepcopy(ordered)
        diagonal['word'].insert(0,dict(kind='diagonal_field',fields=['r','r','r']))
        full=replay.evaluate(binary,case,diagonal,'diagonal');requests['diagonal-enumerated']=diagonal
        factored=copy.deepcopy(diagonal);factored['map_output']='factored';other=replay.evaluate(binary,case,factored,'diagonal-factored');requests['diagonal-factored']=factored
        matched=verify_pair(full,other);independent+=matched['complete_map_coordinates']+matched['source_binding_coordinates']
        # Independent namespace lineages require explicit maps, even with matching keys.
        foreign=copy.deepcopy(query)
        for lineage in foreign['lineages']:lineage['kind']='independent_sources'
        ends=result['result']['endpoints']
        for t,step in enumerate(foreign['correspondences']):
            for k,declared in enumerate(step):
                declared['kind']='explicit';declared['entries']=[]
                for a in ends[t]['declaration']['source_grades'][k]:
                    for b in ends[t+1]['declaration']['source_grades'][k]:
                        if a['key']==b['key']:declared['entries'].append(dict(row=b['id'],column=a['id'],value='1'))
        replay.evaluate(binary,case,foreign,'independent');requests['independent']=foreign
        path=replay.save(case/'native-queries.json',requests)
        result=subprocess.run([str(consumer),str(case/'study.json'),str(path),str(case/'native-results')],capture_output=True,text=True)
        (case/'native.log').write_text(result.stdout+result.stderr);assert result.returncode==0,result.stderr;native.append(json.loads(result.stdout))
    receipt=dict(status='PASS',independent_checks=independent,native_cases=native)
    (root/'contracts.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':main(*map(lambda value:Path(value).resolve(),sys.argv[1:]))
