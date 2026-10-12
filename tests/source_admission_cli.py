"""Versioned span/quote, retained-source, budget and malformed-table contracts."""
import copy,json,random,shutil,subprocess,sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'examples/bioinformatics/source-admission'))
import replay
binary,contracts,root=map(lambda x:Path(x).resolve(),sys.argv[1:])
shutil.rmtree(root,ignore_errors=True);root.mkdir()
spec_path=replay.authored(root/'original',7)
result=replay.compare(binary,spec_path,root/'comparison')
original=json.loads(spec_path.read_text())
spans=root/'comparison/spans/study'
corpus=[]
def representable(value):
 if isinstance(value,dict):return all('\0' not in key and representable(item) for key,item in value.items())
 if isinstance(value,list):return all(map(representable,value))
 return True
def json_case(source):
 try:
  value=json.loads(source);admit=representable(value)
 except (ValueError,RecursionError):
  value=None;admit=False
 case=dict(source=source,admit=admit)
 if admit:case['value']=value
 corpus.append(case)
# An independent standard decoder decides every valid token's actual value and
# whether any nested member name contains an unrepresentable NUL.
for token in ['null','true','false','0','-17','9223372036854775807','""','"\\u0000"','[]','{}']:
 for whitespace in ['', ' ', '\t\r\n']:
  json_case(whitespace+token+whitespace)
for slashes in range(1,17):
 for whitespace in ['', ' ', '\t\r\n']:
  json_case('{"prefix'+'\\'*slashes+'u0000suffix"'+whitespace+':"value"}')
for name in ['', 'schema', 'β', 'quote"slash\\', 'name\\u0000']:
 for position in range(len(name)+1):
  key=name[:position]+'\0'+name[position:]
  for wrap in [lambda x:x,lambda x:[x],lambda x:dict(outer=[dict(inner=x)])]:
   json_case(json.dumps(wrap({key:'original'}),ensure_ascii=True))
 for value in ['\0', 'before\0after', '\\u0000', '" : "\\u0000', 'β🙂\0tail']:
  json_case(json.dumps({name:value},ensure_ascii=True))
for source in ['{"schema\\u0000suffix":"control"}',
               '{"\\u0073chema\\u0000suffix":"control"}',
               '{"schema":"original","schema\\u0000suffix":"replacement"}',
               '{"value":"\\\\u0000","other":[{"name":"\\u0000"}]}',
               '{"quote\\\" : \\"\\u0000suffix":1}',
               '{"value":"raw\0bytes"}', '{"raw\0key":1}',
               '{"value":1}\0', '{"value":1}\0{"other":2}',
               '', 'tru', '1e', '[', '{', '{"name":"\\u000',
               '{"name":"\\z"}', '{"name":1} false', '/* comment */ {}']:
 json_case(source)
for depth in [31,32,33,64,96]:
 for leaf in ['{"name":"\\u0000"}','{"name\\u0000suffix":1}']:
  json_case('['*depth+leaf+']'*depth)
rng=random.Random(971031)
names=['', 'schema', 'value', 'β', '🙂', 'quoted"key', '\\u0000', 'nul\0name']
values=[None,True,False,0,-17,971031,'', '\0','before\0after','quoted" : \\"text','β🙂','\\u0000']
def authored(depth):
 if depth==0 or rng.randrange(3)==0:return rng.choice(values)
 if rng.randrange(2):return [authored(depth-1) for _ in range(rng.randrange(4))]
 return {key:authored(depth-1) for key in rng.sample(names,rng.randrange(5))}
for _ in range(512):
 json_case(json.dumps(authored(rng.randrange(1,7)),ensure_ascii=rng.choice([True,False]),
                      separators=rng.choice([(',',':'),(', ',': ')])))
corpus_path=replay.save(root/'json-corpus.json',corpus)
subprocess.run([str(contracts),str(spans/'study.json'),str(corpus_path)],check=True)
refusals=[]
def refuse(args):
 p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
 assert p.returncode,(args,p.stdout)
 refusals.append(p.stderr.strip())
for i,change in enumerate([dict(evidence_profile='invented'),dict(intake=dict(max_source_bytes=1)),dict(intake=dict(max_record_bytes=1)),dict(intake=dict(max_token_offsets=1)),dict(intake=dict(max_token_offsets=-1)),dict(intake=dict(other=1))]):
 q={**original,**change};path=replay.save(spec_path.parent/f'refuse-{i}.json',q);dest=root/f'bad-{i}'
 refuse(['import-table','--spec',path,'--out',dest]);assert not dest.exists()
 assert not list(root.glob('.gadget-table-*'))
# Exact caller thresholds admit the whole source rather than truncate it.
raw=(spec_path.parent/'table.tsv').read_bytes();q={**original,'evidence_profile':'document_spans_v2','intake':dict(max_source_bytes=len(raw),max_record_bytes=max(map(len,raw.splitlines())),max_token_offsets=7*5)}
p=replay.save(spec_path.parent/'exact-budget.json',q)
replay.run(binary,'import-table','--spec',p,'--out',root/'exact-budget')
# Legacy declarations still require quotes; version-two supplied quotes are checked.
declaration=json.loads((spans/'study.json').read_text())
for i,change in enumerate(['legacy','quote','unknown','bounds','empty','utf8','hash']):
 d=copy.deepcopy(declaration);span=d['cells'][0]['evidence'][0]
 if change=='legacy':d['schema']='gadget.scientific.study.v1'
 if change=='quote':span['text']='invented quotation'
 if change=='unknown':span['document']='not-admitted'
 if change=='bounds':span['end']=len(raw)+1
 if change=='empty':span['end']=span['start']
 if change=='utf8':span['start']=raw.index('β'.encode())+1
 if change=='hash':d['documents'][0]['sha256']='sha256:'+'0'*64
 path=replay.save(spans/f'bad-{i}.json',d)
 refuse(['source-index','--study',path])
# Quoted input remains admissible in v2 and matches actual document bytes.
d=copy.deepcopy(declaration);span=d['cells'][0]['evidence'][0];span['text']=raw[span['start']:span['end']].decode()
replay.run(binary,'source-index','--study',replay.save(spans/'mixed.json',d))
# Unquoted dialect has explicit refusal rather than silently discarding records.
for i,bad in enumerate([b'identity\tgroup\tx\ty\tz\n',raw+b'\n',raw.replace(b'0/3',b'"0/3"',1),raw.replace(b'0/3',b'0\0/3',1),raw.replace(b'0/3',b'0\t/3',1)]):
 (spec_path.parent/'bad.tsv').write_bytes(bad);q={**original,'source':{**original['source'],'path':'bad.tsv','sha256':replay.digest(bad)}}
 dest=root/f'malformed-{i}';refuse(['import-table','--spec',replay.save(spec_path.parent/'bad-table.json',q),'--out',dest]);assert not dest.exists()
# UTF-8 validation covers the whole retained document, including unreferenced bytes.
bad_dictionary=(spec_path.parent/'dictionary.txt').read_bytes()+b'\xff'
(spec_path.parent/'bad-dictionary.txt').write_bytes(bad_dictionary)
q=copy.deepcopy(original);q['evidence_profile']='document_spans_v2';q['dictionary']=dict(path='bad-dictionary.txt',sha256=replay.digest(bad_dictionary))
for column in q['columns']:column['dictionary_end']=10
refuse(['import-table','--spec',replay.save(spec_path.parent/'invalid-utf8-dictionary.json',q),'--out',root/'invalid-utf8-dictionary'])
assert not (root/'invalid-utf8-dictionary').exists()
# External artifact paths and source hash mismatches cannot be admitted.
for i,source in enumerate([{**original['source'],'sha256':'sha256:'+'0'*64},{**original['source'],'path':'../table.tsv'}]):
 dest=root/f'source-{i}';refuse(['import-table','--spec',replay.save(spec_path.parent/'bad-source.json',{**original,'source':source}),'--out',dest]);assert not dest.exists()
assert not list(root.glob('.gadget-table-*'))
refuse(['import-table','--spec',spec_path,'--out',spans])
# Missing final newline and CRLF preserve original byte spans and exact fields.
for name,data in [('no-final-newline',raw.rstrip(b'\n')),('crlf',raw.replace(b'\n',b'\r\n'))]:
 (spec_path.parent/f'{name}.tsv').write_bytes(data);q={**original,'source':{**original['source'],'path':f'{name}.tsv','sha256':replay.digest(data)}}
 replay.compare(binary,replay.save(spec_path.parent/f'{name}.json',q),root/name)
replay.save(root/'contracts.json',dict(**result,refusals=refusals,exact_caller_budgets=True,crlf_and_eof=True,immutable_original_bytes=True,
 json_corpus=dict(cases=len(corpus),admitted=sum(case['admit'] for case in corpus),refused=sum(not case['admit'] for case in corpus),
                 original_bytes_retained=True,native_nul_names_preserved=True)))
print('PASS: mapped retained source, versioned evidence, exact native parity, caller budgets and source refusals')
