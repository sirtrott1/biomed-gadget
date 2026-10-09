"""Measure complete table admission on retained original input bytes."""
import argparse,json,os,statistics,subprocess,sys,time
from pathlib import Path
import replay


def measure_command(command, stdout=subprocess.DEVNULL):
    # Keep the measuring parent small and independent of input/binary reads.
    # wait4 also sees pre-exec child residency, so a large orchestrator's heap
    # must not become the apparent peak of a small importer.
    if not hasattr(os, 'wait4'):
        raise RuntimeError('per-process resource measurement requires POSIX wait4')
    start=time.perf_counter()
    child=subprocess.Popen(list(map(str,command)),stdout=stdout)
    _,status,usage=os.wait4(child.pid,0)
    child.returncode=os.waitstatus_to_exitcode(status)
    seconds=time.perf_counter()-start
    if child.returncode:
        raise subprocess.CalledProcessError(child.returncode,child.args)
    rss=usage.ru_maxrss / 1024 if sys.platform=='darwin' else usage.ru_maxrss
    return dict(seconds=seconds,peak_rss_kib=rss)


def measure_child(executable, packet, destination):
    return measure_command([executable,'import-table','--spec',packet,'--out',destination])


def measure(binary,spec_path,out,repeat=3,baseline=None,label='source'):
    if repeat <= 0:
        raise ValueError('repeat must be positive')
    out.mkdir(parents=True,exist_ok=False)
    original=json.loads(spec_path.read_text())
    source=out/'original';source.mkdir()
    for kind in ['source','dictionary']:
        relative=Path(original[kind]['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('input artifacts must stay relative to their specification')
        target=source/relative;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((spec_path.parent/relative).read_bytes())
    replay.save(source/'spec.json',original)
    profiles=[('mapped-quoted',binary,'quoted_v1'),('mapped-spans',binary,'document_spans_v2')]
    if baseline:
        profiles.insert(0,('previous-quoted',baseline,None))
    records=[]
    for profile,executable,kind in profiles:
        samples=[]
        for i in range(repeat):
            declaration=dict(original)
            if kind:
                declaration['evidence_profile']=kind
            packet=replay.save(source/f'measured-{profile}.json',declaration)
            destination=out/f'{profile}-{i}'
            with (out/f'{profile}-{i}.stderr').open('w') as error:
                worker=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--record-child','--binary',str(executable),'--spec',str(packet),'--out',str(destination)],stdout=subprocess.PIPE,stderr=error,text=True,check=True)
            observed=json.loads(worker.stdout)
            seconds,rss=observed['seconds'],observed['peak_rss_kib']
            samples.append(dict(seconds=float(seconds),peak_rss_kib=int(rss),native_sha256=replay.digest((destination/'study.rfcxs').read_bytes()),declaration_bytes=(destination/'study.json').stat().st_size))
        records.append(dict(case=label,profile=profile,runs=len(samples),median_seconds=statistics.median(s['seconds'] for s in samples),median_peak_rss_kib=statistics.median(s['peak_rss_kib'] for s in samples),binary_sha256=replay.digest(executable.read_bytes()),samples=samples))
    assert len({s['native_sha256'] for x in records for s in x['samples']})==1
    report=dict(status='PASS',scope='Complete importer command, including native construction/admission and catalogue output; named index counts are not RSS; sequential samples do not establish a universal scaling or kernel benchmark',measurements=records)
    replay.save(out/'performance.json',report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',required=True,type=Path);p.add_argument('--spec',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--repeat',type=int,default=3);p.add_argument('--baseline-binary',type=Path);p.add_argument('--label',default='source');p.add_argument('--record-child',action='store_true',help=argparse.SUPPRESS);a=p.parse_args()
    if a.record_child:
        print(json.dumps(measure_child(a.binary,a.spec,a.out)))
        return
    print(json.dumps(measure(a.binary.resolve(),a.spec.resolve(),a.out.resolve(),a.repeat,a.baseline_binary.resolve() if a.baseline_binary else None,a.label),indent=2))


if __name__=='__main__':main()
