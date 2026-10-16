#!/usr/bin/env python3
"""Display complete native assay sections, source meanings and exact defects."""
import argparse
import hashlib
import html
import json
import os
import subprocess
from pathlib import Path
from urllib.parse import quote


def digest(path):
    with path.open('rb') as stream:
        return 'sha256:' + hashlib.file_digest(stream, 'sha256').hexdigest()


def render(binary, operation, output):
    with operation.open('rb') as stream:
        if stream.read(5) != b'RFCXS':
            raise ValueError('An owning native operation is required')
    before = digest(operation)
    output.mkdir(parents=True, exist_ok=False)
    with (output / 'view.json').open('wb') as stdout, (output / 'view.log').open('wb') as stderr:
        process = subprocess.run([str(binary), 'record-view', '--input', str(operation),
                                  '--member', 'result', '--format', 'compact'], stdout=stdout, stderr=stderr)
    if process.returncode:
        raise ValueError('Native reading refused; see view.log')
    if digest(operation) != before:
        raise ValueError('The native operation changed while reading')
    view = json.loads((output / 'view.json').read_text())
    if view['schema'] != 'gadget.assay.result.v1' or view['query']['operation'] != 'assay_sections':
        raise ValueError('The operation must own an assay_sections result')
    result = view['result']
    sections = result['sections']
    if not sections:
        raise ValueError('The native result has no declared sections')
    data = json.dumps(view, ensure_ascii=False).replace('<', '\\u003c')
    native_link = html.escape(quote(os.path.relpath(operation, output), safe='/'), quote=True)
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Native assay sections</title><style>
body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#172e35}
table{border-collapse:collapse;width:100%;margin:1rem 0}th,td{text-align:left;padding:.5rem;border-bottom:1px solid #ccd7d9;overflow-wrap:anywhere}
code{overflow-wrap:anywhere}canvas{width:100%;height:220px}select{max-width:100%;padding:.5rem}pre{white-space:pre-wrap;overflow-wrap:anywhere}
</style><h1>Native assay sections</h1>
<p>Each section retains its original observation identity, source fields, declared interface maps and metrics.
Compatible fields and residuals use those declarations. A declared cross assay cost does not establish biochemical calibration,
clinical patient identity or causal treatment effects. Missing measurements remain absent.</p>
<p><a href="''' + native_link + '''">Complete owning native operation</a> · <a href="view.json">Native display view</a> ·
<a href="receipt.json">Source and byte binding</a></p><p>Native SHA256: <code>''' + before + '''</code></p>
<h2>All declared observations</h2><p>Every section is shown once in its original source order. The horizontal coordinate is an observation index,
not elapsed time or an inferred spatial position. Charts convert exact rational values to floating point for display.</p>
<label>Constraint <select id="constraint"></select></label><p id="chart-note"></p>
<figure><figcaption>Observed signed constraint defect</figcaption><canvas id="defect"></canvas></figure>
<figure><figcaption>Complete residual quadrance in the declared metric</figcaption><canvas id="quadrance"></canvas></figure>
<h2>Exact observation reading</h2><label>Observation <select id="section"></select></label><p id="summary"></p>
<table><thead><tr><th>View</th><th>Source coordinate</th><th>Units</th><th>Observed</th><th>Compatible</th><th>Residual</th></tr></thead><tbody id="fields"></tbody></table>
<h3>Constraint support and exact defects</h3><table><thead><tr><th>Row</th><th>Assessed</th><th>Observed defect</th><th>Compatible defect</th><th>Missing coordinates</th></tr></thead><tbody id="constraints"></tbody></table>
<details><summary>Selected original field identities, support and source meanings</summary><pre id="source"></pre></details>
<details><summary>Complete preparation, coefficient families, maps and metrics</summary><pre id="declarations"></pre></details>
<p>The complete preparation retains every joint free direction and row dependency that was requested.
The source operation owns the original field carriers and byte provenance. Equal coordinate dimensions do not identify independent assays.</p>
<script id="data" type="application/json">''' + data + '''</script><script>
const view=JSON.parse(document.getElementById('data').textContent),r=view.result;
const select=document.getElementById('section'),constraint=document.getElementById('constraint');
const exact=x=>x===null||x===undefined?'Absent':String(x);
function option(select,value,label){const o=document.createElement('option');o.value=String(value);o.textContent=label;select.append(o)}
function row(body,items){const tr=document.createElement('tr');for(const item of items){const td=document.createElement('td');td.textContent=exact(item);tr.append(td)}body.append(tr)}
r.sections.forEach((s,i)=>option(select,i,s.id));r.constraint_coordinates.forEach((c,i)=>option(constraint,i,c.interface+': '+c.observable));
document.getElementById('declarations').textContent=JSON.stringify({source_digest:view.source_digest,views:r.views,interfaces:r.interfaces,metric:r.metric,coboundary:r.coboundary,coordinates:r.coordinates,constraint_coordinates:r.constraint_coordinates,prepared_groups:r.prepared_groups},null,2);
function update(){const s=r.sections[Number(select.value)];document.getElementById('summary').textContent=s.id+'; complete: '+s.complete+'; reconstruction verified: '+s.reconstruction_verified+'; exact residual quadrance: '+exact(s.residual_quadrance)+'; orthogonality: '+exact(s.orthogonality);
const fields=document.getElementById('fields');fields.replaceChildren();for(const f of s.fields)f.ids.forEach((id,i)=>row(fields,[f.view,id,f.meaning.units,f.input[i],f.compatible[i],f.residual[i]]));
const constraints=document.getElementById('constraints');constraints.replaceChildren();s.constraints.forEach(c=>row(constraints,[c.row,c.assessed,c.input_defect,c.compatible_defect,c.missing_coordinates.join(', ')]));document.getElementById('source').textContent=JSON.stringify(s,null,2);}
function number(q){if(q===null||q===undefined)return NaN;const a=String(q).split('/').map(Number);return a.length===2?a[0]/a[1]:a[0]}
function draw(){let unavailable=0;for(const kind of ['defect','quadrance']){const canvas=document.getElementById(kind);canvas.width=1100;canvas.height=220;const ctx=canvas.getContext('2d');ctx.clearRect(0,0,canvas.width,canvas.height);
const values=r.sections.map(s=>kind==='quadrance'?s.residual_quadrance:s.constraints.find(c=>c.row===Number(constraint.value))?.input_defect).map(number);
let max=0;for(const value of values)if(Number.isFinite(value))max=Math.max(max,Math.abs(value));else unavailable++;
const middle=kind==='defect'?110:210,scale=max?(kind==='defect'?100:200)/max:0;ctx.strokeStyle='#ccd7d9';ctx.beginPath();ctx.moveTo(0,middle);ctx.lineTo(canvas.width,middle);ctx.stroke();
values.forEach((value,i)=>{if(!Number.isFinite(value))return;ctx.fillStyle=value<0?'#246f83':'#ad542b';ctx.fillRect(i*canvas.width/values.length,middle-value*scale,Math.max(1,canvas.width/values.length),value*scale||1)});}
document.getElementById('chart-note').textContent=r.sections.length+' original observations; quadrance units: '+r.metric.units+'. Absent or nonfinite display coordinates: '+unavailable+'. Exact values remain in the tables.';}
select.addEventListener('change',update);constraint.addEventListener('change',draw);
for(const id of ['defect','quadrance'])document.getElementById(id).addEventListener('mousemove',event=>{const rect=event.target.getBoundingClientRect();const i=Math.max(0,Math.min(r.sections.length-1,Math.floor((event.clientX-rect.left)*r.sections.length/rect.width)));select.value=String(i);update()});
update();draw();
</script></html>'''
    (output / 'report.html').write_text(page)
    receipt = dict(schema='gadget.assay.sections.display.v1', source_digest=view['source_digest'],
        native_operation=os.path.relpath(operation, output), native_sha256=before,
        display_sha256=digest(output / 'view.json'), sections=len(sections),
        source_coordinates=len(result['coordinates']), constraint_coordinates=len(result['constraint_coordinates']),
        displayed_sections='Every original section, in source order',
        scope='Native record display only; no new mathematical calculation or biological linkage')
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--operation', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(render(args.gadget.resolve(strict=True), args.operation.resolve(strict=True), args.out.resolve())))
