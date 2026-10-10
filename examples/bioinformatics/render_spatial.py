#!/usr/bin/env python3
"""Render verified spatial section tables on their original coordinate carrier."""
import argparse
import csv
import hashlib
import html
from fractions import Fraction
from pathlib import Path


def rows(path):
    with path.open() as source:
        return list(csv.DictReader(source, delimiter='\t'))


def render(args):
    observations, features, summary = rows(args.observations), rows(args.features), rows(args.summary)
    symbols = {row['id']: row['symbol'] for row in features}
    selected = sorted({row['field'] for row in summary})
    expected = {(str(i + 1), field) for i in range(len(observations)) for field in selected}
    actual = {(row['anchor'], row['field']) for row in summary}
    if len(actual) != len(summary) or actual != expected:
        raise ValueError('Presentation requires one verified row per original spot and requested field')
    options = ''.join('<option value="' + html.escape(field, quote=True) + '">' +
                      html.escape(symbols[field.removeprefix('all.')]) + '</option>' for field in selected)
    statistics = '<table><thead><tr><th>Gene</th><th>Assigned UMI count</th><th>Nonzero spots</th><th>Charts with varying expression</th><th>Varying charts with zero anchor residual</th></tr></thead><tbody>'
    for field in selected:
        measured = [row for row in summary if row['field'] == field]
        total = sum((Fraction(row['count']) for row in measured), Fraction(0))
        nonzero = sum(Fraction(row['count']) != 0 for row in measured)
        varying = sum(Fraction(row['residual_quadrance']) != 0 for row in measured)
        zero_anchor = sum(Fraction(row['residual_quadrance']) != 0 and
                          Fraction(row['anchor_residual']) == 0 for row in measured)
        statistics += '<tr><td>' + html.escape(symbols[field.removeprefix('all.')]) + \
                      f'</td><td>{total}</td><td>{nonzero}</td><td>{varying}</td><td>{zero_anchor}</td></tr>'
    statistics += '</tbody></table>'
    axis_text = 'id\tx\ty\n' + ''.join(row['id'] + '\t' + row['x'] + '\t' + row['y'] + '\n'
                                           for row in observations)
    table_text = args.summary.read_text()
    if any('</script' in value.lower() for value in [axis_text, table_text]):
        raise ValueError('Unsupported presentation text')
    with args.publication.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    page = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Spatial expression and local compatibility</title>
<style>
body{margin:32px auto;padding:0 20px;max-width:1280px;font:16px/1.6 system-ui;color:#182b32;background:#fafaf7}
h1,h2{line-height:1.2}p{max-width:1000px}select{font:inherit;padding:6px}a{color:#165b75}
.maps{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}
figure{margin:0;background:white;border:1px solid #d6dfdf;padding:14px}canvas{width:100%;display:block}
figcaption{font-weight:600}.legend{font-size:13px;font-weight:400}#detail{min-height:3em;white-space:pre-wrap}
code{overflow-wrap:anywhere;font-size:13px}@media(max-width:700px){.maps{grid-template-columns:1fr}}
table{border-collapse:collapse;font-size:14px}th,td{text-align:right;padding:5px 14px;border-bottom:1px solid #d6dfdf}th:first-child,td:first-child{text-align:left}
</style>
<h1>Spatial expression and local compatibility</h1>
<p>The complete released breast tissue section contains @SPOTS@ tissue spots and @GENES@ genes.
These figures show each spot at its original image coordinates. Select a gene to compare its measured
assigned UMI count, its local compatible section and its residual. All selected fields use the same
declared scalar restrictions and identity metric.</p>
<label>Gene <select id="gene">@OPTIONS@</select></label>
<div class="maps">
<figure><figcaption>Measured count <span class="legend" id="count-scale"></span></figcaption><canvas id="count"></canvas></figure>
<figure><figcaption>Anchor residual <span class="legend" id="anchor_residual-scale"></span></figcaption><canvas id="anchor_residual"></canvas></figure>
<figure><figcaption>Local residual quadrance <span class="legend" id="residual_quadrance-scale"></span></figcaption><canvas id="residual_quadrance"></canvas></figure>
<figure><figcaption>Anchor Green response <span class="legend" id="anchor_green-scale"></span></figcaption><canvas id="anchor_green"></canvas></figure>
</div><p id="detail">Move the pointer over a spot to inspect its original barcode and exact values.</p>
<h2>Reading the figures</h2>
@STATISTICS@
<p>The table counts each original anchor once. A chart has varying expression when its exact
residual quadrance is nonzero. An all zero chart is compatible with the scalar identity restrictions;
that compatibility records equal observed values. The last column identifies charts whose anchor
equals the compatible mean while neighboring measurements differ. Their anchor residual and anchor
Hodge value vanish, but the full residual field has nonzero quadrance.</p>
<p>Each local chart contains one anchor, its incident physical neighborhood relations and their
complete boundaries. The radius is 75,413 squared pixels in the published image frame. No higher
cell is inferred. Each two occurrence relation has oriented coefficients minus one and plus one.</p>
<p>For one gene, an identity restriction declares equal values to be compatible across a relation.
The exact compatible section on a connected chart is its mean. The anchor residual is the measured
anchor count minus that mean. Blue indicates a negative residual and red a positive residual. Local
residual quadrance sums the squared residuals at every chart observation. Overlapping charts share
the original observations; their separately projected sections can disagree on that overlap.</p>
<p>The local Green response solves (I + L)y = x with the chart's native grade zero Hodge operator.
It is a response to that declared operator. It remains separate from the original count field.
The independent verifier checks every projected value, residual, coboundary defect, Hodge action,
Green equation and constant section sum using exact rational arithmetic.</p>
<p>Each map uses a linear color scale for the selected gene. Colors use floating point display
coordinates; the underlying result and pointer values retain exact fractions. Quadrance has squared
UMI units and depends on chart size. These maps describe one measured section and the stated
compatibility conditions. Tissue annotations and treatment response require separately supplied evidence.</p>
<h2>Reproduction and source</h2>
<p>Use <code>examples/bioinformatics/spatial-expression.md</code> and its preparation script to
download the pinned original archives, import every gene and tissue spot, construct native fields
and charts, run the native scientific sequence and verify it. This page is a presentation of the
verified table, with every selected spot retained.</p>
<p>Data source: <a href="https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0">10x Genomics human breast cancer block A section 1</a>.
The publisher supplies this dataset under CC BY 4.0.</p>
<p>Native result SHA256: <code>@DIGEST@</code>.</p>
<script id="axis" type="text/plain">@AXIS@</script>
<script id="table" type="text/plain">@TABLE@</script>
<script>
const parse=id=>{const lines=document.getElementById(id).textContent.trim().split('\\n');
const keys=lines.shift().split('\\t');return lines.map(line=>Object.fromEntries(line.split('\\t').map((v,i)=>[keys[i],v])))};
const fraction=s=>{const p=s.split('/').map(Number);return p.length===2?p[0]/p[1]:p[0]};
const axis=parse('axis'),table=parse('table'),fields=['count','anchor_residual','residual_quadrance','anchor_green'];
const xs=axis.map(p=>fraction(p.x)),ys=axis.map(p=>fraction(p.y));
const xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
const scale=460/Math.max(xmax-xmin,ymax-ymin,1);
const positions=axis.map((p,i)=>[20+(xs[i]-xmin)*scale,20+(ys[i]-ymin)*scale]);
let current=[];
function draw(){const gene=document.getElementById('gene').value;
const byAnchor=new Map(table.filter(row=>row.field===gene).map(row=>[row.anchor,row]));
current=axis.map((p,i)=>byAnchor.get(String(i+1)));
for(const field of fields){const canvas=document.getElementById(field);canvas.width=500;canvas.height=500;
const ctx=canvas.getContext('2d');ctx.fillStyle='#ffffff';ctx.fillRect(0,0,500,500);
const values=current.map(row=>fraction(row[field]));const signed=field==='anchor_residual';
const max=Math.max(...values.map(v=>signed?Math.abs(v):v),0);
document.getElementById(field+'-scale').textContent=(signed?'−':'0 to ')+max.toPrecision(4)+(signed?' to '+max.toPrecision(4):'');
values.forEach((value,i)=>{const t=max?Math.min(1,Math.abs(value)/max):0;
ctx.fillStyle=signed?(value<0?`rgb(${Math.round(245-220*t)},${Math.round(245-135*t)},${Math.round(245-60*t)})`:
`rgb(${Math.round(245-45*t)},${Math.round(245-210*t)},${Math.round(245-210*t)})`):
`rgb(${Math.round(240-220*t)},${Math.round(246-110*t)},${Math.round(247-100*t)})`;
ctx.beginPath();ctx.arc(...positions[i],2.5,0,2*Math.PI);ctx.fill()});}}
for(const field of fields){document.getElementById(field).onmousemove=event=>{const box=event.target.getBoundingClientRect();
const x=(event.clientX-box.left)*500/box.width,y=(event.clientY-box.top)*500/box.height;
let best=Infinity,index=0;positions.forEach((p,i)=>{const d=(p[0]-x)**2+(p[1]-y)**2;if(d<best){best=d;index=i}});
if(best>64)return;const row=current[index];document.getElementById('detail').textContent=
axis[index].id+'  '+row.context+'  '+row.local_cells+' local observations\\n'+
'Count '+row.count+'   Compatible mean '+row.local_mean+'   Residual '+row.anchor_residual+'\\n'+
'Residual quadrance '+row.residual_quadrance+'   Hodge '+row.anchor_hodge+'   Green '+row.anchor_green;};}
document.getElementById('gene').onchange=draw;draw();
</script></html>'''
    for token, value in {'SPOTS': str(len(observations)), 'GENES': str(len(features)), 'OPTIONS': options,
                         'DIGEST': digest, 'AXIS': axis_text, 'TABLE': table_text,
                         'STATISTICS': statistics}.items():
        page = page.replace('@' + token + '@', value)
    args.out.write_text(page)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('observations', 'features', 'summary', 'publication', 'out'):
        parser.add_argument('--' + name, required=True, type=Path)
    render(parser.parse_args())
