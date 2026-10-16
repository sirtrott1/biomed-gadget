# Finite homology diagrams

`gadget homology-diagram` accepts a finite ordered family of original native
publications and an explicitly oriented chain map between each consecutive pair.
Each vertex names every actual source grade carrier. Equal dimensions or equal
carrier numbers in independent publications do not establish identity.

Use the standalone `homology-diagram` and `replay-homology` commands for these
products. The `bioinformatics` catalogue describes assay workspace services.
The diagram schema declares its own publication family and complete maps.

The `rational_zigzag` product computes rational homology from the original cycle
and boundary quotients and retains a complete interval decomposition, including
generators that reconstruct every original arrow. The `integral_modules` product
retains integral cycle and boundary lattices, Smith actions, free coordinates,
torsion and induced homomorphisms. Integral output is a module diagram without
a general barcode classification. GADGET delegates both products to the
independent canonical SDK.

Interval endpoints are inclusive positions in the supplied finite vertex order.
They do not assign elapsed physical time to a vertex index.

This authored example assumes `four.rfcxs` and `two.rfcxs` have one vertex, one
loop and one attached face. Their face boundaries are respectively four and two
times the loop, with source grade carriers `10`, `11` and `12`. The declared map
is the identity on vertices and loops and multiplication by two on faces. It
induces the quotient homomorphism from Z/4 to Z/2 in grade one. These are exact
controls, with no biological measurement claim.

Save this declaration as `diagram.json` alongside those native sources:

```json
{
  "schema": "gadget.homology.diagram.v1",
  "product": "integral_modules",
  "grade": 1,
  "vertices": [
    {"source": "four.rfcxs", "carriers": [10, 11, 12]},
    {"source": "two.rfcxs", "carriers": [10, 11, 12]}
  ],
  "arrows": [{
    "forward": true,
    "maps": [
      {"rows": 1, "columns": 1, "entries": [[0, 0, "1"]]},
      {"rows": 1, "columns": 1, "entries": [[0, 0, "1"]]},
      {"rows": 1, "columns": 1, "entries": [[0, 0, "2"]]}
    ]
  }],
  "evidence": "Authored exact chain map on the declared integral lattices"
}
```

```sh
gadget homology-diagram --spec diagram.json --out diagram.rfcxs
gadget replay-homology --record diagram.rfcxs --out replay.rfcxs
```

Source paths resolve relative to the specification file, including when the
command runs from another directory. Replay uses complete archived native
sources and needs no original input path or model. C++ callers can use
`homology_diagram_record`, `run_homology_diagram`, `replay_homology_diagram` and
`homology_diagram_view` from `<gadget/homology_diagram.hpp>`.

`forward: false` declares a map from the next vertex to the current vertex.
Every grade map is required, including maps with zero rows or columns. Matrix
entries are distinct row, column and exact coefficient triples. The entire
rectangle is supported; omitted coefficients are supported zero. Original
source occurrences, missing support, inactive contributions, repeated rows,
metrics, cutoffs and raw byte owners remain in their original publications.
Rational boundaries use the source cell scales. Integral boundaries use their
declared lattice coefficients and require integral maps. No metric is used by
homology; supplied metric channels are retained unchanged.

Output contains native exact tensor fields, every original source archive, the
request and the canonical provider certificate. Replay recomputes from source
and checks the certificate. It never imports a rendered numerical result into
a new calculation. Malformed maps, failed boundary or chain identities, missing
carriers, nonintegral integer maps, exhausted caller policies and cancellation
refuse. The C++ API accepts a shared `CpuResources` and an optional canonical
execution lease. Exact elimination fill and complete output have their actual
costs; container policies do not measure process RSS or the GMP heap.

For further native calculation, the canonical
`decode_native_homology_diagram` function reads the result's `view.provider`
byte span and returns a validated source and carrier bound SDK diagram. Its
typed accessors expose the full chain, quotient or lattice maps and interval
generators with their actual publication owners. The numerical field view uses
separate declared coordinate axes for each displayed matrix.
