# Sequencing and assay recipes

GADGET separates physical records, assignment choices, measured fields and
relational structure. These recipes prepare fields for the same native analysis
interface used by bulk expression, cell-level assays and longitudinal studies.

## Build with genomic codecs

SAM, BAM, CRAM, VCF and BCF use HTSlib 1.20 or later. Compression is decoded by
HTSlib; GADGET retains the original compressed bytes and their hash. Enable it
with `-DGADGET_ENABLE_HTSLIB=ON`. Core builds without this optional dependency;
their capability catalog omits these formats and rejects their admission.
H5AD is enabled separately with `-DGADGET_ENABLE_HDF5=ON`.

A genomic resource declares its assembly, hash of the reference FASTA, contig
dictionary and zero-based half-open coordinate convention. Its source hash binds
the complete input file. CRAM also requires `reference_source`, a path/hash record
for that exact FASTA. GADGET supplies a private reference image and index to the
decoder. Resource admission retains that image for portable reimport.

Decoded binary records are identified by file hash and ordered record index.
Their decoded SAM/VCF representation accompanies the original file; record
indices are not invented byte offsets. Binary end markers, CIGAR bounds, query
lengths, read-group dictionaries, GT allele indices and reference lengths are
checked. Record, genotype and CIGAR-piece counts follow the complete decoded
source, with checked size arithmetic. The retained file, workspace manifest and
certificate use their actual byte/container representations. There is no
application record-count or fixed MiB ceiling on these adapters.

## Counting and splice evidence

`alignment_evidence` returns aligned blocks, CIGAR `N` introns, junction support
and the exact selected/excluded record indices under the declared filter.
Mapping quality 255 means unavailable; the filter explicitly includes or excludes
it. Junction counts count alignment records. Multiple mappings remain separate
evidence rather than being silently deduplicated.

`feature_counts` groups retained mappings by read group, query name and mate
identity. This is a **sequenced segment**, so paired ends are separate units.
Annotation grouping uses an explicit GTF/GFF3 attribute, such as `gene_id`,
`transcript_id` or `Parent`; BED uses its record ID. Repeated intervals belonging
to one feature are unioned, while original records remain identifiable. Shared
GFF3 parents remain shared. Strand policies are `ignore`, `same` or `opposite`
relative to the alignment strand; these are not automatic library-protocol calls.

```json
{
  "schema": "gadget.assays.query.v1",
  "operation": "feature_counts",
  "resource": "alignments",
  "annotation": "genes",
  "feature_type": "exon",
  "group_attribute": "gene_id",
  "unit": "sequenced_segment",
  "overlap": "any_aligned_base",
  "strand_policy": "ignore",
  "filter": {
    "exclude_flags": 1024,
    "minimum_mapping_quality": 0,
    "unknown_mapping_quality": "exclude"
  },
  "candidate_search": {
    "method": "The recorded upstream aligner and parameters",
    "completeness": "Assignments are conditional on the retained mappings",
    "reference_sha256": "sha256:<reference FASTA digest>"
  }
}
```

Only `M`, `=` and `X` positions enter aligned-base overlap. Deletions and skipped
introns do not. The alternative rule `all_aligned_bases` requires each reported
aligned base to belong to the feature's interval union. Every overlap includes
its base count and source-record witness. A retained alternative outside every
selected feature adds an unassigned candidate, even when another mapping overlaps
a gene. Filters and upstream candidate search determine the resulting universe.

RexGraph retains candidate classes as primary multiway relations and computes
sharp marginal assignment bounds through exact count maps. A segment contributes
to a feature's lower bound when every candidate belongs to that feature; it
contributes to its upper bound when any candidate does. The output is a bound
field, with assignment uncertainty preserved. It is not an estimated abundance.
Paired-template, UMI and barcode deduplication require additional declared models
and are not inferred from a query name or optional tag.

## Coverage, variants, quality and sequence pairs

| Method | Result |
|---|---|
| `read_quality` | Read lengths, stored Phred+33 score counts by cycle, caller thresholds and IUPAC GC bounds |
| `annotation_features` | Attribute-grouped annotation interval unions, physical frames and original record indices |
| `chromatin_contrast` | Common physical intersection segments, complete scaled coverage fields, matched differences, integrals, quadrances and threshold components |
| `variant_genotypes` | Complete called allele-copy fields, ordered allele dictionaries, ploidy, phasing separators and original variant records |
| `exact_sequence_alignment` | A complete edit grid, ordered optimal path, all optimal-step support, and native prefix/suffix primal–dual certificates |

Chromatin samples declare positive scales with their producing evidence, named
treated/control pairs, strand handling and thresholds. The caller explicitly
declares measured coverage zeros within the joint interval union. All memberships
are retained. Thus one can compare treatment/control structure across replicates
without choosing a statistical peak model. Regions outside the retained union
are outside this calculation. Threshold components describe observed coverage;
they are not calibrated peak significance.

Variant `filter_policy: pass` selects only `PASS`; a dot means filters were not
applied. `all` retains every input call. Allele-copy zero requires a complete GT;
partial and missing GT remain absent, including mixed-ploidy calls. REF/ALT and
upstream likelihood/filter fields remain in the source representation. Variant
calling, haplotype reconstruction and pathogenicity annotation are separate tasks.

The sequence profiles use DNA/RNA IUPAC bases, including U; reverse-complement
search requires DNA tokens. Pair alignment uses FASTA records and the complete
`(query_length + 1) × (reference_length + 1)` global edit grid, with nonnegative
rational match, mismatch and linear gap costs. `literal` compares tokens;
`possible` accepts a nonempty IUPAC intersection; the pair method's `definite`
rule requires equal singleton sets. The result names this compatibility rule.

RexGraph's native routes produce prefix and suffix costs. An attained minimum
route also attains the minimum on each of its subroutes, so GADGET retains each
native route once and references its positions for previously unfilled source
coordinates. Native primal–dual checks verify both complete potential fields,
and native boundary actions verify tightness along every referenced route.
`route_witnesses`, `prefix_suffix_witnesses` and `native_route_calls` expose this
reuse. The chosen ordered route, its realized field and every step on any
minimum route remain separate outputs. All source relations and the complete
cost, chosen-path and optimal-support fields retain supported zeros.

Sequence lengths and grid size have no separate policy ceiling. Checked native
cell, occurrence and integer domains govern representation. The complete grid,
route histories and result stay in memory; this method retains their actual
construction and routing cost. It implements global linear-gap alignment.
[The portable example](sequence-witnesses/README.md) provides both inputs and a
pipeline request.

## Cells, spatial structure and assay integration

`aggregate_assay` takes explicit feature and observation maps. Each group names
original identities and rational coefficients. This serves specimen-level sums
of single-cell counts, supplied gene-set totals, declared region aggregates and
other linear assay questions. Arithmetic means require explicit coefficients.
Any unsupported input with a nonzero coefficient makes the corresponding output
absent. Equal shapes do not establish assay correspondence, and grouping does
not automatically transfer the source topology.

Selected panels, rational weighted maps and count sums read the already admitted
native field pins. Each stored packet is checked against its disjoint range in the
original axis dictionary, original identities and actual publication. RexGraph's exact
matrix actions perform the maps; its tensor window handles complete unit sums
along the stored axis. Unsupported coordinates retain their dependency semantics.
Execution receipts retain every packet's source digest, namespace, cutoff and
coordinate range without parsing the raw matrix again. Complete unit sums use
native window actions within packets and an explicit exact matrix action across
their results. Dictionary-based designs compile the same declared maps into this
route. Storage ranges address measurement coordinates; they do not glue source
topology or define biological groups.

Imports retain the complete carrier by default within the native source index
domain. An explicit `partition_size` requests disjoint storage ranges; the
numerical results keep the same original-coordinate interpretation.

`spatial_relations` reads retained observation coordinates with their explicit
`coordinate_frame` and `coordinate_units` columns. A positive diagonal native
metric and caller radius squared define physical neighbors without a fitted
embedding. Each pair becomes an oriented grade-one relation; coordinate order
sets the sign convention. Neighborhoods describe the supplied physical geometry.
Higher attachments and biological interactions need separate source declarations.
The existing `field_transfer` method compares assay fields under an explicit
measurement map, retaining its rank, kernel, residual and coordinate identities.

## Continue with native relational analysis

`derive-assays` materializes feature-count bounds, genotype-copy fields, chromatin
contrasts, declared aggregates, physical spatial relations or complete edit grids into a new canonical
assay workspace. Its declaration names the generating recipe, new namespace,
study design and independent GET channel. Source files, the complete producing
certificate, map declarations and original measurement metadata are retained.
Spatial derivation appends physical relations to authored structure, preserving
existing higher attachments and their fields; duplicate identities reject.
Edit grids retain prefix/suffix potential fields at grade zero and independent
cost, chosen-path and optimal-support fields at grade one. Hodge and Green
therefore act on the actual retained edit relations. Alternative edit histories
are kept; sharing endpoints does not automatically create a higher face.

`replay-derivation` reconstructs the used source inputs from a derived workspace,
runs its original generating method again, and compares the complete mathematical
result. The original source workspace is unnecessary. Changed workspace and
build receipts are compared separately from numerical payloads; the native
dependency fingerprint remains checked.

```sh
gadget derive-assays --workspace source/workspace.json \
  --request derivation.json --out derived
gadget assay-query --workspace derived/workspace.rfcxs --request native-query.json
gadget replay-derivation --workspace derived/workspace.rfcxs --assay experiment \
  --out regenerated
gadget workflow --workspace derived/workspace.rfcxs --assay experiment \
  --request workflow.json --models roles.json
```

The resulting studies use the existing controls, paired fields, Hodge, Green,
local sections, field transport and temporal recipes. Local llama.cpp models can
request or observe supported native analyses through the same bounded workflow
interface; model prose remains a separate observation. CPU reservations and
native work graphs retain ownership of execution. Packed CPU/GPU capabilities
apply to their existing native packed operations, not to rational metric solves.

The independent suite test checks per-base chromatin fields, strand and CIGAR
semantics, ambiguous assignments, phased and missing genotypes, known spatial
distances, support propagation, archive relocation and exhaustive tiny edit-path
oracles. These authored inputs qualify the contracts. Biological experiments
must supply their own designs and sources; the airway expression demonstration
provides a separately retained real assay example.

`reference_sequence_matches` returns every accepted ungapped placement in the
retained FASTA contigs, on the declared strands, under an integer mismatch bound.
It also searches every placement across one explicitly supplied junction with
the declared minimum overhang. FASTA and FASTQ query records retain their
identities and source spans, including unmatched records. Different junction
assertions retain their separate identities even at equal physical coordinates.

For a literal ACGT query of length `m` with mismatch bound `k < m`, at least one
of `k + 1` disjoint seeds is unchanged in every acceptable placement. The existing
injective two-bit index retains every requested seed hit. A seed is shortened
when needed to fit its intrinsic 64-bit representation; this preserves the
necessary condition. Full query-length native mismatch actions verify every
candidate. IUPAC queries, nonliteral rules and `k >= m` use the full contiguous
placement universe. Every declared-junction placement is checked in full.

`literal` means equal tokens and `possible` means a nonempty IUPAC intersection.
This method's established `definite` rule requires the query's IUPAC set to be
contained in the reference set; `compatibility_rule` records that directed rule.
It differs from the singleton equality used by pair alignment. Result ordering
follows query, strand, reference and start, then junction declaration and split.
`candidate_universe` records all admitted placements; `tested_candidates` and
`base_comparisons` record actual native verification work. These exact counts
use decimal strings when they exceed signed JSON integer representation.

Query, reference, mismatch, overhang, comparison, junction and accepted-result
counts have no separate policy ceiling. Source, index and complete outputs
remain resident. Repetitive seeds and IUPAC/junction fallback incur their full
work. The declared universe contains contiguous and single-junction placements;
indels and undeclared junctions belong to other methods.

`sequence_assembly_graph` constructs the complete observed literal k-mer
transition complex. Prefix and suffix identities are shared only when the tokens
are equal. Each retained k-mer relation has two independently keyed signed
occurrences; a self-loop keeps both even though its boundary cancels. The
occurrence-count field is computed by the native count map. Every read's walk,
sequence offsets, original record span and reconstructed sequence are retained.
All branches remain visible. An optional declared pair of endpoints requests
one native minimum-hop walk. Its sequence is a graph walk, rather than an inferred
biological consensus. Materialization admits the graph, independent grade-one
counts and grade-zero outgoing-vocabulary counts for native Hodge, Green and
sheaf analyses. `k >= 2`, and every read must be at least k bases long. Source
read, k-mer, node, occurrence and k counts have no separate policy ceiling.
The native exact count map contains one membership per original occurrence;
`source_occurrences` and `count_operator_nonzeros` expose that complete support.
Cycle dimension is computed by native exact boundary rank. The complete source,
walks and outputs stay in memory; native rank retains its actual elimination
cost. Native cell, identity, occurrence and codec domains govern representation.

`constrained_allocations` admits finite lower/upper bounds, exact rational linear
constraints (`eq`, `le`, `ge`) and an explicit minimizing or maximizing objective.
RexGraph solves the equalities and fixed bounds as `x = x0 + K z`, with a complete
native kernel frame. Free parameters retain the corresponding original variable
coordinates. Reduced inequality and bound supports determine independent
parameter groups; a native contrast-kernel frame certifies their partition.
Each original variable and inequality belongs to at most one group. The original
rows, coefficients, bounds, objective and names remain available for verification.
The reported `affine_dimension` is the dimension of the equality and fixed-bound
solution space; inequalities may restrict it further. Arithmetic is exact
rational arithmetic. Linear feasibility and optimization require no metric.
Authored constraints are native grade-one relations; every dense source term,
including a zero coefficient, remains an existent occurrence. Arity does not
create a higher-grade attachment.

`bounded_integer` enumerates each independent free-coordinate box and checks
integrality of reconstructed coordinates. `bounded_rational` enumerates distinct
active hyperplanes within each group, using native exact normalization, rank,
solve and residual checks. Scaled duplicate rows retain their original evidence
but share active search work. Inconsistent equalities receive a native left-null
contradiction witness; constant bound and integrality contradictions are also
reported. Degenerate feasible polytopes remain supported.

Bounds describe the complete feasible and optimal sets in original coordinates.
Explicit output retains every optimal integer point or rational vertex. Set
`witness_representation: factorized` to request the complete affine map and every
optimal tuple per group, with their Cartesian product count. Rational vertex
products span the complete optimal face by convex hull. Set
`factor_certificate: true` to retain these factors alongside explicit witnesses. Counts exceeding the
signed 64-bit JSON integer representation use decimal strings. Materialized
assays use the same sharp bound endpoints in either representation.

There are no application ceilings on variable/constraint counts, candidates,
feasible points or optimal witnesses. Native cell indices and source/container
representations govern admission. Genuinely coupled enumeration and complete
explicit output incur their actual costs. Every variable still requires declared
finite bounds. [The coupled-allocation example](coupled-allocation-factors.json)
shows declared class totals, a shared capacity and an independent ambiguity fiber.

`relational_persistence` acts on one admitted assay partition and its complete
source complex. Supply increasing rational `levels`, a birth-level index for
every `{grade,index}` coordinate, the filtration meaning and evidence. Every
level must be closed under all retained boundary occurrences, including inactive
ones. The method preserves authored higher relations, verifies all chain
compositions and native inclusion maps, and records cycle/boundary frames and
every induced homology image rank. Rational interval multiplicities must
reconstruct the full rank invariant. Endpoints are birth-inclusive and
death-exclusive. Surviving intervals are censored at the final supplied level;
the input does not establish their future lifetime. All source cells and supplied
levels are retained. Integral torsion requires a separately declared lattice; this profile
uses rational chain homology.

Run the independent examples and inspect their requests and certificates:

```sh
python3 tests/finite_services_cli.py build/gadget/gadget results/finite-services
gadget assay-query --workspace results/finite-services/workspace/workspace.json \
  --request results/finite-services/persistence-request.json
gadget replay-derivation \
  --workspace results/finite-services/derived-assembly/workspace.json \
  --assay derived-assembly --out results/regenerated-assembly
```

The native query `chain_certificate` takes a context, grade and either a source
field or a complete vector. `coefficient_domain: "rational"` uses the declared
rational coordinates. `"retained_lattice"` accepts integer vectors in the source's
coordinated lattice basis. A source field always uses rational coordinates.
The method retains every incidence, including inactive assertions, and records
its signed contribution and lower-stalk sum. It reports cycle and boundary
membership, integral free rank and torsion, higher-chain filling witnesses,
gradient/curl/harmonic fields, their metric quadrances, and complete field
reconstruction. Source size follows the complete native complex and checked
coordinate/container representations.

```json
{"operation":"chain_certificate","context":"all","grade":1,
 "coefficient_domain":"rational","field":"all.kmer_occurrence_count"}
```

`realized_path_comparison` accepts named paths, each with ordered
`steps: [{"index": 0, "coefficient": "1"}]`, on one positive source grade.
It returns full realization and boundary matrices, both inherited metric forms,
their ranks, and the Hodge sectors of every pair difference. Repeated steps remain
in the request and receipt. Matching boundaries do not identify equal fields.
Declared realizations need not be simple or forward directed paths.

`relational_persistence` also accepts `lattice_snapshots: true`. Every stage then
reports free rank and torsion from its retained integral lattice. Its barcode and
inclusion ranks remain over Q; stage homology groups alone do not classify an
integral persistence module.

For local structured readers, a role can declare `task: "claim_reader"`. The
model returns `gadget.model.claims.v1` with the source digest, the supplied
certificate hash, and typed paths and values. A workflow `check_claims` node with
`claims_from` verifies that binding and compares each claim with the complete
certificate. Its contradiction count can select a reinspection branch. Model
output remains prospective until checked. The structural certificate also provides
a deterministic explanation derived from the checked incidence and field facts.

## Indexed sequence search and read reconstruction

`indexed_sequence_alignment` searches every nonempty reference substring within
a declared unit edit budget. It retains one exact minimizing edit history for
each admitted placement, its complete declared route carrier, signed boundary
check and cost check. The index uses injective two-bit DNA keys for the requested
seeds. With at most `e` edits, one of `e+1` disjoint seeds must survive unchanged;
its reference position locates the substring start within `e` bases. Every
admitted path lies inside `|query_prefix-reference_prefix| <= e`. These bounds
permit complete pruning before `rfc::path::dijkstra_exact` solves each endpoint.

```json
{
  "schema": "gadget.assays.query.v1",
  "operation": "indexed_sequence_alignment",
  "query_resource": "reads",
  "reference_resource": "reference",
  "strands": ["+", "-"],
  "scoring": "unit_edit_cost",
  "maximum_edits": 2,
  "evidence": "The admitted FASTA reference and retained read records"
}
```

The edit budget must be strictly less than every query length, which ensures
that the seed partition has `e+1` nonempty parts. The index retains every hit
for each requested seed. Seed keys encode up to 32 bases injectively in a
64-bit word; longer reads use seeds within that representation. Candidate
starts, band carriers, native route calls and placements follow the complete
source and declared budget. Native cell coordinates and checked machine sizes
bound the representation. Every returned placement includes one minimizing
ordered edit history, its realized field, boundary and cost checks.

Reference ambiguity tokens cannot match a literal query base at zero cost.
Queries with ambiguity require an explicitly different question, such as the
existing `reference_sequence_matches` compatibility rules. Source intake and
complete certificates remain in memory; route work and output size grow with
the declared search universe. This implementation does not yet provide a
streaming chromosome-scale alignment service.

`sequence_splice_discovery` requires no annotation or supplied junction list.
It indexes every query overhang using an exact seed within the native 64-bit
key domain, verifies the full overhang and retains every matching exon extension.
It declares transitions satisfying query order, reference order, intron lengths
and transcript-oriented boundary motifs. Native routes minimize junction count;
the result retains all admitted complete exon chains, including unspliced
placements. Supply `query_resource`, `reference_resource`, `strands`,
`minimum_exon_bases`, `minimum_intron_bases`, `maximum_intron_bases`,
`maximum_junctions`, `motifs` such as `["GT-AG", "GC-AG", "AT-AC"]`, and
`evidence`. Coordinates are zero-based half-open in reference orientation. On
the negative strand, the reported motif is reverse-complemented into transcript
orientation. Exons match literally, with zero exon edits. The caller's junction
bound defines the universe; there is no additional four-junction policy.
Query prefix and exon coordinates index the query in reference orientation;
on `-` they index its reverse complement. Original read intervals are obtained
as `[m-end, m-begin)`. Unresolved reference tokens cannot form a literal DNA
boundary motif, and are skipped without removing valid chains elsewhere.
Reference-start interval lookup excludes only transitions outside the declared
intron range. An explicit traversal stack enumerates every admitted chain in
source order. Query, exon, transition, comparison, junction and match-output
counts have no separate policy ceiling. Complete carriers and output retain
their actual costs; repeated source patterns can produce many valid chains.

`read_overlap_reconstruction` minimizes the length of a suffix-prefix cover of
every source read, each used once in its retained orientation. Its request
contains `resource`, positive `minimum_overlap`,
`objective: "minimum_suffix_prefix_read_cover"`, and `evidence`. A native
relational state records the used source-read identities and last read. Every
transition costs the number of added bases. For fixed source and target states,
a shorter admissible overlap is strictly dominated by the longest one; read
order is still optimized globally. All original reads remain distinct, including
equal spellings. Native directed reachability identifies impossible starts: a
complete cover starting at read a must be able to reach every other read in the
full overlap graph. `overlap_carrier` and `start_reachability` retain the carrier,
attained paths and unreachable-target witnesses for those decisions.

The derived search carrier declares every reachable prefix state from the starts
that survive this necessary condition. This avoids constructing unreachable
members of the full subset universe. `source_reads`, `source_read_mask` and
`last_read` retain the exact original coordinates of each state; source-sized
native packed bits support masks beyond one machine word. This carrier is a
restriction of the formal subset-state search universe, and its retained steps
include every feasible complete cover. It is not the earlier enumeration of
all formal subsets, and its native carrier digest and arena addresses can change.
Original read variables, overlap alternatives and minimum layouts are preserved.

Native minimum-route histories certify complete prefix costs and attained suffix
costs. Dead-end suffix costs are absent. The native suffix dual is checked on the
boundary-closed viable vertex/edge restriction, retaining the original coordinate
maps. The complete prefix dual is checked on the full derived carrier. During
heap-based layout traversal, native suffix witnesses exclude only prefixes that
cannot attain the global minimum. Source, overlap, layout, spelling and visit
counts have no separate policy ceiling. Dense overlap comparisons, branching
state spaces and complete requested outputs retain their actual costs.

The result retains every minimum layout and every distinct minimum sequence.
Each layout supplies original read placements and native categorical section
gluing with complete overlap coverage; restricting the glued sequence recovers
every source read. Equal read strings retain separate source identities.
`identifiable_sequence` means that all minimum layouts spell one sequence under
these declared constraints. Interior containment, orientation changes and error
correction require additional declarations. The shared-field Gram compares
representative layouts on the read-cover transition carrier with identity
metric. GADGET checks each representative's native endpoint field once; RexGraph
computes the complete inherited Gram. Equal endpoint fields certify every pair's
closed difference by linearity, without repeating the boundary calculation for
each Gram entry. This is a carrier of reconstruction states, rather than a physical
genome coordinate metric.

[The portable composite example](composite-sequences/README.md) contains the
source inputs and a pipeline selecting these three methods.

`query_batch` joins the declared independent workspace queries through
`rfc::runtime::WorkGraph`, with positive `concurrent_jobs` admitted against the
observed native CPU capacity and
`queries: [{"id": "search", "query": {...}}]`. Each child is a complete
`gadget.assays.query.v1` request. Source-query jobs receive one native CPU worker;
a `native_query` receives its declared worker demand. The result follows declared
job order after all owned jobs join. Recursive batches reject. Jobs retain only
source-order receipts in the graph, with complete certificates staged under one
private owner. The joined response reads every certificate in declaration order.
The final complete JSON response still occupies memory and uses json-c's actual
representation; no fixed certificate MiB ceiling is imposed here.

These operations use the same command as other assay queries:

```sh
gadget assay-query --workspace workspace/workspace.json --request search.json
```

The general `expression_allocations` recipe joins raw transcript matching, exon-compatible alignment evidence or unweighted equivalence classes to exact expression bounds and sample contrasts. See [the RNA evidence example](rna-evidence/README.md) for its full contract, coupled constraints and source-only replay.
