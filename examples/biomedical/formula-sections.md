# Comparing original reported chemical formulas

Build the examples, then select a frozen database snapshot containing one ChEMBL
`molecule` response and one PubChem `compound` response. This adapter interprets
the original formula strings. It never reads a numerical result view as input.

```sh
python3 tools/prepare_chemical_formula_sections.py \
  --snapshot /absolute/path/to/snapshot \
  --chembl-response chembl-molecule --pubchem-response pubchem-compound \
  --binary /absolute/path/to/build/presets/cpu/gadget \
  --collection-tool /absolute/path/to/build/presets/cpu/gadget_biomedical_collection \
  --out '/absolute/path/to/new formula comparison'
```

The producer reported InChIKeys must agree. This is an explicit conditional
correspondence; it does not independently establish molecular structure,
stereochemistry or equivalence of other identifiers. ChEMBL PubChem Substance
identifiers remain distinct from PubChem Compound identifiers.

The adapter retains the complete original response and HTTP header bytes, the
snapshot declaration, the original formula strings and every parsed token with
its string offsets. SHA256 values and paths are checked before output creation.
Ordinary element symbols and positive integer counts are supported. Charges,
isotopes, brackets and mixtures refuse rather than being simplified.

Every element in either complete formula becomes an axis. An element absent
from a complete formula has a supported zero count. A missing formula refuses.
The two source studies retain distinct namespaces and carriers. Declared element
identity restrictions connect them through a common interface; full identity
cost forms define the discrepancy. The existing native section provider retains
the input, compatible and residual fields, exact defects, all compatible
directions and row dependencies. The cost is a declared source comparison, not
an estimate of experimental precision or a clinical conclusion.

`workspace.rfcxs` owns both original source studies and their provenance.
`result.rfcxs` owns the complete calculation and its workspace. `result.json`
is a researcher view. Reopen the workspace and query without the intake files:

```sh
gadget assay-query --workspace workspace.rfcxs --request query.json \
  --out-record replay.rfcxs
gadget record-view --input replay.rfcxs --member result
```

`gadget_biomedical_collection WORKSPACE REQUEST NEW_RECORD` also works with other
admitted assay studies. Its request follows
`schemas/biomedical-collection.schema.json`. Source paths resolve from the
process working directory. Every selected study retains the same declared
original byte owners; its existing field accessions and scientific digest stay
unchanged. Output creation refuses an existing destination.

# Joining retained source assertions

Biomedical adapter version 3 adds reported UniProt disease annotations and
ChEMBL `drug_indication` rows as source relationships. Import it explicitly:

```sh
gadget import-evidence --snapshot snapshot --out evidence --adapter-version 3
```

Versions 1 and 2 keep their prior interpretations; version 2 remains the default.
A disease comment without one accession remains a response local statement.
Reported phase values stay in their original assertions and do not become
binding measurements. Cross references remain assertions of their source.

The `role_join` operation uses the canonical native role join provider. Supply
the evidence digest, both assertion kinds and roles, activity selection and an
explicit `max_matches` execution budget. Kinds are `binding_observation`,
`identifier_cross_reference` and `source_relationship`; roles are
`compound_or_left`, `target_or_right`, `assay`, `document` and `response`.
The schema is `schemas/biomedical-query.schema.json`.

Each match retains the shared participant identity, both original assertion
identities, record identifiers and occurrence positions. The full source
publication retains the original records. A join budget refusal reports the
required result count; it does not truncate results. Matching uses explicitly
admitted identifiers within this publication. Separate identifiers or ambiguous
cross references do not become an inferred identity relation.
