# Reading source records in pages

An admitted scientific study exposes source cells, occurrences, fields and
operators through paged views. Each page preserves the selected source
identities and evidence spans. This runner reassembles the pages and checks
returned bytes against their original owners.

Use its source and producer options to select the publication:

```sh
python3 examples/bioinformatics/source-streaming/replay.py --help
```

Paging bounds the returned view. Admission and subsequent mathematical
operations have their own memory requirements. The [native data guide](../../native-data/README.md)
describes source ownership, while [dataset selection](../../../tools/rcdb/DATA.md)
covers bounded reads from working records.
