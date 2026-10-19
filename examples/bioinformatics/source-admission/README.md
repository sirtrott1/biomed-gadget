# Table admission and source evidence

A strict table declaration assigns field meanings, coordinate identities and
observation support to original records. Admission stores shared source bytes
and the spans supporting each measurement or assertion. Repeated source
occurrences keep their own coordinates.

This runner checks admitted fields and category membership operations. It
creates an authored table by default; `--spec` selects an existing specification
beside its source files:

```sh
python3 examples/bioinformatics/source-admission/replay.py --help
```

For bounded dataset import and selection, see the
[Python data API](../../../tools/rcdb/DATA.md). Native annotation and byte-span
bindings are described in [study declarations](../../native-studies/README.md).
