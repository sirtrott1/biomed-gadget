# Historical field execution

A temporal calculation can select different cutoffs from retained source
histories and apply a declared coupled observation law. Each endpoint keeps
its own carrier and field identity. The execution record binds the ordered
queries to the source versions they used.

This runner requires an explicit history fixture and writes to a new output
directory:

```sh
python3 examples/bioinformatics/temporal-execution/replay.py --help
```

[Versioned field histories](../versioned-history/README.md) describe the
endpoint maps. [Scientific query sequences](../../scientific/README.md) cover
execution of several requests through one admitted session.
