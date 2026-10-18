"""Native dataset intake, exploration and selected materialization."""
import argparse
import json
from pathlib import Path

from . import DataStore, Limits, SourceResolver


def _coordinates(value):
    from .catalog import CoordinateSpec
    result = []
    for term in value.split(","):
        if ":" in term:
            parts = term.split(":")
            if len(parts) not in (2, 3) or any(not p for p in parts):
                raise argparse.ArgumentTypeError("coordinates require explicit start:stop[:step]")
            result.append(range(*map(int, parts)))
        else:
            i = int(term)
            result.append(range(i, i + 1))
    return CoordinateSpec(tuple(result))


def _display(value):
    try:
        return json.dumps(value)
    except TypeError:
        from rexgraph.value_codec import pack_value
        return json.dumps({"native_value": pack_value(value).hex(), "codec": "rexgraph.value"})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record-mib", type=int, default=16)
    parser.add_argument("--decode-mib", type=int, default=64)
    parser.add_argument("--selection-mib", type=int, default=64)
    parser.add_argument("--catalog-mib", type=int, default=4)
    commands = parser.add_subparsers(dest="command", required=True)
    dense = commands.add_parser("import-dense", help="Read an original NPY matrix or selected HDF5 dataset")
    dense.add_argument("source", type=Path); dense.add_argument("output", type=Path)
    dense.add_argument("--reader", choices=("npy", "hdf5"), required=True)
    dense.add_argument("--dataset"); dense.add_argument("--name", required=True)
    dense.add_argument("--field", default="measurements")
    dense.add_argument("--quantity", required=True); dense.add_argument("--units", required=True)
    dense.add_argument("--missing", choices=("nan",))
    dense.add_argument("--tile-rows", type=int); dense.add_argument("--tile-columns", type=int)
    coordinate = commands.add_parser("import-coordinate", help="Retain original Matrix Market occurrences")
    coordinate.add_argument("source", type=Path); coordinate.add_argument("output", type=Path)
    coordinate.add_argument("--name", required=True); coordinate.add_argument("--field", default="measurements")
    coordinate.add_argument("--quantity", required=True); coordinate.add_argument("--units", required=True)
    coordinate.add_argument("--implicit", choices=("missing", "observed_zero"), required=True)
    coordinate.add_argument("--block-entries", type=int, default=4096)
    records = commands.add_parser("import-records", help="Read through a native DatasetDeclaration")
    records.add_argument("source", type=Path); records.add_argument("output", type=Path)
    records.add_argument("--declaration", type=Path, required=True); records.add_argument("--name", required=True)
    records.add_argument("--grouping", choices=("contiguous",))
    table = commands.add_parser("import-table", help="Read declared numeric columns as exact native measurements")
    table.add_argument("source", type=Path); table.add_argument("output", type=Path)
    table.add_argument("--declaration", type=Path, required=True); table.add_argument("--name", required=True)
    table.add_argument("--measurement", action="append", required=True)
    table.add_argument("--label"); table.add_argument("--context", action="append", default=[])
    table.add_argument("--quantity", required=True); table.add_argument("--units", required=True)
    table.add_argument("--index-labels", action="store_true"); table.add_argument("--index-context", action="append", default=[])
    find = commands.add_parser("find-axis", help="Stream original axis labels and declared context matches")
    find.add_argument("store", type=Path); find.add_argument("--axis", required=True)
    find.add_argument("--label"); find.add_argument("--where", default="{}")
    index = commands.add_parser("build-index", help="Index literal labels on independent native axes")
    index.add_argument("output", type=Path)
    index.add_argument("--source", type=Path, action="append", required=True)
    index.add_argument("--axis", action="append", required=True)
    index.add_argument("--context-role", action="append", default=[])
    lookup = commands.add_parser("lookup", help="Stream matching original positions without joining sources")
    lookup.add_argument("index", type=Path); lookup.add_argument("literal")
    lookup.add_argument("--source", type=Path, action="append", required=True)
    bundle = commands.add_parser("export-bundle", help="Export native dataset dependencies and optional original files")
    bundle.add_argument("store", type=Path); bundle.add_argument("output", type=Path)
    bundle.add_argument("--dependency", action="append", type=Path, default=[])
    bundle.add_argument("--source-root", type=Path)
    for name, help_text in (("explore", "Inspect source, axis and field declarations without decoding field blocks"),
                            ("verify", "Check every native record and optional original sources"),
                            ("select", "Publish a selected native dataset with original coordinate maps"),
                            ("provenance", "Map selected coordinates to native records and original sources"),
                            ("export-native", "Export a selected field to a GADGET native publication")):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("store", type=Path)
        if name in ("select", "export-native", "provenance"):
            if name != "provenance":
                command.add_argument("output", type=Path)
            command.add_argument("--field", required=True)
            command.add_argument("--rows", type=_coordinates, required=True)
            command.add_argument("--columns", type=_coordinates, required=True)
        if name == "verify":
            command.add_argument("--source-root", type=Path)
        if name == "export-native":
            command.add_argument("--importer", type=Path, required=True)
            command.add_argument("--namespace", type=int, required=True)
    args = parser.parse_args(argv)
    limits = Limits(args.record_mib * 1024**2, args.decode_mib * 1024**2,
                    args.selection_mib * 1024**2, args.catalog_mib * 1024**2)
    if args.command == "import-dense":
        from .readers import import_dense
        result = import_dense(args.source, args.output, reader=args.reader, name=args.name, field=args.field,
            dataset=args.dataset, quantity=args.quantity, units=args.units, missing=args.missing,
            tile=None if args.tile_rows is None and args.tile_columns is None else
                (32 if args.tile_rows is None else args.tile_rows, 256 if args.tile_columns is None else args.tile_columns), limits=limits)
    elif args.command == "import-coordinate":
        from .readers import import_coordinate
        result = import_coordinate(args.source, args.output, name=args.name, field=args.field,
            quantity=args.quantity, units=args.units, implicit=args.implicit,
            block_entries=args.block_entries, limits=limits)
    elif args.command in ("import-records", "import-table"):
        from rexgraph.io.declaration import DatasetDeclaration
        from .readers import import_records
        declaration = DatasetDeclaration.from_bytes(args.declaration.read_bytes())
        if args.command == "import-records":
            result = import_records(args.source, args.output, name=args.name, limits=limits,
                                    declaration=declaration, grouping=args.grouping)
        else:
            from .readers import import_measurements
            result = import_measurements(args.source, args.output, name=args.name, limits=limits,
                declaration=declaration, measurements=args.measurement, label=args.label,
                context=args.context, quantity=args.quantity, units=args.units,
                index_labels=args.index_labels, indexed_context=args.index_context)
    elif args.command == "find-axis":
        from ._records import UNSET
        with DataStore(args.store, limits=limits) as store:
            for match in store.find_axis(args.axis, label=UNSET if args.label is None else json.loads(args.label),
                                         context=json.loads(args.where)):
                print(_display(match))
        return
    elif args.command in ("build-index", "lookup"):
        from contextlib import ExitStack
        from .index import AxisIndex, build_axis_index
        with ExitStack() as stack:
            sources = [stack.enter_context(DataStore(path, limits=limits)) for path in args.source]
            if args.command == "build-index":
                if len(sources) != len(args.axis):
                    raise ValueError("each source needs one corresponding axis declaration")
                result = build_axis_index(zip(sources, args.axis, strict=True), args.output,
                                          context_roles=args.context_role, limits=limits)
            else:
                index = stack.enter_context(AxisIndex(args.index, limits=limits))
                for match in index.find(json.loads(args.literal), datasets=sources):
                    print(_display(match))
                return
    elif args.command == "export-bundle":
        from contextlib import ExitStack
        from .export import export_bundle
        with ExitStack() as stack:
            stores = [stack.enter_context(DataStore(path, limits=limits)) for path in (args.store, *args.dependency)]
            resolvers = None if args.source_root is None else {
                store.identity: SourceResolver(args.source_root) for store in stores}
            result = export_bundle(stores[0], args.output, dependencies=stores[1:],
                resolvers=resolvers, include_originals=args.source_root is not None)
    else:
        with DataStore(args.store, limits=limits) as store:
            if args.command == "explore":
                result = store.catalog()
            elif args.command == "verify":
                count = store.verify()
                originals = None if args.source_root is None else {
                    key: str(path) for key, path in store.verify_sources(SourceResolver(args.source_root)).items()}
                result = {"native_records": count, "originals": originals, "reference": store.reference}
            elif args.command == "provenance":
                result = store.provenance(args.field, rows=args.rows, columns=args.columns)
            else:
                if args.command == "select":
                    result = store.materialize(args.field, args.output, rows=args.rows, columns=args.columns)
                else:
                    from .bridge import export_gadget
                    selected = store.select(args.field, rows=args.rows, columns=args.columns)
                    result = {"publication": str(export_gadget(selected, args.output, importer=args.importer,
                              namespace=args.namespace, max_bytes=limits.decode_bytes))}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
