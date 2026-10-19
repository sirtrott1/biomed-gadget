"""Admit original source readings into bounded native working records."""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import subprocess
import tempfile

from gadget_data import DataWriter, Limits
from gadget_data._records import natural, graph


@contextmanager
def events(reader, profile, source, request, limits):
    """Only original input tokens cross this compatibility interface."""
    with tempfile.TemporaryDirectory(prefix="gadget-source-") as temporary:
        control = Path(temporary) / "request.json"
        control.write_text(json.dumps(request) + "\n")
        with tempfile.TemporaryFile() as errors:
            process = subprocess.Popen([str(reader), profile, str(source), str(control)],
                                       stdout=subprocess.PIPE, stderr=errors)
            def rows():
                while True:
                    line = process.stdout.readline(limits.decode_bytes + 1)
                    if not line:
                        status = process.wait()
                        if status:
                            errors.seek(0)
                            raise ValueError(errors.read(4096).decode(errors="replace"))
                        break
                    limits.admit(len(line) * 8, "decode_bytes")
                    if not line.endswith(b"\n"):
                        raise ValueError("unterminated original source event")
                    yield json.loads(line)
            try:
                yield rows()
                process.stdout.close()
                status = process.wait()
                if status:
                    errors.seek(0)
                    raise ValueError(errors.read(4096).decode(errors="replace"))
            finally:
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                process.wait()


def coefficient(token, maximum, *, integer=False):
    """Admit decimal allocation before the canonical exact converter runs."""
    from rexgraph.value import Absent, convert_number
    if token is None:
        return Absent
    if integer and token.isascii() and token.lstrip("+").isdecimal():
        if len(token) > maximum:
            raise MemoryError("source integer exceeds its coefficient admission")
        if len(token) <= 18:
            return int(token)
    mantissa, separator, exponent = token.lower().partition("e")
    digits = mantissa.lstrip("+-").replace(".", "")
    if not digits.strip("0"):
        return 0
    fractional = len(mantissa.partition(".")[2])
    # Saturate before int conversion, including very long exponent spellings.
    scale = 0
    for digit in exponent.lstrip("+-"):
        scale = scale * 10 + int(digit)
        if scale > maximum:
            raise MemoryError("source decimal exponent exceeds exact coefficient admission")
    if exponent.startswith("-"):
        scale = -scale
    demand = len(digits) + abs(scale - fractional)
    if demand > maximum:
        raise MemoryError("source exact coefficient exceeds its byte admission")
    value = convert_number(token, "decimal_exact")
    if integer:
        if value.denominator != 1:
            raise ValueError("declared integer count is fractional")
        return int(value)
    return value


def binding_context(row, header, dataset):
    columns, fields = header["source_columns"], row["fields"]
    positions = {name: index for index, name in enumerate(columns)}
    def token(name):
        return fields[positions[name]] if name in positions else ""
    monomer = token("BindingDB MonomerID").strip()
    reactant = token("BindingDB Reactant_set_id").strip()
    if not monomer.isascii() or not monomer.isdecimal() or not monomer.strip("0") or not reactant.isascii() or not reactant.isdecimal() or not reactant.strip("0"):
        raise ValueError("positive original monomer and reactant identifiers required")
    result = {"compound": "bindingdb:monomer:" + monomer.lstrip("0"),
        "reactant_set_id": reactant, "reported_construct": (dataset, row["index"]),
        "target_name": token("Target Name"),
        "organism": token("Target Source Organism According to Curator or DataSource")}
    for name in ("Ligand SMILES", "Ligand InChI", "Ligand InChI Key", "BindingDB Ligand Name", "pH", "Temp (C)", "Curation/DataSource", "Article DOI", "BindingDB Entry DOI", "PMID", "PubChem AID", "PubChem SID", "PubChem CID"):
        if name in positions:
            result[name] = token(name)
    start = header.get("chain_start")
    for position, value in enumerate(fields):
        name = columns[position] if position < len(columns) else columns[start + (position - start) % 12]
        if "Target Chain" in name or "Number of Protein Chains" in name or start is not None and position >= start:
            result["chain_field_" + str(position)] = (position, name, value)
    return result


def import_measurements(source, output, *, reader, profile, request, name, limits=None,
                        tile_rows=16, tile_columns=256):
    """Retain complete logical fields while admitting bounded physical records."""
    limits = limits or Limits()
    natural(tile_rows, "row block", positive=True)
    natural(tile_columns, "column block", positive=True)
    if profile not in ("gtex", "bindingdb"):
        raise ValueError("measurement source profile must be explicit")
    with DataWriter(output, name=name, limits=limits) as writer:
        writer.source("measurements", source)
        if profile == "gtex":
            for role in ("sample_attributes", "subject_phenotypes"):
                writer.source(role, request[role])
        declaration = {key: value for key, value in request.items()
                       if key not in ("sample_attributes", "subject_phenotypes")}
        profile_owner = graph("source profile", ((writer.identity, "profile"),))
        profile_owner.attach_metadata(1, 0, "source_profile", {"profile": profile, "declaration": declaration,
            "interpretation": "Original count or signal readings" if profile == "gtex" else
                "Reported affinity thresholds in nM with separate comparison codes; target constructs retain all source fields"})
        writer.native("profile", profile_owner)
        with events(reader, profile, source, request, limits) as stream:
            header = next(stream)
            if header.get("event") != "header" or min(header["rows"], header["columns"]) <= 0:
                raise ValueError("source requires complete positive declared axes")
            nrows, ncolumns = header["rows"], header["columns"]
            limits.admit(ncolumns * 256, "decode_bytes")
            row_block = min(tile_rows, max(1, limits.decode_bytes // (ncolumns * 512)))
            column_block = min(tile_columns, ncolumns)
            origin = lambda row: ({"source": "measurements", "selector":
                {"row": row["index"], "offset": row["offset"], "length": row["length"]}},)
            if profile == "gtex":
                for start in range(0, ncolumns, column_block):
                    samples = [next(stream) for _ in range(min(column_block, ncolumns - start))]
                    if any(row.get("event") != "sample" or row["index"] != start + i for i, row in enumerate(samples)):
                        raise ValueError("sample coordinates changed")
                    context = [{"donor": row["donor"], "tissue": row["tissue"],
                        **{"sample/" + name: value for name, value in zip(header["sample_columns"], row["sample_fields"], strict=True)},
                        **{"donor/" + name: value for name, value in zip(header["donor_columns"], row["donor_fields"], strict=True)}} for row in samples]
                    if start == 0:
                        schema_owner = graph("context declaration", ((writer.identity, "metadata columns"),))
                        schema_owner.attach_metadata(1, 0, "columns", {"sample": header["sample_columns"],
                            "donor": header["donor_columns"]})
                        metadata_schema = writer.native("metadata columns", schema_owner)
                    origins = tuple({"source": role, "selector": {"sample": row["index"],
                        "offset": row[prefix + "_offset"], "length": row[prefix + "_length"]}}
                        for row in samples for role, prefix in (("sample_attributes", "sample"), ("subject_phenotypes", "donor")))
                    writer.axis("observations", [row["id"] for row in samples], offset=start,
                        context=context, origins=origins,
                        index_labels=True, indexed_context=("donor", "tissue"))
                quantity, units = request["measurement"], "TPM" if request["measurement"] == "tpm" else "count"
            else:
                columns = tuple(header["source_columns"])
                schema_owner = graph("context declaration", ((writer.identity, "source columns"),))
                schema_owner.attach_metadata(1, 0, "columns", columns)
                metadata_schema = writer.native("source columns", schema_owner)
                for start in range(0, ncolumns, column_block):
                    endpoints = header["endpoints"][start:start + column_block]
                    writer.axis("observations", [row["endpoint"] for row in endpoints], offset=start,
                        context=[{"source_column": row["column"], "units": "nM"} for row in endpoints])
                quantity, units = "reported affinity threshold", "nM"
                writer.field("comparison", axes=("features", "observations"), shape=(nrows, ncolumns),
                    layout="exact", number_rule="exact", quantity="reported comparison code: 1 =, 2 <, 3 <=, 4 >, 5 >=, 6 approximate", units="dimensionless")
            writer.field("measurements", axes=("features", "observations"), shape=(nrows, ncolumns),
                layout="exact", number_rule="decimal_exact", quantity=quantity, units=units)
            for start in range(0, nrows, row_block):
                rows = [next(stream) for _ in range(min(row_block, nrows - start))]
                if any(row.get("event") != "row" or row["index"] != start + i or
                       len(row["values"]) != ncolumns for i, row in enumerate(rows)):
                    raise ValueError("original measurement coordinates changed")
                origins = tuple(item for row in rows for item in origin(row))
                if profile == "gtex":
                    labels = [row["id"] for row in rows]
                    context = [{"description": row["description"]} for row in rows]
                else:
                    # An occurrence remains distinct even when all source identifiers repeat.
                    labels = [row["index"] for row in rows]
                    context = [binding_context(row, header, writer.identity) for row in rows]
                writer.axis("features", labels, offset=start, context=context, origins=origins,
                            index_labels=True, indexed_context=("compound", "reactant_set_id", "target_name", "organism") if profile == "bindingdb" else ())
                for column in range(0, ncolumns, column_block):
                    stop = min(column + column_block, ncolumns)
                    values = [[coefficient(token, limits.decode_bytes // (16 * len(rows) * (stop - column)),
                                         integer=profile == "gtex" and request["measurement"] == "integer_counts")
                               for token in row["values"][column:stop]] for row in rows]
                    if profile == "gtex" and request["measurement"] == "integer_counts":
                        import numpy as np
                        if all(type(value) is int and 0 <= value <= 2**63 - 1 for row in values for value in row):
                            values = np.asarray(values, dtype="<i8")
                    writer.exact("measurements", values, offset=(start, column), origins=origins)
                    if profile == "bindingdb":
                        from rexgraph.value import Absent
                        writer.exact("comparison", [[Absent if value is None else value for value in row["comparisons"][column:stop]]
                            for row in rows], offset=(start, column), origins=origins)
            if next(stream, None) is not None:
                raise ValueError("source emitted undeclared rows")
        return writer.finish()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=("gtex", "bindingdb"))
    parser.add_argument("source", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--reader", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True); parser.add_argument("--name", required=True)
    parser.add_argument("--record-mib", type=int, default=16); parser.add_argument("--decode-mib", type=int, default=64)
    args = parser.parse_args()
    result = import_measurements(args.source, args.output, reader=args.reader, profile=args.profile,
        request=json.loads(args.request.read_text()), name=args.name,
        limits=Limits(record_bytes=args.record_mib * 1024**2, decode_bytes=args.decode_mib * 1024**2))
    print(json.dumps(result))


if __name__ == "__main__":
    main()
