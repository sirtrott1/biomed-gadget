"""RCQL accessions for native dataset field records."""
import numpy as np

from ._records import declaration
from .selection import block_values, coordinates


def field(source, columns):
    """Bind selected measurement coordinates to the original native row carrier."""
    from rexgraph.chain_map import CoordinateComplex
    from rexgraph.tensor_field import FieldSource, TensorField
    from rexgraph.type_accession import CoordinateSpace
    from rexgraph.value import convert_number
    from rcql.execution_trace import current_binding, record_method
    value = declaration(source)
    if value["kind"] != "field" or value["layout"] not in ("dense", "exact"):
        raise ValueError("native dense or exact field record required")
    selected = coordinates(columns, value["shape"][1])
    raw, present = block_values(source, tuple(range(value["shape"][0])), selected)
    if not present.all():
        raise ValueError("unobserved selected coordinates require a declared missing data operation")
    coefficients = np.asarray([convert_number(v, value["number_rule"]) for v in raw.flat], dtype=object).reshape(raw.shape)
    binding = current_binding()
    reference = FieldSource(source, None if binding is None else binding.ref.record_id,
                            None if binding is None else binding.ref.record_version)
    result = TensorField(CoordinateComplex.from_rex(source).spaces[1], coefficients,
        (CoordinateSpace(value["axes"][1], tuple(
            f"{value['dataset']}:{value['axes'][1]}:{value['offset'][1] + i}" for i in selected)),),
        reference, 1, "chain")
    record_method("gadget-native-data-field", dataset=value["dataset"], field=value["field_id"],
                  original_rows=(value["offset"][0], value["offset"][0] + value["shape"][0]),
                  original_columns=tuple(value["offset"][1] + i for i in selected))
    return result


def install():
    """Install one declared application operator through existing RCQL registries."""
    from rcql import Domain, Exactness, OperatorSignature, RCType, TypePattern, ValueKind
    from rcql.operators import get_operator, register as implementation
    from rcql.signatures import catalogued, lookup, register
    name, key = "GADGET_DATA_FIELD", "gadget.native-data-field.v1"
    signature = OperatorSignature(name=name, source_kind=ValueKind.REX,
        inputs=(TypePattern("columns"),),
        result=RCType("TensorField", kind=ValueKind.TENSOR_FIELD, grade=1,
                      exactness=Exactness.RATIONAL, domain=Domain.RATIONAL),
        implementation_key=key, requires=frozenset({"read"}), memoizable=True,
        preconditions=("version bound native field; explicit local columns; observed selected support",))
    try:
        operator = get_operator(name)
    except KeyError:
        operator = None
    if name in catalogued():
        if lookup(name) != signature or operator is None or getattr(operator.fn, "_gadget_data_binding", None) != key:
            raise ValueError("GADGET_DATA_FIELD already has a different provider")
        return
    if operator is not None:
        raise ValueError("GADGET_DATA_FIELD implementation has no matching declaration")
    field._gadget_data_binding = key
    register(signature)
    implementation(name)(field)
