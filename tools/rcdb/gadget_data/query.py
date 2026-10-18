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
    """Install declared application operations through existing RCQL registries."""
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
    else:
        if operator is not None:
            raise ValueError("GADGET_DATA_FIELD implementation has no matching declaration")
        field._gadget_data_binding = key
        register(signature)
        implementation(name)(field)
    from rcql.arguments import EXPRESSION_ARGUMENTS
    from rcql.names import insert_unique
    def arguments_for(name, arguments):
        expected = (arguments, ())
        if name in EXPRESSION_ARGUMENTS:
            if EXPRESSION_ARGUMENTS[name] != expected:
                raise ValueError(name + " has a conflicting expression declaration")
        else:
            insert_unique(EXPRESSION_ARGUMENTS, name, expected)
    arguments_for(name, ("columns",))
    for name, function, arguments, result in (
        ("GADGET_RECORD_OWNER", record_owner, ("digest",),
         RCType("FieldSource", kind=ValueKind.UNKNOWN, exactness=Exactness.STRUCTURAL)),
        ("GADGET_RELATION_FAMILY", relation_family, ("field",),
         RCType("SectionFamily", kind=ValueKind.SECTION_FAMILY, domain=Domain.RATIONAL,
                exactness=Exactness.STRUCTURAL)),
        ("GADGET_RETAINED_FIELD", retained_field, ("owners",),
         RCType("TensorField", kind=ValueKind.TENSOR_FIELD, domain=Domain.RATIONAL,
                exactness=Exactness.RATIONAL))):
        signature = OperatorSignature(name=name, source_kind=ValueKind.REX,
            inputs=tuple(TypePattern(arg) for arg in arguments), result=result,
            implementation_key="gadget." + name.lower() + ".v1", requires=frozenset({"read"}), memoizable=True)
        if name in catalogued():
            if lookup(name) != signature or get_operator(name).fn is not function:
                raise ValueError(name + " already has a different provider")
        else:
            try:
                get_operator(name)
            except KeyError:
                pass
            else:
                raise ValueError(name + " implementation has no matching declaration")
            register(signature)
            implementation(name)(function)
        arguments_for(name, arguments)


def record_owner(source, digest):
    """Bind only the expected native state and actual RCQL record version."""
    from rexgraph.tensor_field import FieldSource
    from rcql.execution_trace import current_binding
    binding = current_binding()
    if binding is None:
        raise ValueError("version bound query required")
    return FieldSource(source, binding.ref.record_id, binding.ref.record_version, digest)


def relation_family(source, field):
    """Ask the canonical exact provider for every relation compatible with a field."""
    from rexgraph.coordinate_map import CoordinateMap
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    if not isinstance(field, TensorField) or len(field.axes) != 1:
        raise TypeError("one retained measurement axis on an actual tensor field required")
    field.bind(source).check_state()
    equations = CoordinateMap(field.space, field.axes[0],
        tuple((j, i, value) for (i, j), value in np.ndenumerate(field.values) if value))
    rhs = TensorField(equations.codomain, [0] * equations.shape[0], source=field.source,
                      dependencies=field.dependencies)
    return SectionFamily.solve(equations, rhs, field.coefficient_digest)


def retained_field(source, owners):
    """Retrieve a stored tensor through already admitted native owner references."""
    from rexgraph.tensor_field import FieldSource, TensorField
    value = declaration(source)
    if value["kind"] != "native" or value["native_type"] != "TensorField":
        raise TypeError("retained tensor field record required")
    field = source.get_metadata(1, 0, "native_value")
    if not isinstance(field, TensorField):
        raise TypeError("actual native tensor field required")
    supplied = tuple(owners)
    bound = []
    for expected in (field.source, *field.dependencies):
        candidates = [ref for ref in supplied if isinstance(ref, FieldSource) and ref.matches(expected)]
        if len(candidates) != 1 or candidates[0].source is None:
            raise ValueError("each retained field owner needs its actual native version")
        candidates[0].check()
        bound.append(candidates[0])
    from dataclasses import replace
    return replace(field, source=bound[0], dependencies=tuple(bound[1:]))
