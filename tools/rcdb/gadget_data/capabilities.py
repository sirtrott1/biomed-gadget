"""Inspect installed RCQL contracts without granting access or executing queries."""


def operation_capabilities(names):
    """Report adapter contracts separately from execution and qualification.

    Actual argument types, source owners and permissions remain RCQL's execution
    boundary. A registered name supplies no platform or biomedical qualification.
    """
    from rcql.inventory import operator_inventory
    names = tuple(names)
    if any(type(name) is not str or not name for name in names):
        raise TypeError("operation names require nonempty strings")
    requested = set(names)
    rows = {}
    offset = 0
    while requested - rows.keys():
        page = operator_inventory(limit=1000, offset=offset)
        if not page:
            break
        for row in page:
            if row["name"] in requested:
                rows[row["name"]] = row
        offset += len(page)
    result = []
    for name in names:
        contract = rows.get(name)
        result.append({"name": name, "contract": contract,
                       "status": "unavailable" if contract is None else "declared",
                       "execution": "requires admitted arguments and source permissions",
                       "platform_qualification": "unverified",
                       "biomedical_qualification": "unverified"})
    return result


def native_primitive_capabilities():
    """Return the provider's scoped constructor and action ledger unchanged."""
    from importlib.util import find_spec
    if find_spec("rcql.c_workflow") is None:
        return {"available": False, "primitives": ()}
    from rcql.c_workflow import primitive_ledger
    return {"available": True, "primitives": primitive_ledger()}
