"""Original field windows verified against explicitly admitted complete parents."""
from dataclasses import dataclass

from ._records import natural, text


@dataclass(frozen=True)
class QueryScope:
    """Declare dependency scope and distinct evidence and parameter revisions."""
    extent: str
    evidence_revision: int
    parameter_revision: int

    def __post_init__(self):
        if self.extent not in ("connected", "whole-declared-source"):
            raise ValueError("query scope must be connected or whole-declared-source")
        natural(self.evidence_revision, "evidence revision")
        natural(self.parameter_revision, "parameter revision")

    def as_record(self):
        return {"extent": self.extent, "evidence_revision": self.evidence_revision,
                "parameter_revision": self.parameter_revision}


class SourceWindow:
    """Retain measured fields, support and original coordinates for a source window.

    This restricts source data before a model is built. It does not project a
    completed model family. A restored declaration needs its actual parent store
    for verification and supplies no execution authority by itself.
    """
    def __init__(self, record):
        from rexgraph.tensor_field import FieldSource
        declaration = record.get_metadata(1, 0, "gadget_source_window")
        if (type(declaration) is not dict or declaration.get("schema") != "gadget.source-window.v1"
                or set(declaration) != {"schema", "parent", "field", "rows", "columns", "time_axis", "scope"}):
            raise ValueError("complete original source window declaration required")
        QueryScope(**declaration["scope"])
        self._record = record
        self._source = FieldSource(record)

    @classmethod
    def create(cls, parent, field, *, rows, columns, scope, time_axis=None):
        if not isinstance(scope, QueryScope):
            raise TypeError("source window requires an explicit QueryScope")
        text(field, "source field")
        selected = parent.select(field, rows=rows, columns=columns)
        return cls._selected(parent, field, selected, scope, time_axis)

    @classmethod
    def _selected(cls, parent, field, selected, scope, time_axis):
        if time_axis is not None and time_axis not in selected.specification["axes"]:
            raise ValueError("time axis must name an original declared field axis")
        record = selected.native()
        record.attach_metadata(1, 0, "gadget_source_window", {
            "schema": "gadget.source-window.v1", "parent": dict(parent.reference),
            "field": field, "rows": selected.rows, "columns": selected.columns,
            "time_axis": time_axis, "scope": scope.as_record()})
        return cls(record)

    @property
    def coefficient_digest(self):
        self._source.check()
        return self._source.state_digest

    def to_record(self):
        from rexgraph.state import from_state, to_state
        self._source.check()
        return from_state(to_state(self._record))

    @classmethod
    def from_record(cls, record):
        return cls(record)

    def verify(self, parent):
        """Reconstruct all selected values and support from the actual full parent."""
        self._source.check()
        d = self._record.get_metadata(1, 0, "gadget_source_window")
        if parent.reference != d["parent"]:
            raise ValueError("source window requires its original complete parent")
        selected = parent.select(d["field"], rows=d["rows"], columns=d["columns"])
        rebuilt = type(self)._selected(parent, d["field"], selected,
            QueryScope(**d["scope"]), d["time_axis"])
        if rebuilt.coefficient_digest != self.coefficient_digest:
            raise ValueError("source window differs from its original fields or support")
        return selected
