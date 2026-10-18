"""Literal lookup across independent versioned native dataset axes."""
from ._records import declaration, graph, index_term, text
from .store import DataStore, DataWriter
from .workflow import check_link, record_link


def _axis_metadata(link, stores):
    check_link(link)
    source = stores[link["dataset"]["dataset"]]
    if source.reference != link["dataset"] or link["owner_digest"] != link["record"]["state_digest"]:
        raise ValueError("axis locator has a foreign native owner")
    return source._axis_metadata(link["record"])


def _vocabulary(labels, contexts, roles):
    terms = {index_term("literal", label) for label in labels}
    for row in contexts:
        terms.update(index_term("literal", row[role]) for role in roles if role in row)
    return terms


def _bucket(term):
    # Physical candidates are checked against complete literal terms and owners.
    return "gadgetliteral" + term[-1]


def build_axis_index(inputs, output, *, context_roles=(), limits=None):
    """Index bounded original axis blocks without copying measurement arrays."""
    roles = tuple(context_roles)
    if len(set(roles)) != len(roles):
        raise ValueError("distinct context roles required")
    for role in roles:
        text(role, "lookup context role")
    with DataWriter(output, name="native axis index", limits=limits) as writer:
        from rexgraph.value_codec import pack_value
        for source, axis in inputs:
            text(axis, "indexed axis")
            for block in source.manifest["axes"][axis]:
                labels, contexts = source._axis_metadata(block)
                link = record_link(source, block)
                links = writer._links((link,))
                name = "lookup/" + source.identity + "/" + axis + "/" + str(block["offset"])
                literal_terms = _vocabulary(labels, contexts, roles)
                vocabulary = {_bucket(term) for term in literal_terms}
                lookup_bytes = len(pack_value(tuple(vocabulary))) + len(name.encode()) + 1024
                writer.limits.admit(lookup_bytes, "catalog_bytes")
                entry = writer._put(name, graph("axis_locator", ((source.identity, axis, block["offset"]),)),
                    {"kind": "axis_locator", "axis": axis, "offset": block["offset"],
                     "count": block["count"], "context_roles": roles, "dependencies": links,
                     "literal_terms": tuple(sorted(literal_terms))},
                    vocabulary=vocabulary)
                writer.objects.append(dict(entry, lookup_bytes=lookup_bytes))
        return writer.finish()


class AxisIndex(DataStore):
    """Stream original positions and matching roles through admitted dataset stores."""
    def find(self, literal, *, datasets):
        from rexgraph.value_codec import pack_value
        stores = {store.identity: store for store in datasets}
        term = index_term("literal", literal)
        bucket = _bucket(term)
        largest = max((block["lookup_bytes"] for block in self.manifest["objects"]), default=1)
        self.limits.admit(largest, "catalog_bytes")
        cap = max(1, self.limits.catalog_bytes // largest)
        found = self.store.query(limit=cap, labels_any=[bucket])
        candidates = {row.id for row in found} if len(found) < cap else None
        for block in self.manifest["objects"]:
            if candidates is not None and block["id"] not in candidates:
                current = self.store.get_record(block["id"])
                if current is not None and current.version == block["version"]:
                    continue
            if candidates is None:
                row = next(row for row in self.store.history(block["id"]) if row.version == block["version"])
                if bucket not in row.meta["vertex_labels"]:
                    continue
            locator = declaration(self.read(block))
            if locator["kind"] != "axis_locator":
                raise ValueError("axis index contains a different native record kind")
            if term not in locator["literal_terms"]:
                continue
            link, = locator["dependencies"]
            labels, contexts = _axis_metadata(link, stores)
            for i, (label, row) in enumerate(zip(labels, contexts, strict=True)):
                roles = (("label",) if pack_value(label) == pack_value(literal) else ()) + tuple(
                    role for role in locator["context_roles"]
                    if role in row and pack_value(row[role]) == pack_value(literal))
                if roles:
                    from copy import deepcopy
                    result = {"coordinate": locator["offset"] + i, "roles": roles,
                              "label": label, "axis": locator["axis"], "source": link}
                    self._admit_metadata(result, 0)
                    yield deepcopy(result)

    def verify_index(self, *, datasets):
        """Compare every declared lookup vocabulary against original native axes."""
        stores = {store.identity: store for store in datasets}
        count = positions = 0
        for block in self.manifest["objects"]:
            locator = declaration(self.read(block))
            if locator["kind"] != "axis_locator":
                raise ValueError("native axis locator required")
            link, = locator["dependencies"]
            labels, contexts = _axis_metadata(link, stores)
            if (locator["axis"], locator["offset"], locator["count"]) != (
                    link["record"]["axis"], link["record"]["offset"], len(labels)):
                raise ValueError("locator coordinates differ from the original axis")
            row = next(row for row in self.store.history(block["id"]) if row.version == block["version"])
            terms = _vocabulary(labels, contexts, locator["context_roles"])
            if (set(locator["literal_terms"]) != terms or
                    set(row.meta["vertex_labels"]) != {_bucket(term) for term in terms}):
                raise ValueError("lookup vocabulary differs from complete original labels and roles")
            count += 1
            positions += len(labels)
        return {"locators": count, "original_positions": positions}
