"""Explicit native dataset owners for checked query and program execution."""
from contextlib import contextmanager
from threading import RLock

from .workflow import native_owner, resolve_link


class NativeBindings:
    """Retain admitted owners without opening paths from stored declarations.

    RexGraph owners can supply RCQL sources. C++ owners supply the original
    dependencies of model parameters through the canonical result binder.
    Filesystem actions and training require their existing separate admission.
    Bindings and active executions retain their owners until released.
    """
    def __init__(self):
        self._entries = {}
        self._generation = 0
        self._lock = RLock()

    def bind(self, name, store, reference, *, dependencies=(), permissions=("read", "identity")):
        from rcql import SourcePolicy
        from .store import DataStore
        if type(name) is not str or not name.isascii() or not name.isidentifier() or name == "db":
            raise ValueError("native source requires an ASCII identifier other than db")
        if not isinstance(store, DataStore):
            raise TypeError("native binding requires an admitted DataStore")
        dependencies = tuple(dependencies)
        if any(not isinstance(item, DataStore) for item in dependencies):
            raise TypeError("native dependencies require admitted DataStore values")
        permissions = frozenset(permissions)
        if permissions - {"read", "identity"}:
            raise ValueError("native bindings permit read and identity only")
        with self._lock:
            if name in self._entries:
                raise ValueError("native source alias is already bound")
            store.check_reference(reference)
            owner = next((owner for _, (existing, owner, parents, _) in self._entries.items()
                          if existing is store and owner.link["record"] == reference
                          and parents == dependencies), None)
            if owner is None:
                retained = tuple({id(entry[1]): entry[1] for entry in self._entries.values()}.values())
                owner = native_owner(store, reference, dependencies=dependencies, owners=retained)
            owner.value
            self._entries[name] = (store, owner, dependencies, SourcePolicy.allow(*permissions))
            self._generation += 1
        return owner

    def remove(self, name):
        """Release a named owner and invalidate executions using its binding set."""
        with self._lock:
            del self._entries[name]
            self._generation += 1

    @contextmanager
    def scope(self):
        """Check current stores and actual owners around one execution."""
        with self._lock:
            scope = _BindingScope(self, self._generation, tuple(self._entries.items()))
        scope.check()
        yield scope
        scope.check()

    def capture(self, program, output, *, parameters, limits=None):
        from .execution import capture_program
        with self.scope() as scope:
            return capture_program(program, output, sources=scope.sources(),
                parameters=scope.parameters(parameters), contributors=scope.contributors,
                limits=limits, check_sources=scope.check, owners=scope.owners)

    def replay(self, path, *, parameters, limits=None):
        from .execution import replay_program
        with self.scope() as scope:
            return replay_program(path, sources=scope.sources(), parameters=scope.parameters(parameters),
                datasets=tuple(scope.stores.values()), limits=limits, check_sources=scope.check, owners=scope.owners)

    def journal(self, program, output, *, base_parameters=None, limits=None):
        from .journal import ProgramJournal
        with self.scope() as scope:
            return ProgramJournal(output, program, sources=scope.sources(), contributors=scope.contributors,
                base_parameters=scope.parameters(base_parameters or {}), limits=limits,
                check_sources=scope.check, owners=scope.owners)

    def verify_journal(self, path, *, base_parameters=None, limits=None):
        from .journal import verify_journal
        with self.scope() as scope:
            return verify_journal(path, sources=scope.sources(), datasets=tuple(scope.stores.values()),
                base_parameters=scope.parameters(base_parameters or {}), limits=limits,
                check_sources=scope.check, owners=scope.owners)


class _BindingScope:
    def __init__(self, bindings, generation, entries):
        from rcql import SourcePolicy
        self._bindings, self._generation, self.entries = bindings, generation, entries
        self.stores, self._owners = {}, {}
        for _, (store, owner, dependencies, policy) in entries:
            for dataset in (store, *dependencies):
                if dataset.identity in self.stores and self.stores[dataset.identity] is not dataset:
                    raise ValueError("one admitted store per native dataset required")
                self.stores[dataset.identity] = dataset
            key = (store.identity, owner.link["record"]["id"], owner.link["record"]["version"])
            if key in self._owners:
                previous = self._owners[key]
                policy = SourcePolicy.intersection(previous[2], policy)
            self._owners[key] = (store, owner, policy)
        self.contributors = tuple((store, owner.link["record"])
                                  for store, owner, policy in self._owners.values() if policy.permits("read"))
        self.owners = tuple(owner for _, owner, policy in self._owners.values() if policy.permits("read"))

    def check(self):
        with self._bindings._lock:
            if self._bindings._generation != self._generation:
                raise ValueError("native source bindings changed during execution")
        owners = tuple(owner for _, owner, _ in self._owners.values())
        for _, owner, _ in self._owners.values():
            resolve_link(owner.link, self.stores, owners)

    def sources(self):
        from rcql import BoundSource
        from rcql.types import SourceRef
        from rexgraph import RexGraph
        self.check()
        result = {}
        for name, (_, owner, _, policy) in self.entries:
            value = owner.value
            if isinstance(value, RexGraph):
                pin = owner.link["record"]
                ref = SourceRef(name, state_digest=owner.link["owner_digest"], policy_digest=policy.digest,
                                record_id=pin["id"], record_version=pin["version"])
                result[name] = BoundSource(value, policy, ref=ref)
        return result

    def parameters(self, parameters):
        from rcql import BoundSource
        from rcql.result_codec import bind_sources
        from rcql.source_context import field_references
        from rexgraph.source_binding import CppSourceReference
        self.check()
        if not self.entries:
            return parameters
        references = {}
        for reference in field_references(parameters):
            matches = [(owner, policy) for _, owner, policy in self._owners.values()
                if owner.link["owner_digest"] == reference.state_digest and (
                    isinstance(reference, CppSourceReference) or reference.record_id is None or
                    (reference.record_id, reference.version) ==
                        (owner.link["record"]["id"], owner.link["record"]["version"]))]
            known = any(owner.link["owner_digest"] == reference.state_digest
                        for _, owner, _ in self._owners.values())
            if not matches and not known and reference.source is not None:
                continue
            if len(matches) != 1:
                raise ValueError("native parameter needs one original admitted owner")
            owner, policy = matches[0]
            BoundSource(owner.value, policy).require("read")
            references[reference.state_digest, reference.record_id, reference.version] = owner.value
        return bind_sources(parameters, references)
