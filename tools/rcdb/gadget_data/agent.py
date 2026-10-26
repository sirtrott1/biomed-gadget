"""GADGET publications and native queries through the optional RexGraph Agent runtime."""
from collections.abc import Mapping


class GadgetSession:
    """Bind GADGET's existing providers to an existing scoped Agent session.

    Supply an executable and a publication byte budget explicitly. Native
    execution, replay and field selection require explicit capabilities. The object
    is also a callable Hive worker accepting ``query`` and its named arguments.
    GADGET owns native admission and execution; RCDB owns record publication.
    """

    def __init__(self, session, *, gadget, max_bytes, capabilities=(), selector=None):
        try:
            from agent.relational import RelationalSession
        except ImportError as exc:
            raise ImportError("Agent integration requires a matching RexGraph Agent installation "
                              "with RelationalSession support; install the RexGraph packages "
                              "before biomed-gadget-data[agent]") from exc
        if not isinstance(session, RelationalSession):
            raise TypeError("GADGET requires an existing RelationalSession")
        capabilities = frozenset(capabilities)
        if capabilities - {"execute", "replay", "select"}:
            raise ValueError("GADGET capabilities are execute, replay and select")
        if type(max_bytes) is not int or max_bytes < 1:
            raise ValueError("GADGET requires a positive publication byte budget")
        if "select" in capabilities and selector is None:
            raise ValueError("native field selection requires an explicit selector executable")
        try:
            import gadget_rcdb
        except ImportError as exc:
            raise ImportError("GADGET requires the optional biomed-gadget-data package; "
                              "install biomed-gadget-data[agent]") from exc
        from rcql import BoundSource, SourcePolicy

        self.session = session
        self.max_bytes = max_bytes
        self._adapter = gadget_rcdb
        from pathlib import Path
        self.gadget = str(Path(gadget).resolve(strict=True))
        self.codec = gadget_rcdb.install_codec(gadget=self.gadget, max_bytes=max_bytes)
        gadget_rcdb.install_rcql_operator()
        permissions = {"read", "identity", "records"}
        if "execute" in capabilities:
            gadget_rcdb.install_rcql_executor(gadget=self.gadget)
            gadget_rcdb.install_rcql_recipes(gadget=self.gadget, max_bytes=max_bytes)
            permissions.update(("files", "gadget_execute"))
        if "replay" in capabilities:
            gadget_rcdb.install_rcql_scientific_replay(gadget=self.gadget)
            permissions.update(("files", "gadget_replay"))
        if "select" in capabilities:
            gadget_rcdb.install_rcql_field_selector(selector=selector, gadget=self.gadget)
            permissions.update(("files", "gadget_field_select"))
        self._permissions = frozenset(permissions)
        from .bindings import NativeBindings
        from .query import install
        install()
        self._bindings = NativeBindings()
        self._database = BoundSource(session.store, SourcePolicy.allow(*permissions))

    @property
    def bindings(self):
        return self._bindings

    def bind_native(self, name, store, reference, *, dependencies=(), permissions=("read", "identity")):
        """Bind actual dataset owners through trusted caller code, outside worker arguments."""
        return self.bindings.bind(name, store, reference, dependencies=dependencies, permissions=permissions)

    def unbind_native(self, name):
        self.bindings.remove(name)

    def _sources(self, scope):
        return {"db": self._database, **scope.sources()}

    def _native_runtime(self, scope):
        from agent.rcql_runtime import RCQLRuntime
        runtime = RCQLRuntime()
        for name, value in self._sources(scope).items():
            runtime.register(name, value)
        return runtime

    def operation_capabilities(self, names):
        """Inspect current contracts and declared grants without opening field data."""
        from .capabilities import operation_capabilities
        from .query import install
        install()
        return {"operations": operation_capabilities(names),
                "declared_permissions": sorted(self._permissions)}

    def recipe_worker(self, hive, model, *, workspace, intent, execute=True, generation=None):
        """Propose catalog goals through Hive and submit them to native admission."""
        from .agent_models import GadgetRecipeWorker
        return GadgetRecipeWorker(self, hive, model, workspace=workspace, intent=intent,
                                  execute=execute, generation=generation)

    def publish(self, record_id, publication, *, expected_version, meta=None):
        """Admit complete RFCXS bytes and publish under the caller's store scope."""
        from agent.relational import _integer
        _integer(expected_version, "expected_version")
        if type(publication) is not bytes or len(publication) > self.max_bytes:
            raise ValueError("GADGET publication requires bytes within the declared budget")
        value = self._adapter.validate_publication(publication, gadget=self.gadget)
        return self.session.store.put_record(record_id, value, codec=self.codec,
            expected_version=expected_version, meta=meta)

    def publish_file(self, record_id, path, *, expected_version, meta=None):
        """Use GADGET's existing bounded file admission and RCDB publication."""
        from agent.relational import _integer
        _integer(expected_version, "expected_version")
        return self._adapter.put_publication(self.session.store, record_id, path,
            expected_version=expected_version, meta=meta)

    def query(self, query, *, parameters=None, record_parameters=None,
              save_as=None, expected_version=0):
        """Execute a typed read/action query with explicitly granted capabilities."""
        from dataclasses import replace

        from rcql import Query, parse
        query = parse(query) if isinstance(query, str) else query
        if not isinstance(query, Query):
            raise TypeError("GADGET query requires an ordinary RCQL Query")
        supplied, references = self.session._parameters(parameters, record_parameters)
        with self.bindings.scope() as scope:
            result = self._native_runtime(scope).execute(query, params=scope.parameters(supplied))
        if references:
            result = replace(result, native_plan={**(result.native_plan or {}),
                                                 "agent_record_parameters": references})
        if save_as is not None:
            self.session.save_result(save_as, result, expected_version=expected_version)
        return result

    def run_program(self, record_id, *, version, parameters=None, record_parameters=None,
                    save_as=None, expected_version=0):
        """Run a stored read program; filesystem actions retain their own boundary."""
        from rcql import Program

        from agent.relational import _integer
        _integer(version, "version", 1)
        snapshot = self.session.read(record_id, version)
        program = Program.from_record(snapshot.value)
        supplied, references = self.session._parameters(parameters, record_parameters)
        with self.bindings.scope() as scope:
            output = self._native_runtime(scope).execute_program(program, params=scope.parameters(supplied))
        result = self.session._program_result(snapshot, program, output, references)
        if save_as is not None:
            self.session.save_result(save_as, result, expected_version=expected_version)
        return result

    def capture_program(self, record_id, output, *, version, parameters=None, record_parameters=None, limits=None):
        """Retain a stored read program and complete outputs under current native bindings."""
        from rcql import Program
        from agent.relational import _integer
        from .execution import capture_program
        _integer(version, "version", 1)
        program = Program.from_record(self.session.read(record_id, version).value)
        supplied, _ = self.session._parameters(parameters, record_parameters)
        with self.bindings.scope() as scope:
            return capture_program(program, output, sources=self._sources(scope),
                parameters=scope.parameters(supplied), contributors=scope.contributors,
                limits=limits, check_sources=scope.check, owners=scope.owners)

    def replay_program(self, path, *, parameters=None, record_parameters=None, limits=None):
        """Replay retained read outputs with current session and native dataset access."""
        from .execution import replay_program
        supplied, _ = self.session._parameters(parameters, record_parameters)
        with self.bindings.scope() as scope:
            return replay_program(path, sources=self._sources(scope), parameters=scope.parameters(supplied),
                datasets=tuple(scope.stores.values()), limits=limits, check_sources=scope.check, owners=scope.owners)

    def __call__(self, arguments):
        """Run as an ordinary capability registered with ``Hive.add_worker``."""
        if (not isinstance(arguments, Mapping) or "query" not in arguments
                or set(arguments) - {"query", "parameters", "record_parameters", "save_as", "expected_version"}):
            raise ValueError("GADGET worker requires query and declared query arguments")
        return self.query(**arguments)
