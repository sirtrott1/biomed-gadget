"""Stream physical occurrences while retaining each complete native request once."""
from hashlib import sha256

from ._records import graph
from .execution import native_digest, parameter_digest
from .store import DataStore, DataWriter, Limits
from .workflow import record_link, require_owner_links
from ._runtime import ProgramRuntime


class ProgramJournal:
    """Reuse complete read requests under current RCQL source and policy checks.

    The request index lives in RCDB. Pending physical rows have a separate byte
    admission and flush to bounded native blocks. Native inputs supplied through
    base_parameters remain explicit caller bindings during cold replay.
    """
    def __init__(self, output, program, *, sources, contributors, base_parameters=None, limits=None,
                 check_sources=None, owners=()):
        from rcql import Program
        if not isinstance(program, Program):
            raise TypeError("journal requires an RCQL Program")
        self._check_sources = check_sources
        if self._check_sources is not None:
            self._check_sources()
        self.limits = limits or Limits()
        self.owners = tuple(owners)
        self.program, self.sources = program, sources
        self.base_parameters = dict(base_parameters or {})
        self._base_digest = parameter_digest(self.base_parameters, self.limits)
        self._implementation = ProgramRuntime()
        inputs = {}
        for store, reference in contributors:
            link = record_link(store, reference, owners=self.owners)
            key = (store.identity, reference["id"], reference["version"])
            if key in inputs and inputs[key][1] != link:
                raise ValueError("conflicting journal contributor")
            inputs[key] = (store, link)
        self.links = tuple(link for _, link in inputs.values())
        self.datasets = tuple(dict.fromkeys(store for store, _ in inputs.values()))
        self.writer = DataWriter(output, name=program.name, limits=self.limits)
        try:
            self.definition = self.writer.native("program", program, dependencies=self.links)
        except BaseException:
            self.writer.__exit__(*__import__("sys").exc_info())
            raise
        self.rows = self.requests = self.blocks = 0
        self._pending, self._pending_bytes = [], 0
        self._finished = False

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        return self.writer.__exit__(*exception)

    def _parameters(self, parameters):
        if set(parameters) & self.base_parameters.keys():
            raise ValueError("request parameters cannot replace shared native inputs")
        if parameter_digest(self.base_parameters, self.limits) != self._base_digest:
            raise ValueError("journal shared parameters changed")
        return dict(self.base_parameters, **parameters)

    def append(self, identity, *, parameters, occurrence=None):
        """Keep every supplied physical occurrence, including repeated identities."""
        from rcql import Executor
        from rexgraph.value_codec import pack_value, unpack_value
        if self._finished:
            raise ValueError("journal is sealed")
        if self._check_sources is not None:
            self._check_sources()
        supplied = self._parameters(parameters)
        executor = Executor(sources=self.sources, params=supplied)
        topology = self.program.topology(executor)
        if any(not node["reusable"] for node in topology.declaration["nodes"]):
            raise ValueError("journal reuse requires RCQL reusable read operations")
        implementation = self._implementation.check()
        digest = parameter_digest(supplied, self.limits)
        key = sha256(pack_value((self.program.coefficient_digest, topology.coefficient_digest,
                                  digest, implementation))).hexdigest()
        identifier = "objects/requests/" + key
        row = self.writer.store.get_record(identifier)
        if row is None:
            result = executor.execute_program(self.program)
            require_owner_links(result.dependencies, self.links)
            if (parameter_digest(supplied, self.limits) != digest
                    or self.program.topology(executor).coefficient_digest != topology.coefficient_digest
                    or self._implementation.check() != implementation):
                raise ValueError("journal request changed during execution")
            exports = tuple(self.writer.native("values/" + key + "/" + str(i), value,
                           dependencies=self.links) for i, value in enumerate(result.values))
            trace = graph("execution", ("execution",))
            trace.attach_metadata(1, 0, "rcql_execution",
                tuple((name, step.execution, step.provenance) for name, step in result.steps))
            trace_ref = self.writer.native("traces/" + key, trace, dependencies=self.links)
            plan = self.writer.native("topologies/" + key, topology.to_record(), dependencies=self.links)
            request = graph("request", (key,))
            request.attach_metadata(1, 0, "gadget_journal_request", {
                "key": key, "parameters": parameters, "parameter_digest": digest,
                "topology": plan, "topology_digest": topology.coefficient_digest,
                "trace": trace_ref, "trace_digest": native_digest(trace), "outputs": exports,
                "output_digests": tuple(native_digest(value) for value in result.values)})
            reference = self.writer.native("requests/" + key, request, dependencies=self.links)
            self.requests += 1
        else:
            reference = {"id": row.id, "version": row.version,
                         "state_digest": row.envelope.object_digest,
                         "native_bytes": row.meta["native_bytes"], "store_id": self.writer.store.store_id}
        if self._check_sources is not None:
            self._check_sources()
        physical = {"ordinal": self.rows, "identity": identity,
                    "occurrence": occurrence, "request": reference}
        raw = pack_value(physical)
        physical = unpack_value(raw)
        size = len(raw)
        maximum = min(self.limits.record_bytes, self.limits.decode_bytes) // 4
        if size > maximum:
            raise MemoryError("journal occurrence exceeds its record byte admission")
        if self._pending and (len(self._pending) == 128 or self._pending_bytes + size > maximum):
            self._flush()
        self._pending.append(physical)
        self._pending_bytes += size
        self.rows += 1
        return reference

    def _flush(self):
        if not self._pending:
            return
        record = graph("occurrences", ("block",))
        record.attach_metadata(1, 0, "gadget_journal_rows", tuple(self._pending))
        self.writer.native("occurrences/" + str(self.blocks), record, dependencies=self.links)
        self.blocks += 1
        self._pending, self._pending_bytes = [], 0

    def finish(self):
        """Seal the program, unique complete requests and ordered physical rows."""
        if self._finished:
            raise ValueError("journal is sealed")
        if self._check_sources is not None:
            self._check_sources()
        self._flush()
        record = graph("journal", ("journal",))
        record.attach_metadata(1, 0, "gadget_journal", {
            "schema": "gadget.program-journal.v1", "program": self.definition,
            "program_digest": self.program.coefficient_digest,
            "base_parameter_digest": self._base_digest,
            "implementation": self._implementation.check(),
            "rows": self.rows, "requests": self.requests, "blocks": self.blocks})
        reference = self.writer.native("journal", record, dependencies=self.links)
        if self._check_sources is not None:
            self._check_sources()
        self.writer.finish()
        self._finished = True
        return reference


def journal_occurrences(store, *, dependencies=(), owners=()):
    """Read bounded blocks in their original publication and occurrence order."""
    ordinal = 0
    for reference in store.manifest["objects"]:
        if not reference["id"].startswith("objects/occurrences/"):
            continue
        record = store.native(reference, dependencies=dependencies, owners=owners)
        for row in record.get_metadata(1, 0, "gadget_journal_rows"):
            if row["ordinal"] != ordinal:
                raise ValueError("journal physical occurrence order differs")
            ordinal += 1
            yield row


def verify_journal(path, *, sources, datasets, base_parameters=None, limits=None, check_sources=None, owners=()):
    """Reexecute every unique complete request and check all physical row links."""
    from rcql import Executor
    from rcql.plan_topology import PlanTopology
    from rexgraph.value_codec import pack_value
    if check_sources is not None:
        check_sources()
    limits = limits or Limits()
    datasets = tuple(datasets)
    owners = tuple(owners)
    base_parameters = dict(base_parameters or {})
    with DataStore(path, limits=limits) as store:
        summary = next(ref for ref in store.manifest["objects"] if ref["id"] == "objects/journal")
        header = store.native(summary, dependencies=datasets, owners=owners).get_metadata(1, 0, "gadget_journal")
        if header["schema"] != "gadget.program-journal.v1":
            raise ValueError("native program journal required")
        implementation = ProgramRuntime()
        if implementation.check() != header["implementation"]:
            raise ValueError("journal implementation differs; explicit migration required")
        if parameter_digest(base_parameters, limits) != header["base_parameter_digest"]:
            raise ValueError("journal shared parameters differ")
        program = store.native(header["program"], dependencies=datasets, owners=owners)
        if program.coefficient_digest != header["program_digest"]:
            raise ValueError("journal program differs")
        requests = 0
        for reference in store.manifest["objects"]:
            if not reference["id"].startswith("objects/requests/"):
                continue
            if check_sources is not None:
                check_sources()
            request = store.native(reference, dependencies=datasets, owners=owners).get_metadata(1, 0, "gadget_journal_request")
            if base_parameters.keys() & request["parameters"].keys():
                raise ValueError("journal request replaces shared inputs")
            parameters = dict(base_parameters, **request["parameters"])
            executor = Executor(sources=sources, params=parameters)
            topology = program.topology(executor)
            digest = parameter_digest(parameters, limits)
            key = sha256(pack_value((program.coefficient_digest, topology.coefficient_digest,
                                      digest, implementation.check()))).hexdigest()
            if (key != request["key"] or reference["id"] != "objects/requests/" + key
                    or digest != request["parameter_digest"]
                    or topology.coefficient_digest != request["topology_digest"]
                    or PlanTopology.from_record(store.native(request["topology"], dependencies=datasets, owners=owners)).coefficient_digest
                        != topology.coefficient_digest):
                raise ValueError("journal source, scope, topology or request differs")
            result = executor.execute_program(program)
            if native_digest(store.native(request["trace"], dependencies=datasets, owners=owners)) != request["trace_digest"]:
                raise ValueError("journal execution evidence differs")
            if len(result.values) != len(request["outputs"]):
                raise ValueError("journal export count differs")
            for ref, actual, expected in zip(request["outputs"], result.values, request["output_digests"], strict=True):
                if native_digest(actual) != expected or native_digest(store.native(ref, dependencies=datasets, owners=owners)) != expected:
                    raise ValueError("journal complete native output differs")
            requests += 1
        rows = 0
        for row in journal_occurrences(store, dependencies=datasets, owners=owners):
            ref = row["request"]
            if not ref["id"].startswith("objects/requests/"):
                raise ValueError("journal occurrence has no complete request")
            store.read(ref)
            rows += 1
        if (rows, requests) != (header["rows"], header["requests"]):
            raise ValueError("journal occurrence or request count differs")
        implementation.check()
        if check_sources is not None:
            check_sources()
        return {"rows": rows, "requests": requests}
