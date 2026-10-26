"""Capture and replay finite native programs through existing RCQL and RCDB APIs."""
from pathlib import Path

from ._records import graph
from .store import DataStore, DataWriter, Limits
from .workflow import record_link, require_owner_links
from ._runtime import ProgramRuntime


def native_digest(value):
    """Identify a complete native value rather than a rendered observation."""
    from rcql import Program
    from rexgraph import RexGraph
    from rexgraph.object_identity import object_digest
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    from rexgraph.model_codec import MODEL_VALUES
    if isinstance(value, RexGraph):
        return object_digest(value)
    if isinstance(value, (Program, TensorField, SectionFamily, *MODEL_VALUES)):
        return value.coefficient_digest
    from hashlib import sha256
    from rcql import Result
    return sha256(Result((value,)).to_bytes()).hexdigest()


def parameter_digest(parameters, limits):
    from hashlib import sha256
    from rcql import Result
    if not isinstance(parameters, dict) or any(type(key) is not str for key in parameters):
        raise TypeError("program parameters require a mapping with string names")
    raw = Result((dict(sorted(parameters.items())),)).to_bytes(max_bytes=limits.decode_bytes)
    return sha256(raw).hexdigest()


def capture_program(program, output, *, sources, parameters, contributors, limits=None, check_sources=None,
                    owners=()):
    """Execute explicit bindings and retain complete outputs and their owners.

    Contributors are admitted (DataStore, record reference) pairs. RCQL keeps
    responsibility for program effects, source policies and argument contracts.
    """
    from rcql import Executor, Program
    if check_sources is not None:
        check_sources()
    if not isinstance(program, Program):
        raise TypeError("capture requires an RCQL Program")
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    limits = limits or Limits()
    owners = tuple(owners)
    inputs = {}
    for store, reference in contributors:
        key = (store.identity, reference["id"], reference["version"])
        link = record_link(store, reference, owners=owners)
        if key in inputs and inputs[key][1] != link:
            raise ValueError("conflicting contributor identity")
        inputs[key] = (store, link)
    links = tuple(link for _, link in inputs.values())
    datasets = tuple(dict.fromkeys(store for store, _ in inputs.values()))
    params = parameter_digest(parameters, limits)
    runtime = ProgramRuntime()
    executor = Executor(sources=sources, params=parameters)
    topology = program.topology(executor).to_record()
    result = executor.execute_program(program)
    if check_sources is not None:
        check_sources()
    require_owner_links(result.dependencies, links)
    identities = tuple(native_digest(value) for value in result.values)
    if parameter_digest(parameters, limits) != params:
        raise ValueError("program parameters changed during execution")
    if native_digest(program.topology(executor).to_record()) != native_digest(topology):
        raise ValueError("program source or topology changed during execution")
    runtime.check()
    execution = graph("execution", ("execution",))
    execution.attach_metadata(1, 0, "rcql_execution",
        tuple((name, step.execution, step.provenance) for name, step in result.steps))
    with DataWriter(output, name=program.name, limits=limits) as writer:
        definition = writer.native("program", program, dependencies=links)
        plan = writer.native("topology", topology, dependencies=links)
        trace = writer.native("trace", execution, dependencies=links)
        exports = tuple(writer.native("output/" + str(i), value, dependencies=links)
                        for i, value in enumerate(result.values))
        receipt = {"schema": "gadget.native-execution.v1", "program": definition,
                   "program_digest": program.coefficient_digest, "parameter_digest": params,
                   "topology": plan, "topology_digest": native_digest(topology),
                   "trace": trace, "trace_digest": native_digest(execution),
                   "outputs": exports, "output_digests": identities,
                   "sources": result.sources, "aliases": result.aliases,
                   "implementation": runtime.check()}
        record = graph("execution receipt", ("receipt",))
        record.attach_metadata(1, 0, "gadget_execution", receipt)
        reference = writer.native("execution", record, dependencies=links)
        if check_sources is not None:
            check_sources()
        writer.finish()
    with DataStore(output, limits=limits) as stored:
        restored = stored.native(reference, dependencies=datasets, owners=owners)
        if restored.get_metadata(1, 0, "gadget_execution") != receipt:
            raise ValueError("retained execution receipt differs")
        for ref, expected in ((definition, program.coefficient_digest),
                              (plan, native_digest(topology)), (trace, native_digest(execution)),
                              *zip(exports, identities, strict=True)):
            if native_digest(stored.native(ref, dependencies=datasets, owners=owners)) != expected:
                raise ValueError("retained program value differs")
    return reference


def replay_program(path, *, sources, parameters, datasets, limits=None, check_sources=None, owners=()):
    """Restore explicit dependencies and compare every complete native export."""
    from rcql import Executor
    if check_sources is not None:
        check_sources()
    limits = limits or Limits()
    datasets = tuple(datasets)
    owners = tuple(owners)
    with DataStore(path, limits=limits) as stored:
        candidates = [ref for ref in stored.manifest["objects"] if ref["id"] == "objects/execution"]
        if len(candidates) != 1:
            raise ValueError("one retained program execution is required")
        receipt = stored.native(candidates[0], dependencies=datasets, owners=owners).get_metadata(1, 0, "gadget_execution")
        if receipt["schema"] != "gadget.native-execution.v1":
            raise ValueError("native program execution required")
        runtime = ProgramRuntime()
        if runtime.check() != receipt["implementation"]:
            raise ValueError("program implementation differs; explicit migration required")
        if parameter_digest(parameters, limits) != receipt["parameter_digest"]:
            raise ValueError("replay parameters differ")
        program = stored.native(receipt["program"], dependencies=datasets, owners=owners)
        if program.coefficient_digest != receipt["program_digest"]:
            raise ValueError("replay program differs")
        executor = Executor(sources=sources, params=parameters)
        topology = program.topology(executor).to_record()
        if native_digest(topology) != receipt["topology_digest"]:
            raise ValueError("replay source, policy or topology differs")
        for name in ("topology", "trace"):
            if native_digest(stored.native(receipt[name], dependencies=datasets, owners=owners)) != receipt[name + "_digest"]:
                raise ValueError("retained execution evidence differs")
        result = executor.execute_program(program)
        if (result.sources != tuple(receipt["sources"]) or result.aliases != tuple(receipt["aliases"])
                or len(result.values) != len(receipt["outputs"])):
            raise ValueError("replay source bindings or exports differ")
        for ref, value, expected in zip(receipt["outputs"], result.values, receipt["output_digests"], strict=True):
            if (native_digest(value) != expected
                    or native_digest(stored.native(ref, dependencies=datasets, owners=owners)) != expected):
                raise ValueError("replay changed a complete native export")
        runtime.check()
        if check_sources is not None:
            check_sources()
        return result
