"""Complete RCQL provider installation, repeated setup and binding refusals."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

from rcql.operators import get_operator, register as register_operator
from rcql.signatures import catalogued, lookup, register as register_signature

module_path = Path(__file__).resolve().parents[1] / "tools/rcdb/gadget_rcdb.py"
binary, selector = (Path(p).resolve() for p in sys.argv[1:3])


def load_adapter():
    spec = importlib.util.spec_from_file_location("gadget_rcdb", module_path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    return adapter


adapter = load_adapter()
installers = {
    "GADGET_PUBLICATION_BYTES": lambda: adapter.install_rcql_operator(),
    "GADGET_ASSAY_QUERY": lambda: adapter.install_rcql_executor(gadget=binary),
    "GADGET_NATIVE_FIELD": lambda: adapter.install_rcql_field_selector(
        selector=selector, gadget=binary),
}


def rejects(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("Expected provider registration refusal")


if len(sys.argv) > 3:
    mode, name = sys.argv[3:5]
    install = installers[name]

    def foreign(*args):
        raise AssertionError("Foreign implementation must not execute")

    if mode == "function-only":
        register_operator(name)(foreign)
        rejects(install)
        assert name not in catalogued()
        assert get_operator(name).fn is foreign
    elif mode == "signature-only":
        with patch("rcql.operators.register", side_effect=RuntimeError("interrupted setup")):
            try:
                install()
            except RuntimeError:
                pass
            else:
                raise AssertionError("Expected interrupted installation")
        declaration = lookup(name)
        rejects(install)
        assert lookup(name) is declaration
        try:
            get_operator(name)
        except KeyError:
            pass
        else:
            raise AssertionError("Incomplete installation must stay explicit")
    elif mode == "foreign-binding":
        def register_foreign(operator_name):
            def decorate(fn):
                register_operator(operator_name)(foreign)
                return fn
            return decorate

        with patch("rcql.operators.register", register_foreign):
            install()
        declaration = lookup(name)
        rejects(install)
        assert lookup(name) is declaration
        assert get_operator(name).fn is foreign
    elif mode == "different-declaration":
        def register_different(signature):
            return register_signature(replace(
                signature, implementation_key=signature.implementation_key + ".other"))

        with patch("rcql.signatures.register", register_different):
            install()
        declaration, operator = lookup(name), get_operator(name)
        rejects(install)
        assert lookup(name) is declaration and get_operator(name) is operator
    else:
        raise AssertionError("Unknown registration case")
    print(name, mode, "PASS")
    sys.exit(0)

for name, install in installers.items():
    assert install() == name
    declaration, operator = lookup(name), get_operator(name)
    assert install() == name
    assert lookup(name) is declaration and get_operator(name) is operator

adapter = load_adapter()
for name, install in installers.items():
    declaration, operator = lookup(name), get_operator(name)
    assert install() == name
    assert lookup(name) is declaration and get_operator(name) is operator

aliased_binary = binary.parent / ".." / binary.parent.name / binary.name
aliased_selector = selector.parent / ".." / selector.parent.name / selector.name
assert adapter.install_rcql_executor(gadget=aliased_binary) == "GADGET_ASSAY_QUERY"
assert adapter.install_rcql_field_selector(
    selector=aliased_selector, gadget=aliased_binary) == "GADGET_NATIVE_FIELD"
for name, action in [
    ("GADGET_ASSAY_QUERY", lambda: adapter.install_rcql_executor(gadget=selector)),
    ("GADGET_NATIVE_FIELD", lambda: adapter.install_rcql_field_selector(
        selector=binary, gadget=binary)),
    ("GADGET_NATIVE_FIELD", lambda: adapter.install_rcql_field_selector(
        selector=selector, gadget=selector)),
]:
    declaration, operator = lookup(name), get_operator(name)
    rejects(action)
    assert lookup(name) is declaration and get_operator(name) is operator

for mode in ("function-only", "signature-only", "foreign-binding", "different-declaration"):
    for name in installers:
        process = subprocess.run([sys.executable, __file__, str(binary), str(selector), mode, name],
                                 capture_output=True, text=True)
        assert process.returncode == 0, (mode, name, process.stdout, process.stderr)
        print(process.stdout.strip())
print("RCQL repeated setup, resolved bindings, reload and unchanged registry refusals PASS")
