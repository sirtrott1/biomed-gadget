"""Caller-authored scientific switching retains weights and isolates prompt state."""
import json
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    binary, model, fixture, root = sys.argv[1:]
    fixture, root = Path(fixture), Path(root)
    shutil.rmtree(root, ignore_errors=True)
    shutil.copytree(fixture, root)
    declaration = json.loads((root / "study.json").read_text())
    declaration["id"] += "-rebound"
    (root / "study.json").write_text(json.dumps(declaration) + "\n")
    prompt = dict(prompt="Return a query proposal with operation inspect to inspect this study catalog.")
    requests = [prompt, dict(study=str(root / "missing.json")), prompt,
                dict(study=str(root / "study.json")), prompt,
                dict(study=str(root / "study.json"), prompt="ambiguous")]
    completed = subprocess.run([binary, "session", "--study", str(fixture / "study.json"),
                                "--model", model, "--context-tokens", "8192"],
                               input="".join(json.dumps(r) + "\n" for r in requests),
                               capture_output=True, text=True, timeout=300, check=True)
    responses = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(responses) == len(requests)
    first, failed, preserved, bound, switched, ambiguous = responses
    assert first["status"] == preserved["status"] == switched["status"] == "executed", responses
    assert failed["status"] == ambiguous["status"] == "rejected"
    assert first["source_digest"] == preserved["source_digest"] != switched["source_digest"]
    assert first["model"] == preserved["model"] == switched["model"]
    assert preserved["performance"]["prefix_tokens_reused"] > 0
    assert switched["performance"]["prefix_tokens_reused"] == 0
    assert bound["schema"] == "gadget.study.bound.v1" and bound["status"] == "bound"
    assert bound["source_digest"] == switched["source_digest"]
    print("PASS: CLI switch, rejected switch continuity, retained model identity and cleared source-bound prefix")


if __name__ == "__main__":
    main()
