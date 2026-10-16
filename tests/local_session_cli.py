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
    def run_session(items, cache_mib, name):
        completed = subprocess.run(
            [binary, "session", "--study", str(fixture / "study.json"),
             "--model", model, "--context-tokens", "8192",
             "--prompt-cache-mib", str(cache_mib)],
            input="".join(json.dumps(r) + "\n" for r in items).encode("utf-8"),
            capture_output=True, timeout=150, check=False)
        (root / f"{name}.stdout.jsonl").write_bytes(completed.stdout)
        (root / f"{name}.stderr.log").write_bytes(completed.stderr)
        completed.check_returncode()
        return [json.loads(line) for line in completed.stdout.decode("utf-8").splitlines()]

    cache_mib = 512
    responses = run_session(requests, cache_mib, "session")
    assert len(responses) == len(requests)
    first, failed, preserved, bound, switched, ambiguous = responses
    assert first["status"] == preserved["status"] == switched["status"] == "executed", responses
    assert failed["status"] == ambiguous["status"] == "rejected"
    assert first["source_digest"] == preserved["source_digest"] != switched["source_digest"]
    assert first["model"] == preserved["model"] == switched["model"]
    for response in (first, preserved, switched):
        assert 0 < response["performance"]["prompt_cache_bytes"] <= cache_mib * 1024 * 1024
    assert preserved["performance"]["prefix_tokens_reused"] > 0
    assert switched["performance"]["prefix_tokens_reused"] == 0
    assert bound["schema"] == "gadget.study.bound.v1" and bound["status"] == "bound"
    assert bound["source_digest"] == switched["source_digest"]
    bounded = run_session([prompt, prompt], 1, "bounded")
    assert len(bounded) == 2
    for response in bounded:
        assert response["status"] == "executed", bounded
        assert response["source_digest"] == first["source_digest"]
        assert response["performance"]["prompt_cache_bytes"] == 0
        assert response["performance"]["prefix_tokens_reused"] == 0
    print("PASS: CLI switch, rejected switch continuity, retained model identity, "
          "cleared source-bound prefix and bounded cache admission")


if __name__ == "__main__":
    main()
