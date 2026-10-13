"""Installed Agent record retrieval and native admission refusal over loopback."""
import argparse
from contextlib import closing
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time


def adapter_module():
    path = Path(__file__).resolve().parents[1] / "tools/rcdb/gadget_rcdb.py"
    spec = importlib.util.spec_from_file_location("gadget_rcdb", path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    return adapter


def configure(root):
    os.environ.update(REXGRAPH_CONFIG_DIR=str(root / "config"),
                      REXGRAPH_ACTIVITY_JOURNAL=str(root / "activity.jsonl"),
                      REXGRAPH_AUDIT_JOURNAL=str(root / "audit.jsonl"),
                      REXGRAPH_RCDB_URI="local://" + str(root / "remote.rcdb"),
                      REXGRAPH_FRAME_KEY="owned native record contract key")


def serve(args):
    configure(args.root / "server")
    from fastapi import FastAPI
    import uvicorn
    from rcdb import LocalStore
    from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
    from agent import activity
    from agent.server import audit, auth, budget, scope
    from agent.server.routes import rex
    from agent.server.security import add_auth_enforcement
    adapter = adapter_module()
    raw = args.publication.read_bytes()
    attempts = args.root / "native-admissions.txt"
    original = adapter._validate_publication_file

    def observed(*values, **keywords):
        with attempts.open("ab") as file:
            file.write(b"x")
        return original(*values, **keywords)

    adapter._validate_publication_file = observed
    codec = adapter.install_codec(gadget=args.binary, max_bytes=len(raw))
    destination = args.root / "server/remote.rcdb"
    if not destination.exists():
        with closing(LocalStore(destination, compression=BlobCodecSpec(name=args.compression),
                                record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
            for version, workspace in [(1, "alpha"), (2, "beta")]:
                record = store.put_record("literal@2/β?one#%", raw, codec=codec,
                                          expected_version=version - 1,
                                          tx_time=float(version * 10), valid_from=100., valid_to=200.,
                                          meta={"workspace": workspace, "weight": Fraction(2, 7),
                                                "exact": 2**1000, "bytes": b"original\0metadata"},
                                          tags=["native", "β"])
                assert record.version == version
            store.checkpoint()
    manager = auth.get_auth_manager()
    manager.enable_auth()
    ready = {name: manager.create_token(name, [name], role="user") for name in ("alpha", "beta")}
    listener = socket.socket(fileno=args.serve)
    ready["port"] = listener.getsockname()[1]
    app = FastAPI()
    app.include_router(rex.router)
    budget.add_compute_budget(app)
    scope.add_workspace_scope(app)
    add_auth_enforcement(app)
    (args.root / "server-ready.json").write_text(json.dumps(ready))
    server = uvicorn.Server(uvicorn.Config(app, log_level="warning", access_log=False))
    # Uvicorn restores and replays the caller's signal handler after shutdown.
    signal.signal(signal.SIGTERM, lambda *_: None)
    try:
        server.run(sockets=[listener])
    finally:
        from rcdb import reset_default_store
        reset_default_store()
        activity.get_log().close()
        verification = audit.verify()
        (args.root / "server-final.json").write_text(json.dumps(
            {"inflight": budget.inflight(), "audit": verification}))


def run(args):
    args.root.mkdir(parents=True)
    configure(args.root / "client")
    from agent.client import RexClient
    from rcdb import LocalStore, RecordPacket, VALUE_CODEC
    from rcdb.codecs import record_codec
    from rcdb.header import DEFAULT_RECORD_CODECS
    from agent.server.scope import ScopedStore
    import httpx
    adapter = adapter_module()
    raw = args.publication.read_bytes()
    source_digest = hashlib.sha256(raw).hexdigest()
    codec = adapter.install_codec(gadget=args.binary, max_bytes=len(raw))
    try:
        record_codec(codec).cell_counts(raw)
    except ValueError as error:
        assert str(error) == "record codec has no declared cell admission capability"
    else:
        raise AssertionError("native codec invented a cell admission capability")
    with closing(LocalStore(args.root / "local.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as local:
        alpha, beta = ScopedStore(local, "alpha", "alice"), ScopedStore(local, "beta", "bob")
        alpha.put_record("owned", raw, codec=codec, expected_version=0)
        assert alpha.read_record("owned").value == raw and beta.read_record("owned") is None
        try:
            beta.put_record("owned", raw, codec=codec, expected_version=1)
        except PermissionError:
            pass
        else:
            raise AssertionError("another workspace changed a native publication")
        assert local.read_record("owned").record.version == 1
    results = []
    process = None

    def admissions():
        path = args.root / "native-admissions.txt"
        return path.stat().st_size if path.exists() else 0

    def start(generation):
        ready_path = args.root / "server-ready.json"
        ready_path.unlink(missing_ok=True)
        listener = socket.socket()
        listener.bind(("127.0.0.1", 0))
        listener.listen()
        command = [sys.executable, "-I", "-B", str(Path(__file__).resolve()),
                   str(args.binary), str(args.publication), str(args.root),
                   "--compression", args.compression, "--serve", str(listener.fileno())]
        with (args.root / f"server-{generation}.log").open("w") as log:
            child = subprocess.Popen(command, pass_fds=(listener.fileno(),), stdout=log,
                                     stderr=subprocess.STDOUT, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        listener.close()
        try:
            deadline = time.monotonic() + 60
            while time.monotonic() < deadline:
                if child.poll() is not None:
                    raise RuntimeError((args.root / f"server-{generation}.log").read_text())
                if ready_path.exists():
                    ready = json.loads(ready_path.read_text())
                    client = RexClient(f"http://127.0.0.1:{ready['port']}", api_key=ready["alpha"], workspace="alpha")
                    try:
                        hello = client.rex_hello()
                    except httpx.TransportError:
                        pass
                    else:
                        assert hello["auth_enabled"] and hello["signed_frames"]
                        assert hello["workspace"] == "alpha"
                        return child, ready, client
                time.sleep(.05)
            raise RuntimeError("owned Agent listener did not become ready within 60 seconds")
        except BaseException:
            child.terminate()
            child.wait(timeout=30)
            raise

    def stop(child):
        child.terminate()
        assert child.wait(timeout=30) == 0
        final = json.loads((args.root / "server-final.json").read_text())
        assert not final["inflight"] and final["audit"]["valid"]

    def refused(action, status, text=None):
        before = admissions()
        try:
            action()
        except httpx.HTTPStatusError as error:
            assert error.response.status_code == status, error.response.text
            if text:
                assert text in error.response.text
        else:
            raise AssertionError("Agent record request did not refuse")
        assert admissions() == before, "refusal invoked native decoding"
        results.append({"status": status, "native_decode_before_refusal": False})

    def checked(packet, version, workspace):
        snapshot = packet.snapshot()
        assert snapshot.value == raw and snapshot.state_digest == packet.state_digest
        record = snapshot.record
        assert record.id == "literal@2/β?one#%" and record.version == version
        assert record.meta == {"workspace": workspace, "weight": Fraction(2, 7),
                               "exact": 2**1000, "bytes": b"original\0metadata"}
        assert record.signature["tags"] == ["native", "β"]
        assert (record.valid_from, record.valid_to) == (100., 200.)
        assert record.envelope.codec == codec.name
        return snapshot

    try:
        process, ready, alpha = start(1)
        beta = RexClient(alpha.url, api_key=ready["beta"], workspace="beta")
        packet = alpha.rex_fetch_record("literal@2/β?one#%", version=1)
        checked(packet, 1, "alpha")
        checked(alpha.rex_fetch_record("literal@2/β?one#%", as_of=15., valid_at=150.), 1, "alpha")
        checked(beta.rex_fetch_record("literal@2/β?one#%"), 2, "beta")
        for client, selectors in [(alpha, {}), (beta, {"version": 1}),
                                  (alpha, {"version": 3}), (alpha, {"as_of": 25.}),
                                  (alpha, {"as_of": 15., "valid_at": 200.})]:
            refused(lambda: client.rex_fetch_record("literal@2/β?one#%", **selectors), 404)
        before = admissions()
        try:
            alpha.rex_fetch_record("literal@2/β?one#%", version=1, valid_at=150.)
        except ValueError as error:
            assert "exact version cannot be combined with time selectors" in str(error)
        else:
            raise AssertionError("incompatible native record selectors were admitted")
        assert admissions() == before
        foreign = RexClient(alpha.url, api_key=ready["alpha"], workspace="beta")
        refused(lambda: foreign.rex_fetch_record("literal@2/β?one#%", version=1), 403)
        anonymous = RexClient(alpha.url, workspace="alpha")
        refused(lambda: anonymous.rex_fetch_record("literal@2/β?one#%", version=1), 401)
        with closing(LocalStore(args.root / "server/remote.rcdb")) as remote:
            before = remote.state_digest()
        for _ in range(2):
            refused(lambda: alpha.rex_store_record(packet), 400, "no declared cell admission capability")
        unsigned = RexClient(alpha.url, api_key=ready["alpha"], workspace="alpha")
        unsigned.frame_key = None
        refused(lambda: unsigned.rex_store_record(packet), 401)
        from rcdb.packet import PACKET_CONTENT_TYPE
        from rexgraph.protocol import sign
        corrupt = packet.to_bytes()[:-1]
        headers = alpha._headers()
        headers.update({"Content-Type": PACKET_CONTENT_TYPE,
                        "X-Rex-Signature": sign(corrupt, alpha.frame_key)})
        refused(lambda: httpx.post(alpha.url + "/rex/v1/records/store", content=corrupt,
                                  headers=headers, timeout=300).raise_for_status(), 400)
        with closing(LocalStore(args.root / "server/remote.rcdb")) as remote:
            assert remote.state_digest() == before and remote.stats()["n_versions"] == 2
        wrong_key = RexClient(alpha.url, api_key=ready["alpha"], workspace="alpha", frame_key=b"different owned test key")
        try:
            wrong_key.rex_fetch_record("literal@2/β?one#%", version=1)
        except ValueError as error:
            assert "signature" in str(error)
        else:
            raise AssertionError("client returned a packet with a different frame key")
        with closing(LocalStore(args.root / "control.rcdb")) as control:
            control.put_record("declared-value", Fraction(1, 7), codec=VALUE_CODEC, expected_version=0)
            value = RecordPacket.from_snapshot(control.read_record("declared-value"), source_store_id=control.store_id)
            receipt = alpha.rex_store_record(value)
            assert receipt.source_record_id == "declared-value" and receipt.source_version == 1
            assert receipt.source_digest == receipt.destination_digest == value.state_digest
            restored = alpha.rex_fetch_record(receipt.destination_record_id)
            assert restored.snapshot().value == Fraction(1, 7)
        stop(process)
        process = None
        process, _, alpha = start(2)
        checked(alpha.rex_fetch_record("literal@2/β?one#%", version=1), 1, "alpha")
        assert alpha.rex_fetch_record(receipt.destination_record_id).snapshot().value == Fraction(1, 7)
        received = args.root / "received.rfcxs"
        received.write_bytes(packet.snapshot().value)
        assert adapter.validate_publication(received.read_bytes(), gadget=args.binary) == raw
        stop(process)
        process = None
        assert args.publication.read_bytes() == raw
        report = {"status": "PASS", "source_bytes": len(raw), "source_sha256": source_digest,
                  "compression": args.compression, "transport": "actual owned IPv4 loopback socket, installed canonical RexClient and Agent native routes/auth/scope/budget",
                  "refusals": results, "complete_source_bytes_preserved": True,
                  "literal_versions_times_metadata_workspace_and_restart": True,
                  "declared_value_positive_control": True,
                  "native_upload_cell_admission": "Refused before native decoding and destination publication; canonical count capability remains required.",
                  "scope": "Canonical record routes only; no full Agent application, model, other host, native count provider or broad deployment qualification."}
        (args.root / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")
        print("PASS: installed loopback native record retrieval, literal versions, full bytes, scoped refusal, declared value control and restart")
    finally:
        if process is not None and process.poll() is None:
            process.terminate()
            process.wait(timeout=30)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=lambda value: Path(value).resolve())
    parser.add_argument("publication", type=lambda value: Path(value).resolve())
    parser.add_argument("root", type=lambda value: Path(value).resolve())
    parser.add_argument("--compression", choices=("none", "zstd", "zlib"), default="none")
    parser.add_argument("--serve", type=int)
    arguments = parser.parse_args()
    serve(arguments) if arguments.serve is not None else run(arguments)
