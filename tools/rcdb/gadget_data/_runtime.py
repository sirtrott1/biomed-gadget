"""Freeze the installed application and optional native binding implementations."""
from hashlib import sha256
from importlib.util import find_spec
from pathlib import Path


class ProgramRuntime:
    def __init__(self):
        from rcql.program_codec import ImplementationSnapshot
        self._language = ImplementationSnapshot()
        self._roots = [("gadget_data", Path(__file__).parent)]
        spec = find_spec("rexgraph_cpp")
        if spec is not None and spec.origin is not None:
            self._roots.append(("rexgraph_cpp", Path(spec.origin).parent))
        files = self._files()
        self._stamp = self._seal(files)
        digest = sha256(self._language.check().encode())
        for name, path in files:
            digest.update(name.encode())
            with path.open("rb") as stream:
                digest.update(__import__("hashlib").file_digest(stream, "sha256").digest())
        self.identity = digest.hexdigest()
        self.check()

    def _files(self):
        return tuple((name + "/" + path.name, path) for name, root in self._roots
                     for path in sorted(root.iterdir()) if path.suffix in (".py", ".so", ".pyd"))

    @staticmethod
    def _seal(files):
        result = []
        for name, path in files:
            info = path.stat()
            result.append((name, info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns))
        return tuple(result)

    def check(self):
        self._language.check()
        if self._seal(self._files()) != self._stamp:
            raise ValueError("program application or native binding changed during execution")
        return self.identity
