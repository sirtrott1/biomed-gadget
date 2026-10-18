#include <gadget/build.hpp>
#include "gadget_build_fingerprints.hpp"
namespace gadget::build {
const char *source_sha256() noexcept { return source_fingerprint; }
const char *rexgraph_sha256() noexcept { return rexgraph_fingerprint; }
const char *backend_sha256() noexcept { return backend_fingerprint; }
}
