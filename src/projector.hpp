#pragma once

#include <ggml-backend.h>
#include <cstring>

namespace gadget::detail {
inline bool initialize_projector_scratch(ggml_tensor *tensor, bool ask, void *) {
    float parameters[2]{};
    std::memcpy(parameters, tensor->op_params, sizeof(parameters));
    const auto *source = tensor->src[0];
    const bool scratch = tensor->op == GGML_OP_SCALE &&
        parameters[0] == 0.0f && parameters[1] == 0.0f && source &&
        source->op == GGML_OP_NONE && !source->view_src &&
        !(source->flags & (GGML_TENSOR_FLAG_INPUT | GGML_TENSOR_FLAG_PARAM)) &&
        source->buffer &&
        ggml_backend_buffer_get_usage(source->buffer) == GGML_BACKEND_BUFFER_USAGE_COMPUTE;
    if (ask)
        return scratch;
    // A projector clears scratch with a zero scale. Multiplication leaves NaNs
    // unchanged, so initialize its output before the following graph operation.
    if (scratch)
        ggml_backend_tensor_memset(tensor, 0, 0, ggml_nbytes(tensor));
    return true;
}
}
