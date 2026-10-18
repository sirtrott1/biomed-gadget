#include "../src/projector.hpp"
#include <ggml-alloc.h>
#include <ggml-cpu.h>
#include <array>
#include <iostream>
#include <cmath>

int main() {
    auto *context = ggml_init({ggml_tensor_overhead() * 16 + ggml_graph_overhead(), nullptr, true});
    auto *backend = ggml_backend_cpu_init();
    if (!context || !backend)
        return 1;
    auto *scratch = ggml_new_tensor_1d(context, GGML_TYPE_F32, 4);
    auto *bias = ggml_new_tensor_1d(context, GGML_TYPE_F32, 4);
    ggml_set_input(bias);
    auto *cleared = ggml_scale(context, scratch, 0.0f);
    auto *result = ggml_add(context, cleared, bias);
    auto *graph = ggml_new_graph(context);
    ggml_build_forward_expand(graph, result);
    auto *buffer = ggml_backend_alloc_ctx_tensors(context, backend);
    if (!buffer)
        return 1;
    ggml_backend_buffer_set_usage(buffer, GGML_BACKEND_BUFFER_USAGE_COMPUTE);
    ggml_backend_tensor_memset(scratch, 255, 0, ggml_nbytes(scratch));
    const std::array<float, 4> expected{0.25f, -0.5f, 1.0f, 0.0f};
    ggml_backend_tensor_set(bias, expected.data(), 0, sizeof(expected));
    bool good = gadget::detail::initialize_projector_scratch(cleared, true, nullptr);
    good = good && !gadget::detail::initialize_projector_scratch(result, true, nullptr);
    ggml_set_input(scratch);
    good = good && !gadget::detail::initialize_projector_scratch(cleared, true, nullptr);
    scratch->flags = GGML_TENSOR_FLAG_PARAM;
    good = good && !gadget::detail::initialize_projector_scratch(cleared, true, nullptr);
    scratch->flags = 0;
    auto *scaled = ggml_scale(context, scratch, 2.0f);
    auto *offset = ggml_scale_bias(context, scratch, 0.0f, 2.0f);
    good = good && !gadget::detail::initialize_projector_scratch(scaled, true, nullptr);
    good = good && !gadget::detail::initialize_projector_scratch(offset, true, nullptr);
    auto *scheduler = ggml_backend_sched_new(&backend, nullptr, 1, GGML_DEFAULT_GRAPH_SIZE, false, true);
    good = good && ggml_backend_sched_graph_compute(scheduler, graph) == GGML_STATUS_SUCCESS;
    std::array<float, 4> actual{};
    ggml_backend_tensor_get(result, actual.data(), 0, sizeof(actual));
    for (const auto value : actual)
        good = good && std::isnan(value);
    ggml_backend_sched_set_eval_callback(scheduler, gadget::detail::initialize_projector_scratch, nullptr);
    good = good && ggml_backend_sched_graph_compute(scheduler, graph) == GGML_STATUS_SUCCESS;
    ggml_backend_tensor_get(result, actual.data(), 0, sizeof(actual));
    good = good && actual == expected;
    ggml_backend_sched_free(scheduler);
    ggml_backend_buffer_free(buffer);
    ggml_backend_free(backend);
    ggml_free(context);
    if (!good) {
        std::cerr << "projector scratch initialization failed\n";
        return 1;
    }
}
