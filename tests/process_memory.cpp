#include <gadget/process.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#if defined(__linux__)
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
int main() {
    try {
        bool refused = false;
        try { (void)gadget::limit_process_memory(0); }
        catch (const std::invalid_argument &) { refused = true; }
        if (!refused) throw std::runtime_error("zero memory budget admitted");
#if defined(__linux__)
        struct rlimit before{};
        if (getrlimit(RLIMIT_AS, &before)) throw std::runtime_error("getrlimit failed");
        const auto child = fork();
        if (child < 0) throw std::runtime_error("fork failed");
        if (!child) {
            try {
                const auto granted = gadget::limit_process_memory(16 * 1024 * 1024);
                if (!granted || *granted > 16 * 1024 * 1024) _exit(4);
                void *(*volatile allocate)(std::size_t) = std::malloc;
                void *memory = allocate(64 * 1024 * 1024);
                if (memory) { std::free(memory); _exit(5); }
                struct rlimit first{}, second{};
                if (getrlimit(RLIMIT_AS, &first)) _exit(6);
                (void)gadget::limit_process_memory(32 * 1024 * 1024);
                if (getrlimit(RLIMIT_AS, &second) || second.rlim_cur != first.rlim_cur) _exit(7);
                _exit(0);
            } catch (...) { _exit(8); }
        }
        int status{};
        if (waitpid(child, &status, 0) != child || !WIFEXITED(status) || WEXITSTATUS(status))
            throw std::runtime_error("bounded allocation contract child status: " + std::to_string(status));
        struct rlimit after{};
        if (getrlimit(RLIMIT_AS, &after) || before.rlim_cur != after.rlim_cur ||
            before.rlim_max != after.rlim_max)
            throw std::runtime_error("child memory budget changed the caller process");
#endif
        std::cout << "Process memory contracts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
