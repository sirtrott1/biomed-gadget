#include <rfc/storage/file_io.hpp>
#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string_view>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

int main(int argc, char **argv) {
    try {
        if (argc != 5 || std::string_view(argv[1]) != "--out" ||
            std::string_view(argv[3]) != "--bytes")
            throw std::invalid_argument("expected --out PATH --bytes N and exactly N stdin bytes");
        const std::string_view count(argv[4]);
        std::size_t bytes{};
        const auto parsed = std::from_chars(count.data(), count.data() + count.size(), bytes);
        if (count.empty() || parsed.ec != std::errc{} || parsed.ptr != count.data() + count.size())
            throw std::invalid_argument("stdin byte count exceeds its unsigned native representation");
#ifdef _WIN32
        if (_setmode(_fileno(stdin), _O_BINARY) == -1)
            throw std::runtime_error("cannot select binary stdin");
#endif
        rfc::storage::publish_file_exact(argv[2], rfc::storage::FilePublication::CreateNew,
                                        [&](std::ostream &out) {
            std::array<std::byte, 1u << 20> buffer{};
            while (bytes) {
                const auto size = std::min(bytes, buffer.size());
                std::cin.read(reinterpret_cast<char *>(buffer.data()), static_cast<std::streamsize>(size));
                if (!std::cin)
                    throw std::runtime_error("stdin truncated or unreadable before declared byte count");
                rfc::storage::write_stream_exact(out, std::span(buffer).first(size));
                bytes -= size;
            }
            if (std::cin.peek() != std::char_traits<char>::eof() || std::cin.bad())
                throw std::runtime_error("stdin has trailing bytes or is unreadable at EOF");
        });
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
