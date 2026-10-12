#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("usage: native-claims operation.rfcxs claims-request output.rfcxs");
        const auto operation = gadget::NativeRecord::read(argv[1]);
        const auto request = gadget::read_record(argv[2]);
        const auto receipt = gadget::check_certificate_claims_record(operation, request.native());
        receipt.write(argv[3]);
        std::cout << gadget::render_record_view(receipt.at("result")).dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
