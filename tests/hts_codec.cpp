#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
    try {
        if (argc < 4)
            throw std::invalid_argument("codec kind source destination [reference]");
        const std::string kind(argv[1]);
        auto *input = hts_open(argv[2], "r");
        auto *output = hts_open(argv[3], kind == "bam" || kind == "bcf" ? "wb"
                                         : kind == "vcfgz"              ? "wz"
                                                                        : "wc");
        if (!input || !output)
            throw std::runtime_error("codec open");
        int status;
        if (kind == "bam" || kind == "cram") {
            auto *header = sam_hdr_read(input);
            auto *record = bam_init1();
            if (kind == "cram" && (argc != 5 || fai_build(argv[4]) != 0 ||
                                   hts_set_fai_filename(output, argv[4]) != 0))
                throw std::runtime_error("CRAM reference");
            if (!header || !record || sam_hdr_write(output, header) < 0)
                throw std::runtime_error("SAM header");
            while ((status = sam_read1(input, header, record)) >= 0)
                if (sam_write1(output, header, record) < 0)
                    throw std::runtime_error("SAM write");
            sam_hdr_destroy(header);
            bam_destroy1(record);
        } else {
            auto *header = bcf_hdr_read(input);
            auto *record = bcf_init();
            if (!header || !record || bcf_hdr_write(output, header) < 0)
                throw std::runtime_error("VCF header");
            while ((status = bcf_read(input, header, record)) == 0)
                if (bcf_write(output, header, record) < 0)
                    throw std::runtime_error("VCF write");
            bcf_hdr_destroy(header);
            bcf_destroy(record);
        }
        if (status != -1 || hts_close(input) != 0 || hts_close(output) != 0)
            throw std::runtime_error("codec completion");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
