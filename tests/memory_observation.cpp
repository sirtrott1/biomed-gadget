#include <rfc/runtime/hardware_inventory.hpp>
#include "../src/memory_observation.hpp"
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace gadget::hardware_detail;
namespace {
void check(bool value,const char *message) { if(!value)throw std::runtime_error(message); }
template<class E,class F>void rejects(F action,const char *message) {
    bool refused=false;try { action(); }catch(const E &){refused=true;}check(refused,message);
}
}
int main(int argc,char **argv) try {
    if(argc!=2)throw std::invalid_argument("owned test directory required");
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    std::map<std::string,std::string> files{
        {"/proc/meminfo","MemTotal: 16 kB\nMemAvailable: 8 kB\n"},
        {"/proc/self/cgroup","0::/child\n"},
        {"/proc/self/mountinfo","1 0 0:1 / /scope rw - cgroup2 cgroup rw\n"},
        {"/scope/child/memory.max","4000\n"},{"/scope/child/memory.current","1000\n"},
        {"/scope/memory.max","6000\n"},{"/scope/memory.current","5000\n"},
        {"/scope/cpu.max","max 100000\n"},{"/scope/child/cpu.max","200000 100000\n"}};
    std::size_t excluded_reads{},reads{};rfc::runtime::LinuxHardwareProbe probe;
    probe.read_file=[&](const std::string &path)->std::optional<std::string> {
        ++reads;if(!memory_path(path))++excluded_reads;
        const auto it=files.find(path);return it==files.end()?std::optional<std::string>{}:it->second;
    };
    auto compare=[&](std::optional<std::size_t> expected) {
        const auto full=rfc::runtime::probe_linux_hardware(probe).advisory_memory_headroom();
        excluded_reads=0;const auto memory=memory_inventory(probe);
        check(!excluded_reads,"memory observation read CPU or PCI sources");
        check(!memory.affinity_observed && memory.pci_devices.empty(),"memory observation invented CPU or device facts");
        check(full==expected && memory.advisory_memory_headroom()==expected,"canonical memory ancestor result differs");
    };
    compare(1000);
    files["/scope/memory.current"]="6000\n";compare(0);
    files["/scope/memory.current"]="7000\n";compare(0);
    files["/scope/memory.max"]="max\n";compare(3000);
    files["/scope/child/memory.max"]="max\n";compare(8192);
    files.erase("/proc/meminfo");compare({});
    files["/scope/child/memory.max"]="4000\n";compare(3000);
    files["/scope/child/memory.current"]="unavailable\n";compare({});
    files["/proc/meminfo"]="MemAvailable: 8 kB\n";compare(8192);
    files["/proc/meminfo"]="MemAvailable: forged kB\n";compare({});
    rejects<std::invalid_argument>([&]{(void)memory_inventory(rfc::runtime::LinuxHardwareProbe{});},
                                   "missing observation reader was not refused");
    auto bounded=probe;bounded.limits.max_file_bytes=8;
    rejects<std::length_error>([&]{(void)memory_inventory(bounded);},"memory observation bypassed canonical input byte limit");
    auto invalid=probe;invalid.limits.max_cpus=0;
    rejects<std::invalid_argument>([&]{(void)memory_inventory(invalid);},"memory observation bypassed canonical probe policy");
    const auto raw=root/"bounded-observation.bin";
    {std::ofstream file(raw,std::ios::binary);file<<"12345678";}
    check(read_observation(raw.string(),8)=="12345678","bounded original observation bytes changed");
    rejects<std::length_error>([&]{(void)read_observation(raw.string(),7);},"oversized actual observation was not refused");
    check(!read_observation((root/"absent").string(),8),"missing observation became supported empty");
    check(reads>10,"fresh observation reader was cached");
    std::cout<<"PASS canonical memory observations, live changes, ancestor bounds, zero/missing/malformed data and bounded original bytes\n";
    return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
