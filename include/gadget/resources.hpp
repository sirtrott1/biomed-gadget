#pragma once
#include <functional>
#include <gadget/io.hpp>
#include <rfc/runtime/resource_budget.hpp>
#if defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace gadget {
// Shared admission for cooperating application/model calls. The coordinator,
// observations and reservations are canonical RexGraph objects.
class CpuResources {
    rfc::runtime::HardwareInventory initial_;
    std::function<rfc::runtime::HardwareInventory()> refresh_;
    rfc::runtime::CpuResourceCoordinator coordinator_;
    CpuResources(const CpuResources &parent, rfc::runtime::CpuResourceCoordinator assigned);
    void process() const;
#if defined(__linux__) || defined(__APPLE__)
    const pid_t process_id_{::getpid()};
#endif

  public:
    CpuResources(rfc::runtime::HardwareInventory hardware, rfc::runtime::CpuResourcePolicy policy,
                 std::function<rfc::runtime::HardwareInventory()> refresh);
    static std::shared_ptr<CpuResources> local();
    // Share an authentic P6 assigned coordinator, retaining its parent charge.
    std::shared_ptr<CpuResources> assigned(rfc::runtime::CpuResourceCoordinator &coordinator) const;
    rfc::runtime::HardwareInventory readmit() const;
    rfc::runtime::CpuResourceCoordinator &coordinator() {
        process();
        return coordinator_;
    }
    rfc::runtime::CpuResourceCoordinator::Reservation reserve(std::size_t workers);
    Json report() const;
};
} // namespace gadget
