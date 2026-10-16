#include "assay_detail.hpp"
#if GADGET_HAS_HDF5
#include <cmath>
#include <hdf5.h>
#include <mutex>
#endif

namespace gadget::assay_detail {
#if GADGET_HAS_HDF5
namespace {
std::mutex hdf5_mutex;
struct Handle {
    hid_t id;
    herr_t (*close)(hid_t);
    Handle(hid_t h, herr_t (*c)(hid_t)) : id(h), close(c) {
        require(id >= 0, "HDF5 object admission failed");
    }
    ~Handle() { close(id); }
    Handle(const Handle &) = delete;
};
void hard_link(hid_t file, const std::string &path) {
    // Every traversed component is internal to the retained image.
    std::string current;
    for (const auto &part : split(path, '/'))
        if (!part.empty()) {
            current += "/" + part;
            H5L_info2_t info{};
            require(H5Lget_info2(file, current.c_str(), &info, H5P_DEFAULT) >= 0 &&
                        info.type == H5L_TYPE_HARD,
                    "H5AD external/soft links require a separately retained source");
        }
}
std::string attribute(hid_t object, const char *name) {
    Handle a(H5Aopen(object, name, H5P_DEFAULT), H5Aclose), t(H5Aget_type(a.id), H5Tclose);
    require(H5Tget_class(t.id) == H5T_STRING, "H5AD string attribute");
    Handle space(H5Aget_space(a.id), H5Sclose);
    require(H5Sget_simple_extent_npoints(space.id) == 1, "H5AD scalar attribute");
    if (H5Tis_variable_str(t.id) > 0) {
        char *value = nullptr;
        require(H5Aread(a.id, t.id, &value) >= 0 && value, "H5AD attribute read");
        std::string result(value);
        H5free_memory(value);
        return result;
    }
    const auto size = H5Tget_size(t.id);
    require(size > 0 && size <= source_byte_limit(), "H5AD attribute representation");
    std::string result(size, '\0');
    require(H5Aread(a.id, t.id, result.data()) >= 0, "H5AD attribute read");
    const auto zero = result.find('\0');
    if (zero != std::string::npos)
        result.resize(zero);
    return result;
}
void encoding(hid_t object, const char *type, const char *version) {
    require(attribute(object, "encoding-type") == type &&
                attribute(object, "encoding-version") == version,
            "unsupported H5AD encoding type/version");
}
std::vector<hsize_t> shape(hid_t dataset) {
    Handle space(H5Dget_space(dataset), H5Sclose);
    const int rank = H5Sget_simple_extent_ndims(space.id);
    require(rank >= 1 && rank <= 2, "H5AD dataset rank");
    std::vector<hsize_t> dims(rank);
    require(H5Sget_simple_extent_dims(space.id, dims.data(), nullptr) >= 0,
            "H5AD dataset dimensions");
    Handle properties(H5Dget_create_plist(dataset), H5Pclose);
    require(H5Pget_external_count(properties.id) == 0 &&
                H5Pget_layout(properties.id) != H5D_VIRTUAL,
            "H5AD externally stored or virtual datasets require retained dependencies");
    return dims;
}
std::vector<unsigned long long> indices(hid_t file, const std::string &path, std::size_t expected) {
    hard_link(file, path);
    Handle d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
    const auto dims = shape(d.id);
    require(dims.size() == 1 && dims[0] == expected, "H5AD index shape");
    Handle t(H5Dget_type(d.id), H5Tclose);
    require(H5Tget_class(t.id) == H5T_INTEGER && H5Tget_size(t.id) <= 8, "H5AD index dtype");
    std::vector<unsigned long long> out(dims[0]);
    if (H5Tget_sign(t.id) == H5T_SGN_2) {
        std::vector<long long> v(out.size());
        require(H5Dread(d.id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
                "H5AD index read");
        for (std::size_t i = 0; i < v.size(); ++i) {
            require(v[i] >= 0, "negative H5AD sparse index");
            out[i] = v[i];
        }
    } else
        require(H5Dread(d.id, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0,
                "H5AD index read");
    return out;
}
std::vector<std::string> axis(hid_t file, const std::string &path) {
    hard_link(file, path);
    Handle group(H5Gopen2(file, path.c_str(), H5P_DEFAULT), H5Gclose);
    encoding(group.id, "dataframe", "0.2.0");
    const auto index = attribute(group.id, "_index");
    require(!index.empty() && index.find('/') == std::string::npos, "H5AD axis index name");
    hard_link(file, path + "/" + index);
    Handle d(H5Dopen2(group.id, index.c_str(), H5P_DEFAULT), H5Dclose);
    encoding(d.id, "string-array", "0.2.0");
    const auto dims = shape(d.id);
    require(dims.size() == 1, "H5AD axis size");
    (void)rfc::machine::checked_size(dims[0]);
    Handle t(H5Dget_type(d.id), H5Tclose);
    require(H5Tget_class(t.id) == H5T_STRING, "H5AD axis strings");
    std::vector<std::string> out;
    out.reserve(dims[0]);
    if (H5Tis_variable_str(t.id) > 0) {
        std::vector<char *> v(dims[0], nullptr);
        require(H5Dread(d.id, t.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
                "H5AD axis read");
        for (auto s : v)
            out.emplace_back(s ? s : "");
        Handle space(H5Dget_space(d.id), H5Sclose);
        require(H5Dvlen_reclaim(t.id, space.id, H5P_DEFAULT, v.data()) >= 0, "H5AD string release");
    } else {
        const auto width = H5Tget_size(t.id);
        const auto bytes = rfc::machine::checked_size_mul(width, rfc::machine::checked_size(dims[0]));
        require(width && bytes <= source_byte_limit(), "H5AD string representation");
        std::string data(bytes, '\0');
        require(H5Dread(d.id, t.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) >= 0,
                "H5AD axis read");
        for (std::size_t i = 0; i < dims[0]; ++i) {
            auto s = data.substr(i * width, width);
            auto p = s.find('\0');
            if (p != std::string::npos)
                s.resize(p);
            out.push_back(s);
        }
    }
    return out;
}
std::vector<std::optional<Rational>> values(hid_t dataset, const Json &spec) {
    const auto dims = shape(dataset);
    std::size_t count = 1;
    for (auto d : dims) {
        count = rfc::machine::checked_size_mul(count, rfc::machine::checked_size(d));
    }
    Handle t(H5Dget_type(dataset), H5Tclose);
    std::vector<std::optional<Rational>> out(count);
    if (H5Tget_class(t.id) == H5T_INTEGER && H5Tget_size(t.id) <= 8) {
        if (H5Tget_sign(t.id) == H5T_SGN_2) {
            std::vector<long long> v(count);
            require(H5Dread(dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >=
                        0,
                    "H5AD integer read");
            for (std::size_t i = 0; i < count; ++i)
                out[i] = exact_rational(std::to_string(v[i]));
        } else {
            std::vector<unsigned long long> v(count);
            require(H5Dread(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >=
                        0,
                    "H5AD integer read");
            for (std::size_t i = 0; i < count; ++i)
                out[i] = exact_rational(std::to_string(v[i]));
        }
    } else {
        require(H5Tget_class(t.id) == H5T_FLOAT &&
                    (H5Tequal(t.id, H5T_IEEE_F32LE) > 0 || H5Tequal(t.id, H5T_IEEE_F32BE) > 0 ||
                     H5Tequal(t.id, H5T_IEEE_F64LE) > 0 || H5Tequal(t.id, H5T_IEEE_F64BE) > 0),
                "H5AD floating profile requires stored IEEE-754 float32/float64");
        std::vector<double> v(count);
        require(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
                "H5AD floating read");
        for (std::size_t i = 0; i < count; ++i) {
            if (std::isnan(v[i]) && spec.at("nan_policy").text() == "unobserved")
                continue;
            require(std::isfinite(v[i]), "nonfinite H5AD value");
            out[i] = Rational(v[i]);
            out[i]->canonicalize();
        }
    }
    if (is_count_quantity(spec.at("quantity").text()))
        for (const auto &v : out)
            if (v)
                require(*v >= 0 && v->get_den() == 1,
                        "raw H5AD counts require nonnegative integers");
    return out;
}
} // namespace
Matrix read_h5ad(const fs::path &root, const Json &spec) {
    require(spec.at("arithmetic").text() == "ieee754_exact",
            "H5AD arithmetic must declare exact stored IEEE-754 values");
    require(spec.has("nan_policy") && (spec.at("nan_policy").text() == "reject" ||
                                       spec.at("nan_policy").text() == "unobserved"),
            "H5AD NaN policy");
    require(!spec.has("missing") && spec.at("missing_tokens").size() == 0,
            "H5AD missingness is represented by the declared NaN policy");
    Matrix out;
    out.bytes = local_artifact(root, spec.at("source"));
    out.features = read_axis(local_artifact(root, spec.at("features")));
    out.observations = read_axis(local_artifact(root, spec.at("observations")));
    const auto nr = out.features.ids.size(), nc = out.observations.ids.size();
    (void)rfc::numeric::IndexDomain({nr, nc});
    const auto implicit = spec.at("implicit_entries").text();
    require(implicit == "measured_zero" || implicit == "unobserved", "H5AD implicit entries");
    // Open the already hashed bytes, so admission cannot read a different on disk revision.
    std::lock_guard lock(hdf5_mutex);
    Handle access(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
    require(H5Pset_fapl_core(access.id, 65536, false) >= 0 &&
                H5Pset_file_image(access.id, out.bytes.data(), out.bytes.size()) >= 0,
            "H5AD retained image");
    Handle file(H5Fopen("gadget-retained-h5ad-image", H5F_ACC_RDONLY, access.id), H5Fclose);
    encoding(file.id, "anndata", "0.1.0");
    require(axis(file.id, "/obs") == out.observations.ids &&
                axis(file.id, "/var") == out.features.ids,
            "H5AD axis identities/order differ from sidecars");
    hard_link(file.id, "/X");
    H5O_info2_t info{};
    require(H5Oget_info_by_name3(file.id, "/X", &info, H5O_INFO_BASIC, H5P_DEFAULT) >= 0,
            "H5AD X object");
    std::vector<std::vector<rfc::numeric::SparseRationalEntry>> entries(nc);
    std::vector<std::set<std::size_t>> support(nc), missing(nc);
    bool dense = false;
    if (info.type == H5O_TYPE_DATASET) {
        dense = true;
        Handle d(H5Dopen2(file.id, "/X", H5P_DEFAULT), H5Dclose);
        encoding(d.id, "array", "0.2.0");
        require(shape(d.id) == std::vector<hsize_t>{nc, nr}, "H5AD X shape");
        const auto v = values(d.id, spec);
        for (std::size_t c = 0; c < nc; ++c)
            for (std::size_t r = 0; r < nr; ++r) {
                const auto &x = v[c * nr + r];
                if (!x) {
                    missing[c].insert(r);
                    continue;
                }
                support[c].insert(r);
                if (*x != 0)
                    entries[c].push_back({r, *x});
            }
    } else {
        require(info.type == H5O_TYPE_GROUP, "H5AD X storage type");
        Handle group(H5Gopen2(file.id, "/X", H5P_DEFAULT), H5Gclose);
        const auto enc = attribute(group.id, "encoding-type");
        require((enc == "csr_matrix" || enc == "csc_matrix") &&
                    attribute(group.id, "encoding-version") == "0.1.0",
                "H5AD sparse encoding");
        Handle sa(H5Aopen(group.id, "shape", H5P_DEFAULT), H5Aclose),
            ss(H5Aget_space(sa.id), H5Sclose);
        require(H5Sget_simple_extent_npoints(ss.id) == 2, "H5AD sparse shape");
        long long sh[2]{};
        require(H5Aread(sa.id, H5T_NATIVE_LLONG, sh) >= 0 && sh[0] == static_cast<long long>(nc) &&
                    sh[1] == static_cast<long long>(nr),
                "H5AD sparse axis shape");
        hard_link(file.id, "/X/data");
        Handle d(H5Dopen2(file.id, "/X/data", H5P_DEFAULT), H5Dclose);
        require(shape(d.id).size() == 1, "H5AD sparse data rank");
        const auto v = values(d.id, spec);
        const auto major = enc == "csr_matrix" ? nc : nr, minor = enc == "csr_matrix" ? nr : nc;
        const auto ix = indices(file.id, "/X/indices", v.size()),
                   ptr = indices(file.id, "/X/indptr", major + 1);
        require(ptr.size() == major + 1 && ptr[0] == 0 && ptr.back() == v.size() &&
                    ix.size() == v.size(),
                "H5AD sparse pointer shape");
        for (std::size_t j = 0; j < major; ++j) {
            require(ptr[j] <= ptr[j + 1] && ptr[j + 1] <= v.size(), "H5AD sparse pointers");
            std::set<std::size_t> seen;
            for (std::size_t k = ptr[j]; k < ptr[j + 1]; ++k) {
                require(ix[k] < minor && seen.insert(ix[k]).second,
                        "H5AD sparse index range/duplicate");
                const auto c = enc == "csr_matrix" ? j : ix[k], r = enc == "csr_matrix" ? ix[k] : j;
                if (!v[k]) {
                    missing[c].insert(r);
                    continue;
                }
                support[c].insert(r);
                if (*v[k] != 0)
                    entries[c].push_back({r, *v[k]});
            }
        }
    }
    for (std::size_t c = 0; c < nc; ++c) {
        const bool full = missing[c].empty() && (dense || implicit == "measured_zero");
        std::vector<std::size_t> s;
        if (!full)
            for (std::size_t r = 0; r < nr; ++r)
                if (!missing[c].contains(r) &&
                    (dense || implicit == "measured_zero" || support[c].contains(r)))
                    s.push_back(r);
        out.columns.emplace_back(rfc::numeric::IndexDomain::linear(nr),
                                 full ? rfc::numeric::StructuralSupport::full(nr)
                                      : rfc::numeric::StructuralSupport::sparse(nr, s),
                                 rfc::numeric::SparseRationalVector(nr, std::move(entries[c])));
    }
    return out;
}
#else
Matrix read_h5ad(const fs::path &, const Json &) {
    throw std::invalid_argument("H5AD requires a build with GADGET_ENABLE_HDF5=ON");
}
#endif
} // namespace gadget::assay_detail
