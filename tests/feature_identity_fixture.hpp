#pragma once
#include <gadget/native_record.hpp>

namespace gadget::test {
inline NativeRecord repeated_feature_identity_source() {
    using N=NativeRecordValue;
    const std::string label(4096,'X');
    N::Array features;
    for(std::size_t i=0;i<512;++i)
        features.elements.emplace_back(N::Object{{{"id",N(N::ByteSpan{1,0,label.size(),true})}}});
    return NativeRecord::create(N(N::Object{{{"features",N(std::move(features))}}}),
        new_native_record_namespace(),{RawByteSource(1,std::as_bytes(std::span(label)))});
}
}
