#pragma once

// Philox4x32-10 counter-based RNG (Random123 / Salmon, Moraes, Dror, Shaw 2011).
// Official known-answer vectors: DEShawResearch/random123 tests/kat_vectors
//   "philox4x32 10 ..." (not invented).  Copied from the issue-9 CPU engine
// and extended with a stateless (seed, replica, site, stream) uniform.

#include <array>
#include <cstdint>

namespace matching {

inline constexpr uint32_t PHILOX_M4x32_0 = 0xD2511F53u;
inline constexpr uint32_t PHILOX_M4x32_1 = 0xCD9E8D57u;
inline constexpr uint32_t PHILOX_W32_0 = 0x9E3779B9u;
inline constexpr uint32_t PHILOX_W32_1 = 0xBB67AE85u;

inline uint32_t philox_mulhilo32(uint32_t a, uint32_t b, uint32_t& hi) {
    const uint64_t product = static_cast<uint64_t>(a) * static_cast<uint64_t>(b);
    hi = static_cast<uint32_t>(product >> 32);
    return static_cast<uint32_t>(product);
}

inline std::array<uint32_t, 4> philox4x32_10(std::array<uint32_t, 4> ctr,
                                             std::array<uint32_t, 2> key) {
    for (int round = 0; round < 10; ++round) {
        uint32_t hi0 = 0;
        uint32_t hi1 = 0;
        const uint32_t lo0 = philox_mulhilo32(PHILOX_M4x32_0, ctr[0], hi0);
        const uint32_t lo1 = philox_mulhilo32(PHILOX_M4x32_1, ctr[2], hi1);
        const std::array<uint32_t, 4> out = {
            hi1 ^ key[0] ^ ctr[1],
            lo1,
            hi0 ^ key[1] ^ ctr[3],
            lo0,
        };
        ctr = out;
        if (round < 9) {
            key[0] += PHILOX_W32_0;
            key[1] += PHILOX_W32_1;
        }
    }
    return ctr;
}

struct PhiloxKAT {
    std::array<uint32_t, 4> ctr;
    std::array<uint32_t, 2> key;
    std::array<uint32_t, 4> expected;
};

inline std::array<PhiloxKAT, 3> philox4x32_10_official_kats() {
    return {{
        PhiloxKAT{{0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u},
                  {0x00000000u, 0x00000000u},
                  {0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u}},
        PhiloxKAT{{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu},
                  {0xffffffffu, 0xffffffffu},
                  {0x408f276du, 0x41c83b0eu, 0xa20bc7c6u, 0x6d5451fdu}},
        PhiloxKAT{{0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u},
                  {0xa4093822u, 0x299f31d0u},
                  {0xd16cfe09u, 0x94fdccebu, 0x5001e420u, 0x24126ea1u}},
    }};
}

// Stateless uniform in [0, 1).  Independent of thread count and call order.
//   key = (seed_lo, seed_hi)
//   ctr = (index, replica_lo, replica_hi, stream)
// stream 0: site occupation index = cyclic vertex
// stream 1: bond occupation index = packed (src, dx, dy)
inline double counter_uniform(std::uint64_t seed, std::uint64_t replica,
                              std::uint32_t index, std::uint32_t stream = 0) {
    const std::array<uint32_t, 4> ctr{
        index,
        static_cast<uint32_t>(replica),
        static_cast<uint32_t>(replica >> 32),
        stream,
    };
    const std::array<uint32_t, 2> key{
        static_cast<uint32_t>(seed),
        static_cast<uint32_t>(seed >> 32),
    };
    const auto out = philox4x32_10(ctr, key);
    const std::uint64_t bits =
        (static_cast<std::uint64_t>(out[0]) << 32) | static_cast<std::uint64_t>(out[1]);
    return static_cast<double>(bits >> 11) * 0x1.0p-53;
}

inline bool philox_kats_pass() {
    for (const auto& kat : philox4x32_10_official_kats()) {
        if (philox4x32_10(kat.ctr, kat.key) != kat.expected) return false;
    }
    const double first = counter_uniform(17, 23, 5, 0);
    (void)counter_uniform(17, 999, 7, 0);
    if (first != counter_uniform(17, 23, 5, 0) || !(first >= 0.0 && first < 1.0)) {
        return false;
    }
    return true;
}

}  // namespace matching
