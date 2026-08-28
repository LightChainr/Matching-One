#pragma once

// Philox4x32-10 counter-based RNG (Random123 / Salmon, Moraes, Dror, Shaw 2011).
// Official known-answer vectors: DEShawResearch/random123 tests/kat_vectors
//   "philox4x32 10 ..." (not invented).
//
// Application counter layout for production streams:
//   key[0] = global_seed
//   key[1] = batch_id
//   ctr[0] = draw_counter   (increments once per 4-word block)
//   ctr[1] = replica_id
//   ctr[2] = stream_id      (0 = G / shared; 1 = independent G*)
//   ctr[3] = 0
// The same (global_seed, batch_id, replica_id, draw_counter) plus stream_id
// always yields the same block; thread scheduling cannot change values.

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

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

// Random123 tests/kat_vectors: "For each generator, we test: gen(0,0),
// gen(fff,fff) and gen(ctr=digits_of_pi, key=more_digits_of_pi)."
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

struct PhiloxStream {
    uint32_t global_seed = 0;
    uint32_t batch_id = 0;
    uint32_t replica_id = 0;
    uint32_t stream_id = 0;
    uint32_t draw_counter = 0;
    std::array<uint32_t, 4> buf{};
    int buf_i = 4;

    void reset(uint32_t seed, uint32_t batch, uint32_t replica, uint32_t stream) {
        global_seed = seed;
        batch_id = batch;
        replica_id = replica;
        stream_id = stream;
        draw_counter = 0;
        buf_i = 4;
    }

    void refill() {
        const std::array<uint32_t, 4> ctr{draw_counter, replica_id, stream_id, 0u};
        const std::array<uint32_t, 2> key{global_seed, batch_id};
        buf = philox4x32_10(ctr, key);
        buf_i = 0;
        draw_counter += 1u;
    }

    uint32_t next_u32() {
        if (buf_i >= 4) {
            refill();
        }
        return buf[static_cast<std::size_t>(buf_i++)];
    }

    // Lemire nearly-divisionless unbiased sample in 0 .. range-1.
    uint32_t bounded(uint32_t range) {
        if (range <= 1u) {
            return 0u;
        }
        uint32_t x = next_u32();
        uint64_t m = static_cast<uint64_t>(x) * static_cast<uint64_t>(range);
        uint32_t l = static_cast<uint32_t>(m);
        if (l < range) {
            const uint32_t t = (0u - range) % range;
            while (l < t) {
                x = next_u32();
                m = static_cast<uint64_t>(x) * static_cast<uint64_t>(range);
                l = static_cast<uint32_t>(m);
            }
        }
        return static_cast<uint32_t>(m >> 32);
    }
};

inline void fisher_yates(std::vector<int>& perm, PhiloxStream& rng) {
    const int n = static_cast<int>(perm.size());
    for (int i = 0; i < n; ++i) {
        perm[static_cast<std::size_t>(i)] = i;
    }
    for (int i = n - 1; i > 0; --i) {
        const uint32_t j = rng.bounded(static_cast<uint32_t>(i + 1));
        std::swap(perm[static_cast<std::size_t>(i)], perm[static_cast<std::size_t>(j)]);
    }
}

}  // namespace matching
