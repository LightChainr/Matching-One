#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

struct LeftTerm {
    std::int64_t value;
    std::int8_t a0;
    std::int8_t a1;
    std::int8_t a2;
};

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "usage: degree4_fixed_point_screen w1 w2 w3 w4 root_bound near_bound\n";
        return 2;
    }
    std::array<std::int64_t, 4> w{};
    for (int i = 0; i < 4; ++i) w[i] = std::strtoll(argv[i + 1], nullptr, 10);
    const auto root_bound = std::strtoll(argv[5], nullptr, 10);
    const auto near_bound = std::strtoll(argv[6], nullptr, 10);
    constexpr std::int64_t scale = 1000000000000000LL;

    std::vector<LeftTerm> left;
    left.reserve(201ULL * 201ULL * 201ULL);
    for (int a2 = -100; a2 <= 100; ++a2) {
        for (int a1 = -100; a1 <= 100; ++a1) {
            const auto partial = static_cast<std::int64_t>(a1) * w[0] + static_cast<std::int64_t>(a2) * w[1];
            for (int a0 = -100; a0 <= 100; ++a0) {
                left.push_back({partial + static_cast<std::int64_t>(a0) * scale,
                                static_cast<std::int8_t>(a0), static_cast<std::int8_t>(a1),
                                static_cast<std::int8_t>(a2)});
            }
        }
    }
    std::sort(left.begin(), left.end(), [](const LeftTerm& x, const LeftTerm& y) {
        if (x.value != y.value) return x.value < y.value;
        if (x.a0 != y.a0) return x.a0 < y.a0;
        if (x.a1 != y.a1) return x.a1 < y.a1;
        return x.a2 < y.a2;
    });

    std::uint64_t near_count = 0;
    std::uint64_t root_filter_count = 0;
    for (int a4 = 1; a4 <= 100; ++a4) {
        for (int a3 = -100; a3 <= 100; ++a3) {
            const auto right = static_cast<std::int64_t>(a3) * w[2] + static_cast<std::int64_t>(a4) * w[3];
            const auto lo = -right - near_bound;
            const auto hi = -right + near_bound;
            auto first = std::lower_bound(left.begin(), left.end(), lo,
                                          [](const LeftTerm& x, std::int64_t value) { return x.value < value; });
            auto last = std::upper_bound(left.begin(), left.end(), hi,
                                         [](std::int64_t value, const LeftTerm& x) { return value < x.value; });
            for (auto it = first; it != last; ++it) {
                const int a0 = it->a0, a1 = it->a1, a2 = it->a2;
                if (std::gcd(std::gcd(std::gcd(std::abs(a0), std::abs(a1)), std::abs(a2)),
                             std::gcd(std::abs(a3), a4)) != 1) {
                    continue;
                }
                const auto residual = it->value + right;
                const bool root_filter = std::llabs(residual) <= root_bound;
                ++near_count;
                root_filter_count += root_filter;
                std::cout << a0 << ' ' << a1 << ' ' << a2 << ' ' << a3 << ' ' << a4 << ' '
                          << residual << ' ' << (root_filter ? 1 : 0) << '\n';
            }
        }
    }
    std::cerr << "near_candidates=" << near_count << " root_filter_candidates=" << root_filter_count << '\n';
    return 0;
}
