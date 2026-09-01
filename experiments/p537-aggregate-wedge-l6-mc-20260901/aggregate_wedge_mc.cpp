// Held-out square-L6 Bernoulli production for the P537 root-conditioned G4.
// The scientific coordinate is frozen to L^4*G4; no radius/source/minor scan.
#define P537_LIBRARY_ONLY
#include "../p537-aggregate-wedge-l5-20260901/aggregate_wedge_exact.cpp"

#include <cmath>

namespace {

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

double unit(std::uint64_t x) {
    return static_cast<double>(splitmix64(x) >> 11) * 0x1.0p-53;
}

struct McGlobalRow {
    std::uint64_t count = 0;
    std::int64_t sum_q0 = 0, sum_q1 = 0;
    std::int64_t sum_source16_0 = 0, sum_source16_1 = 0;
};

struct McLandingRow {
    std::int64_t signed_count = 0;
    std::int64_t signed_source_mid16 = 0;
    std::uint64_t unsigned_count = 0;
};

void write_mc(const std::string& path, std::uint64_t samples, int shard_index,
              int shard_count, int batches, std::uint64_t seed, double proposal_p,
              std::uint64_t begin, std::uint64_t end,
              const std::vector<std::vector<McGlobalRow>>& global,
              const std::array<std::vector<std::vector<McLandingRow>>,2>& landing) {
    std::ifstream probe(path);
    if (probe.good()) throw std::runtime_error("refusing to overwrite: " + path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open output: " + path);
    out << "# schema=matching-one/p537-aggregate-wedge-mc/v1\n"
        << "# L=6\n# N=36\n# samples=" << samples
        << "\n# shard_index=" << shard_index << "\n# shard_count=" << shard_count
        << "\n# batches=" << batches << "\n# seed=" << seed
        << "\n# proposal_p=" << std::setprecision(17) << proposal_p
        << "\n# begin=" << begin << "\n# end=" << end << '\n';
    out << "batch\tkind\ttransition\tk\tcount\tsum_q0\tsum_q1\t"
           "sum_source16_0\tsum_source16_1\tsigned_count\t"
           "signed_source_mid16\tunsigned_count\n";
    for (int batch = 0; batch < batches; ++batch)
    for (std::size_t k = 0; k < global[batch].size(); ++k) {
        const auto& g = global[batch][k];
        out << batch << "\tglobal\t-\t" << k << '\t' << g.count << '\t'
            << g.sum_q0 << '\t' << g.sum_q1 << '\t' << g.sum_source16_0 << '\t'
            << g.sum_source16_1 << "\t0\t0\t0\n";
        for (int tr = 0; tr < 2; ++tr) {
            const auto& a = landing[tr][batch][k];
            out << batch << "\tlanding\t" << (tr ? "12" : "01") << '\t' << k
                << "\t0\t0\t0\t0\t0\t" << a.signed_count << '\t'
                << a.signed_source_mid16 << '\t' << a.unsigned_count << '\n';
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 10) {
        std::cerr << "usage: aggregate_wedge_mc KERNEL OUTPUT SAMPLES SHARD_INDEX "
                     "SHARD_COUNT BATCHES SEED PROPOSAL_P RESERVED\n";
        return 2;
    }
    try {
        const auto kernel = read_kernel(argv[1]);
        const std::string output = argv[2];
        const std::uint64_t samples = std::stoull(argv[3]);
        const int shard_index = std::stoi(argv[4]);
        const int shard_count = std::stoi(argv[5]);
        const int batches = std::stoi(argv[6]);
        const std::uint64_t seed = std::stoull(argv[7]);
        const double proposal_p = std::stod(argv[8]);
        const std::string reserved = argv[9];
        if (reserved != "frozen-L6-G4") throw std::invalid_argument("reserved token mismatch");
        if (!samples || shard_count <= 0 || shard_index < 0 || shard_index >= shard_count ||
            batches < 20 || !(proposal_p > 0.5 && proposal_p < 0.7))
            throw std::invalid_argument("invalid production arguments");
        Torus torus(6,kernel);
        const int N = torus.size();
        const std::uint64_t begin = samples * std::uint64_t(shard_index) / std::uint64_t(shard_count);
        const std::uint64_t end = samples * std::uint64_t(shard_index+1) / std::uint64_t(shard_count);
        std::vector<std::vector<McGlobalRow>> global(
            batches,std::vector<McGlobalRow>(N));
        std::array<std::vector<std::vector<McLandingRow>>,2> landing{{
            std::vector<std::vector<McLandingRow>>(batches,std::vector<McLandingRow>(N)),
            std::vector<std::vector<McLandingRow>>(batches,std::vector<McLandingRow>(N))}};
        std::vector<unsigned char> occupied(N,0);
        for (std::uint64_t sample = begin; sample < end; ++sample) {
            const int batch = int(sample % std::uint64_t(batches));
            int k = 0;
            for (int v = 1; v < N; ++v) {
                const std::uint64_t key = seed ^ (sample * 0xd6e8feb86659fd93ULL) ^
                                          (std::uint64_t(v) * 0xa0761d6478bd642fULL);
                occupied[v] = unit(key) < proposal_p;
                k += occupied[v];
            }
            occupied[0] = 0;
            const State state0 = torus.evaluate(occupied);
            const int h4 = torus.landing_h4(occupied);
            occupied[0] = 1;
            const State state1 = torus.evaluate(occupied);
            occupied[0] = 0;

            // A single uniformly keyed origin is unbiased for the complete
            // ordered source sum and keeps baseline work O(N), not O(N^2).
            const int origin = int(splitmix64(seed ^ sample ^ 0xe7037ed1a0b428dbULL) % N);
            const auto source0_est = std::int64_t(N) * torus.source16_origin(occupied,state0,origin);
            occupied[0] = 1;
            const auto source1_est = std::int64_t(N) * torus.source16_origin(occupied,state1,origin);
            occupied[0] = 0;
            auto& g = global[batch][k];
            ++g.count; g.sum_q0 += state0.q; g.sum_q1 += state1.q;
            g.sum_source16_0 += source0_est; g.sum_source16_1 += source1_est;

            const int r0 = state0.q + 1, r1 = state1.q + 1;
            int tr = -1;
            if (r0 == 0 && r1 == 1) tr = 0;
            if (r0 == 1 && r1 == 2) tr = 1;
            if (tr >= 0 && h4) {
                // The rare primary landing event receives the exact complete
                // source, avoiding an extra stochastic layer in the wedge.
                const auto source0 = torus.source16(occupied,state0);
                occupied[0] = 1;
                const auto source1 = torus.source16(occupied,state1);
                occupied[0] = 0;
                auto& a = landing[tr][batch][k];
                a.signed_count += h4;
                a.signed_source_mid16 += std::int64_t(h4) * (source0+source1);
                ++a.unsigned_count;
            }
        }
        write_mc(output,samples,shard_index,shard_count,batches,seed,proposal_p,
                 begin,end,global,landing);
        std::cerr << "completed samples=" << (end-begin) << " kernel_rows=" << kernel.rows << '\n';
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
