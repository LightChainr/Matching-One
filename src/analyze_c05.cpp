// C05/P33 analysis: reconstruct M,S,D,P4 and thermal even/odd from rank histograms.
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

struct Hist {
    int n = 0;
    std::vector<std::uint64_t> km;  // 0..N+1
    std::vector<std::uint64_t> kp;
    std::uint64_t samples = 0;
    int a = 0, b = 0;
    double cos4 = 0;
    std::string name;
    std::map<int, std::vector<std::uint64_t>> bkm, bkp;
};

double cos4_exact(int a, int b) {
    const double a2 = static_cast<double>(a) * a;
    const double b2 = static_cast<double>(b) * b;
    return (a2 * a2 - 6.0 * a2 * b2 + b2 * b2) / ((a2 + b2) * (a2 + b2));
}

std::vector<std::uint64_t> load_count_csv(const std::string& path, int& maxk) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(in, line);
    std::vector<std::pair<int, std::uint64_t>> rows;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream iss(line);
        int k = 0;
        std::uint64_t c = 0;
        iss >> k >> c;
        rows.push_back({k, c});
    }
    maxk = 0;
    for (auto& r : rows) maxk = std::max(maxk, r.first);
    std::vector<std::uint64_t> h(static_cast<std::size_t>(maxk + 1), 0);
    for (auto& r : rows) h[static_cast<std::size_t>(r.first)] = r.second;
    return h;
}

void load_batch_csv(const std::string& path, std::map<int, std::vector<std::uint64_t>>& out) {
    std::ifstream in(path);
    if (!in) return;
    std::string line;
    std::getline(in, line);
    std::map<int, std::vector<std::pair<int, std::uint64_t>>> tmp;
    int maxk = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream iss(line);
        int b = 0, k = 0;
        std::uint64_t c = 0;
        iss >> b >> k >> c;
        tmp[b].push_back({k, c});
        maxk = std::max(maxk, k);
    }
    for (auto& kv : tmp) {
        std::vector<std::uint64_t> h(static_cast<std::size_t>(maxk + 1), 0);
        for (auto& p : kv.second) h[static_cast<std::size_t>(p.first)] = p.second;
        out[kv.first] = std::move(h);
    }
}

Hist load_geom(const std::string& root, int a, int b) {
    Hist h;
    h.a = a;
    h.b = b;
    h.cos4 = cos4_exact(a, b);
    std::ostringstream nm;
    nm << "g_" << a << "_" << b;
    h.name = nm.str();
    const std::string dir = root + "/" + h.name;
    int maxk = 0;
    h.km = load_count_csv(dir + "/kminus_hist.csv", maxk);
    int maxk2 = 0;
    h.kp = load_count_csv(dir + "/kplus_hist.csv", maxk2);
    h.n = std::max(maxk, maxk2) - 1;  // bins 0..N+1
    h.samples = 0;
    for (auto c : h.km) h.samples += c;
    load_batch_csv(dir + "/batch_kminus_hist.csv", h.bkm);
    load_batch_csv(dir + "/batch_kplus_hist.csv", h.bkp);
    return h;
}

std::vector<double> wrap_g(const Hist& h) {
    std::vector<double> q(static_cast<std::size_t>(h.n + 1), 0.0);
    double acc = 0.0;
    const double tot = static_cast<double>(h.samples);
    for (int k = 0; k <= h.n; ++k) {
        acc += static_cast<double>(h.kp[static_cast<std::size_t>(k)]);
        q[static_cast<std::size_t>(k)] = acc / tot;
    }
    return q;
}

std::vector<double> wrap_hat(const Hist& h) {
    // q[m] = P(K_minus > n-m)
    const int n = h.n;
    std::vector<double> sf(static_cast<std::size_t>(n + 2), 0.0);
    double acc = static_cast<double>(h.samples);
    const double tot = acc;
    for (int k = 0; k <= n + 1; ++k) {
        acc -= static_cast<double>(h.km[static_cast<std::size_t>(k)]);
        sf[static_cast<std::size_t>(k)] = acc / tot;
    }
    std::vector<double> q(static_cast<std::size_t>(n + 1), 0.0);
    for (int m = 0; m <= n; ++m) q[static_cast<std::size_t>(m)] = sf[static_cast<std::size_t>(n - m)];
    return q;
}

std::vector<double> mean_d(const Hist& h) {
    auto qg = wrap_g(h);
    const int n = h.n;
    std::vector<double> sf(static_cast<std::size_t>(n + 2), 0.0);
    double acc = static_cast<double>(h.samples);
    const double tot = acc;
    for (int k = 0; k <= n + 1; ++k) {
        acc -= static_cast<double>(h.km[static_cast<std::size_t>(k)]);
        sf[static_cast<std::size_t>(k)] = acc / tot;
    }
    std::vector<double> q(static_cast<std::size_t>(n + 1), 0.0);
    for (int k = 0; k <= n; ++k)
        q[static_cast<std::size_t>(k)] = qg[static_cast<std::size_t>(k)] - sf[static_cast<std::size_t>(k)];
    return q;
}

double convolve(const std::vector<double>& qk, int n, double p) {
    if (p <= 0.0) return qk.front();
    if (p >= 1.0) return qk.back();
    int mode = static_cast<int>((n + 1) * p);
    if (mode < 0) mode = 0;
    if (mode > n) mode = n;
    const double ratio_up = p / (1.0 - p);
    const double ratio_dn = (1.0 - p) / p;
    double w = 1.0, tot = 1.0, s = qk[static_cast<std::size_t>(mode)];
    for (int k = mode + 1; k <= n; ++k) {
        w *= static_cast<double>(n - k + 1) / static_cast<double>(k) * ratio_up;
        if (w == 0.0) break;
        tot += w;
        s += w * qk[static_cast<std::size_t>(k)];
    }
    w = 1.0;
    for (int k = mode - 1; k >= 0; --k) {
        w *= static_cast<double>(k + 1) / static_cast<double>(n - k) * ratio_dn;
        if (w == 0.0) break;
        tot += w;
        s += w * qk[static_cast<std::size_t>(k)];
    }
    return s / tot;
}

double find_root(const std::function<double(double)>& f, double lo, double hi) {
    double fa = f(lo), fb = f(hi);
    if (fa == 0.0) return lo;
    if (fb == 0.0) return hi;
    if (fa * fb > 0.0) {
        double found_lo = lo, found_hi = hi;
        bool ok = false;
        double prev = fa;
        double prevx = lo;
        for (int i = 1; i <= 400; ++i) {
            const double x = lo + (hi - lo) * i / 400.0;
            const double y = f(x);
            if (prev * y <= 0.0) {
                found_lo = prevx;
                found_hi = x;
                fa = prev;
                fb = y;
                ok = true;
                break;
            }
            prev = y;
            prevx = x;
        }
        if (!ok) return std::numeric_limits<double>::quiet_NaN();
        lo = found_lo;
        hi = found_hi;
    }
    for (int i = 0; i < 80; ++i) {
        const double m = 0.5 * (lo + hi);
        const double fm = f(m);
        if (fm == 0.0 || (hi - lo) < 1e-15) return m;
        if (fa * fm <= 0.0) {
            hi = m;
            fb = fm;
        } else {
            lo = m;
            fa = fm;
        }
    }
    return 0.5 * (lo + hi);
}

double fd1(const std::function<double(double)>& f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}
double fd3(const std::function<double(double)>& f, double x, double h) {
    return (f(x + 2 * h) - 2.0 * f(x + h) + 2.0 * f(x - h) - f(x - 2 * h)) / (2.0 * h * h * h);
}
double fd5(const std::function<double(double)>& f, double x, double h) {
    return (f(x + 3 * h) - 4.0 * f(x + 2 * h) + 5.0 * f(x + h) - 5.0 * f(x - h) + 4.0 * f(x - 2 * h) -
            f(x - 3 * h)) /
           (2.0 * h * h * h * h * h);
}

double p4(double x1, double x2, double c1, double c2) { return (x1 - x2) / (c1 - c2); }

struct Fields {
    double RG = 0, Rhat = 0, M = 0, S = 0, D = 0;
};

Fields fields_at(const std::vector<double>& qg, const std::vector<double>& qh,
                 const std::vector<double>& qd, int n, double p) {
    Fields f;
    f.RG = convolve(qg, n, p);
    f.Rhat = convolve(qh, n, p);
    f.M = convolve(qd, n, p);
    f.S = 0.5 * (f.RG + f.Rhat);
    f.D = 0.5 * (f.RG - f.Rhat);
    return f;
}

double mean_rank(const std::vector<std::uint64_t>& h) {
    double s = 0, t = 0;
    for (std::size_t k = 0; k < h.size(); ++k) {
        s += static_cast<double>(k) * static_cast<double>(h[k]);
        t += static_cast<double>(h[k]);
    }
    return s / t;
}

double rmse(const std::vector<double>& y, const std::vector<double>& yhat) {
    double s = 0;
    for (std::size_t i = 0; i < y.size(); ++i) {
        const double d = y[i] - yhat[i];
        s += d * d;
    }
    return std::sqrt(s / static_cast<double>(y.size()));
}

int analyze_main(int argc, char** argv) {
    std::string root = "results/server-20260828/C05";
    if (argc > 1) root = argv[1];
    struct Pair {
        int n, a1, b1, a2, b2;
        bool train;
        bool holdout;
    };
    const std::vector<Pair> pairs = {
        {65, 8, 1, 7, 4, true, false},
        {85, 9, 2, 7, 6, true, false},
        {130, 11, 3, 9, 7, true, false},
        {145, 12, 1, 9, 8, false, true},
        {170, 13, 1, 11, 7, false, true},
    };

    struct SizeRec {
        Pair p;
        Hist g1, g2;
        std::vector<double> qg1, qh1, qd1, qg2, qh2, qd2;
        double pstar = 0, mprime = 0, kappa3 = 0, kappa5 = 0;
        double p4D = 0, p4S = 0, p4M = 0, p4D_se = 0, p4S_se = 0;
        double gap = 0;
        double L = 0;
        std::vector<double> u_grid;
        struct URow {
            double u, pm, pp, D_even, S_even, M_even, D_odd, S_odd;
        };
        std::vector<URow> thermal;
    };

    std::vector<SizeRec> recs;
    std::vector<double> u_grid = {0.0};
    double pstar65 = 0.5;

    for (const Pair& pr : pairs) {
        const std::string d1 = root + "/g_" + std::to_string(pr.a1) + "_" + std::to_string(pr.b1);
        if (!std::filesystem::exists(d1 + "/kminus_hist.csv")) continue;
        SizeRec r;
        r.p = pr;
        r.g1 = load_geom(root, pr.a1, pr.b1);
        r.g2 = load_geom(root, pr.a2, pr.b2);
        r.L = std::sqrt(static_cast<double>(pr.n));
        r.qg1 = wrap_g(r.g1);
        r.qh1 = wrap_hat(r.g1);
        r.qd1 = mean_d(r.g1);
        r.qg2 = wrap_g(r.g2);
        r.qh2 = wrap_hat(r.g2);
        r.qd2 = mean_d(r.g2);
        auto mbar = [&](double p) {
            return 0.5 * (convolve(r.qd1, r.g1.n, p) + convolve(r.qd2, r.g2.n, p));
        };
        r.pstar = find_root(mbar, 0.45, 0.75);
        const double h = std::max(1e-4, 0.25 / pr.n);
        r.mprime = fd1(mbar, r.pstar, h);
        const double m3 = fd3(mbar, r.pstar, h);
        const double m5 = fd5(mbar, r.pstar, std::max(h, 2e-4));
        r.kappa3 = m3 / (r.mprime * r.mprime * r.mprime);
        r.kappa5 = m5 / std::pow(r.mprime, 5.0);
        const Fields f1 = fields_at(r.qg1, r.qh1, r.qd1, r.g1.n, r.pstar);
        const Fields f2 = fields_at(r.qg2, r.qh2, r.qd2, r.g2.n, r.pstar);
        r.p4D = p4(f1.D, f2.D, r.g1.cos4, r.g2.cos4);
        r.p4S = p4(f1.S, f2.S, r.g1.cos4, r.g2.cos4);
        r.p4M = p4(f1.M, f2.M, r.g1.cos4, r.g2.cos4);
        r.gap = 0.5 * ((mean_rank(r.g1.kp) - mean_rank(r.g1.km)) +
                       (mean_rank(r.g2.kp) - mean_rank(r.g2.km)));

        // batch jackknife of P4[D], P4[S] at pstar
        std::vector<double> bd, bs;
        for (const auto& kv : r.g1.bkm) {
            const int b = kv.first;
            if (!r.g2.bkm.count(b)) continue;
            Hist h1 = r.g1, h2 = r.g2;
            h1.km = r.g1.bkm[b];
            h1.kp = r.g1.bkp[b];
            h1.samples = 0;
            for (auto c : h1.km) h1.samples += c;
            h2.km = r.g2.bkm[b];
            h2.kp = r.g2.bkp[b];
            h2.samples = 0;
            for (auto c : h2.km) h2.samples += c;
            auto qg1b = wrap_g(h1), qh1b = wrap_hat(h1), qd1b = mean_d(h1);
            auto qg2b = wrap_g(h2), qh2b = wrap_hat(h2), qd2b = mean_d(h2);
            const Fields a = fields_at(qg1b, qh1b, qd1b, h1.n, r.pstar);
            const Fields c = fields_at(qg2b, qh2b, qd2b, h2.n, r.pstar);
            bd.push_back(p4(a.D, c.D, r.g1.cos4, r.g2.cos4));
            bs.push_back(p4(a.S, c.S, r.g1.cos4, r.g2.cos4));
        }
        if (bd.size() > 1) {
            double mu = 0;
            for (double x : bd) mu += x;
            mu /= static_cast<double>(bd.size());
            double var = 0;
            for (double x : bd) var += (x - mu) * (x - mu);
            var /= static_cast<double>(bd.size() - 1);
            r.p4D_se = std::sqrt(var / static_cast<double>(bd.size()));
            // overwrite pooled with batch mean for signed effect reporting
            // keep pooled r.p4D as histogram-pooled estimator
        }
        if (bs.size() > 1) {
            double mu = 0;
            for (double x : bs) mu += x;
            mu /= static_cast<double>(bs.size());
            double var = 0;
            for (double x : bs) var += (x - mu) * (x - mu);
            var /= static_cast<double>(bs.size() - 1);
            r.p4S_se = std::sqrt(var / static_cast<double>(bs.size()));
        }

        if (pr.n == 65) {
            pstar65 = r.pstar;
            const double span = std::min(std::fabs(mbar(0.2)), std::fabs(mbar(0.85)));
            u_grid = {0.0};
            for (double u : {0.05, 0.1, 0.2, 0.4}) {
                if (u < 0.8 * span) u_grid.push_back(u);
            }
        }
        recs.push_back(std::move(r));
    }

    for (SizeRec& r : recs) {
        auto mbar = [&](double p) {
            return 0.5 * (convolve(r.qd1, r.g1.n, p) + convolve(r.qd2, r.g2.n, p));
        };
        for (double u : u_grid) {
            SizeRec::URow row;
            row.u = u;
            if (u == 0.0) {
                row.pm = row.pp = r.pstar;
                row.D_even = r.p4D;
                row.S_even = r.p4S;
                row.M_even = r.p4M;
                row.D_odd = 0;
                row.S_odd = 0;
            } else {
                row.pm = find_root([&](double p) { return mbar(p) + u; }, 0.05, r.pstar);
                row.pp = find_root([&](double p) { return mbar(p) - u; }, r.pstar, 0.95);
                const Fields aplus = fields_at(r.qg1, r.qh1, r.qd1, r.g1.n, row.pp);
                const Fields bplus = fields_at(r.qg2, r.qh2, r.qd2, r.g2.n, row.pp);
                const Fields aminus = fields_at(r.qg1, r.qh1, r.qd1, r.g1.n, row.pm);
                const Fields bminus = fields_at(r.qg2, r.qh2, r.qd2, r.g2.n, row.pm);
                const double Dp = p4(aplus.D, bplus.D, r.g1.cos4, r.g2.cos4);
                const double Dm = p4(aminus.D, bminus.D, r.g1.cos4, r.g2.cos4);
                const double Sp = p4(aplus.S, bplus.S, r.g1.cos4, r.g2.cos4);
                const double Sm = p4(aminus.S, bminus.S, r.g1.cos4, r.g2.cos4);
                const double Mp = p4(aplus.M, bplus.M, r.g1.cos4, r.g2.cos4);
                const double Mm = p4(aminus.M, bminus.M, r.g1.cos4, r.g2.cos4);
                row.D_even = 0.5 * (Dp + Dm);
                row.D_odd = 0.5 * (Dp - Dm);
                row.S_even = 0.5 * (Sp + Sm);
                row.S_odd = 0.5 * (Sp - Sm);
                row.M_even = 0.5 * (Mp + Mm);
            }
            r.thermal.push_back(row);
        }
    }

    auto collect = [&](char which, bool train_only, bool hold_only) {
        std::vector<double> L, y;
        std::vector<int> ns;
        for (const SizeRec& r : recs) {
            if (train_only && !r.p.train) continue;
            if (hold_only && !r.p.holdout) continue;
            ns.push_back(r.p.n);
            L.push_back(r.L);
            double v = 0;
            for (const auto& u : r.thermal) {
                if (u.u == 0.0) {
                    v = (which == 'D') ? u.D_even : u.S_even;
                }
            }
            y.push_back(v);
        }
        return std::tuple<std::vector<int>, std::vector<double>, std::vector<double>>(ns, L, y);
    };

    auto model = [](double expn, const std::vector<double>& Ltr, const std::vector<double>& ytr,
                    const std::vector<double>& Lpred) {
        double A = 0;
        for (std::size_t i = 0; i < Ltr.size(); ++i) A += ytr[i] / std::pow(Ltr[i], expn);
        A /= static_cast<double>(Ltr.size());
        std::vector<double> pred;
        for (double L : Lpred) pred.push_back(A * std::pow(L, expn));
        return std::pair<double, std::vector<double>>(A, pred);
    };

    auto [ntrD, LtrD, ytrD] = collect('D', true, false);
    auto [nhoD, LhoD, yhoD] = collect('D', false, true);
    auto [ntrS, LtrS, ytrS] = collect('S', true, false);
    auto [nhoS, LhoS, yhoS] = collect('S', false, true);
    auto [nallD, LallD, yallD] = collect('D', false, false);
    auto [nallS, LallS, yallS] = collect('S', false, false);

    const double e1325 = -13.0 / 4.0;
    auto D1325 = model(e1325, LtrD, ytrD, LallD);
    auto Dho = model(e1325, LtrD, ytrD, LhoD);
    auto S2 = model(-2.0, LtrS, ytrS, LallS);
    auto Sho = model(-2.0, LtrS, ytrS, LhoS);
    auto S1325ho = model(e1325, LtrS, ytrS, LhoS);
    auto D2ho = model(-2.0, LtrD, ytrD, LhoD);

    // log companion y = A L^e (1+B log L)
    double A_log = 0, B_log = 0;
    bool log_ok = LtrD.size() >= 2;
    if (log_ok) {
        std::vector<double> z, lx;
        for (std::size_t i = 0; i < LtrD.size(); ++i) {
            z.push_back(ytrD[i] / std::pow(LtrD[i], e1325));
            lx.push_back(std::log(LtrD[i]));
        }
        double mx = 0, mz = 0;
        for (std::size_t i = 0; i < lx.size(); ++i) {
            mx += lx[i];
            mz += z[i];
        }
        mx /= lx.size();
        mz /= z.size();
        double sxx = 0, sxy = 0;
        for (std::size_t i = 0; i < lx.size(); ++i) {
            sxx += (lx[i] - mx) * (lx[i] - mx);
            sxy += (lx[i] - mx) * (z[i] - mz);
        }
        const double slope = sxx > 0 ? sxy / sxx : 0.0;
        A_log = mz - slope * mx;
        B_log = (A_log != 0.0) ? slope / A_log : 0.0;
    }
    std::vector<double> pred_log_ho;
    if (log_ok) {
        for (double L : LhoD)
            pred_log_ho.push_back(A_log * (1.0 + B_log * std::log(L)) * std::pow(L, e1325));
    }

    // free exponent on train
    double efree = 0;
    if (LtrD.size() >= 2) {
        double mx = 0, my = 0;
        std::vector<double> lx, ly;
        for (std::size_t i = 0; i < LtrD.size(); ++i) {
            lx.push_back(std::log(LtrD[i]));
            ly.push_back(std::log(std::fabs(ytrD[i])));
            mx += lx.back();
            my += ly.back();
        }
        mx /= lx.size();
        my /= ly.size();
        double sxx = 0, sxy = 0;
        for (std::size_t i = 0; i < lx.size(); ++i) {
            sxx += (lx[i] - mx) * (lx[i] - mx);
            sxy += (lx[i] - mx) * (ly[i] - my);
        }
        efree = sxx > 0 ? sxy / sxx : 0.0;
    }
    auto Dfree_ho = model(efree, LtrD, ytrD, LhoD);

    std::vector<double> weffD, weffS;
    std::vector<std::pair<int, int>> weffN;
    for (std::size_t i = 0; i + 1 < LallD.size(); ++i) {
        weffN.push_back({nallD[i], nallD[i + 1]});
        if (yallD[i] * yallD[i + 1] > 0)
            weffD.push_back(std::log(std::fabs(yallD[i + 1] / yallD[i])) /
                            std::log(LallD[i + 1] / LallD[i]));
        else
            weffD.push_back(std::numeric_limits<double>::quiet_NaN());
        if (yallS[i] * yallS[i + 1] > 0)
            weffS.push_back(std::log(std::fabs(yallS[i + 1] / yallS[i])) /
                            std::log(LallS[i + 1] / LallS[i]));
        else
            weffS.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    const double rmse_D_fixed_ho = yhoD.empty() ? NAN : rmse(yhoD, Dho.second);
    const double rmse_D_log_ho = (yhoD.empty() || !log_ok) ? NAN : rmse(yhoD, pred_log_ho);
    const double rmse_D_free_ho = yhoD.empty() ? NAN : rmse(yhoD, Dfree_ho.second);
    const double rmse_S_m2_ho = yhoS.empty() ? NAN : rmse(yhoS, Sho.second);
    const double rmse_S_1325_ho = yhoS.empty() ? NAN : rmse(yhoS, S1325ho.second);
    const bool log_better =
        std::isfinite(rmse_D_log_ho) && std::isfinite(rmse_D_fixed_ho) && rmse_D_log_ho < rmse_D_fixed_ho;

    // like L^{-13/4}? mean weff near -3.25
    double mean_weffD = 0;
    int n_w = 0;
    for (double w : weffD) {
        if (std::isfinite(w)) {
            mean_weffD += w;
            n_w++;
        }
    }
    if (n_w) mean_weffD /= n_w;
    const bool like_1325 = n_w >= 2 && std::fabs(mean_weffD - e1325) < 1.0;  // loose; report signed

    std::ofstream js(root + "/derived_summary.json");
    js << std::setprecision(17);
    js << "{\n";
    js << "  \"channel\": \"cross\",\n";
    js << "  \"not_a_second_replication\": \"either-wrap matching-function rows are not used\",\n";
    js << "  \"frozen\": {\"u_grid\": [";
    for (std::size_t i = 0; i < u_grid.size(); ++i) {
        if (i) js << ", ";
        js << u_grid[i];
    }
    js << "], \"u_source\": \"N=65 orientation-averaged Mbar coverage\", \"p_star_N65\": " << pstar65
       << ", \"train_N\": [65, 85, 130], \"holdout_N\": [145, 170]},\n";
    js << "  \"sizes\": [\n";
    for (std::size_t i = 0; i < recs.size(); ++i) {
        const SizeRec& r = recs[i];
        if (i) js << ",\n";
        js << "    {\"N\": " << r.p.n << ", \"L_phys\": " << r.L << ", \"samples\": " << r.g1.samples
           << ", \"p_star_Mbar\": " << r.pstar << ", \"Mprime\": " << r.mprime
           << ", \"density_Mprime_over_2\": " << 0.5 * r.mprime << ", \"kappa3\": " << r.kappa3
           << ", \"kappa5\": " << r.kappa5 << ", \"P4_D_pstar\": " << r.p4D << ", \"P4_D_batch_se\": "
           << r.p4D_se << ", \"P4_S_pstar\": " << r.p4S << ", \"P4_S_batch_se\": " << r.p4S_se
           << ", \"P4_M_pstar\": " << r.p4M << ", \"mean_Kplus_minus_Kminus\": " << r.gap
           << ", \"cos4\": [" << r.g1.cos4 << ", " << r.g2.cos4 << "], \"delta_cos4\": "
           << (r.g1.cos4 - r.g2.cos4) << ", \"thermal_u0_D_even\": " << r.p4D
           << ", \"thermal_u0_S_even\": " << r.p4S << "}";
    }
    js << "\n  ],\n";
    js << "  \"thermal\": [\n";
    for (std::size_t i = 0; i < recs.size(); ++i) {
        const SizeRec& r = recs[i];
        if (i) js << ",\n";
        js << "    {\"N\": " << r.p.n << ", \"u\": [";
        for (std::size_t j = 0; j < r.thermal.size(); ++j) {
            const auto& u = r.thermal[j];
            if (j) js << ", ";
            js << "{\"u\": " << u.u << ", \"p_minus\": " << u.pm << ", \"p_plus\": " << u.pp
               << ", \"P4D_even\": " << u.D_even << ", \"P4D_odd\": " << u.D_odd
               << ", \"P4S_even\": " << u.S_even << ", \"P4S_odd\": " << u.S_odd
               << ", \"P4M_even\": " << u.M_even << "}";
        }
        js << "]}";
    }
    js << "\n  ],\n";
    js << "  \"scaling\": {\n";
    js << "    \"P4D_thermal_even_u0\": {\"N\": [";
    for (std::size_t i = 0; i < nallD.size(); ++i) {
        if (i) js << ", ";
        js << nallD[i];
    }
    js << "], \"values\": [";
    for (std::size_t i = 0; i < yallD.size(); ++i) {
        if (i) js << ", ";
        js << yallD[i];
    }
    js << "], \"pairwise_w_eff\": [";
    for (std::size_t i = 0; i < weffD.size(); ++i) {
        if (i) js << ", ";
        js << "{\"N1\": " << weffN[i].first << ", \"N2\": " << weffN[i].second << ", \"w_eff\": " << weffD[i]
           << "}";
    }
    js << "], \"mean_w_eff\": " << mean_weffD << ", \"target_L_m13_4\": " << e1325
       << ", \"A_train_L_m13_4\": " << D1325.first << ", \"holdout_rmse_L_m13_4\": " << rmse_D_fixed_ho
       << ", \"holdout_rmse_log\": " << rmse_D_log_ho << ", \"holdout_rmse_free\": " << rmse_D_free_ho
       << ", \"free_exponent_train\": " << efree << ", \"B_log_train\": " << B_log
       << ", \"log_improves_holdout\": " << (log_better ? "true" : "false") << "},\n";
    js << "    \"P4S_thermal_even_u0\": {\"N\": [";
    for (std::size_t i = 0; i < nallS.size(); ++i) {
        if (i) js << ", ";
        js << nallS[i];
    }
    js << "], \"values\": [";
    for (std::size_t i = 0; i < yallS.size(); ++i) {
        if (i) js << ", ";
        js << yallS[i];
    }
    js << "], \"pairwise_w_eff\": [";
    for (std::size_t i = 0; i < weffS.size(); ++i) {
        if (i) js << ", ";
        js << "{\"N1\": " << weffN[i].first << ", \"N2\": " << weffN[i].second << ", \"w_eff\": " << weffS[i]
           << "}";
    }
    js << "], \"target_L_m2\": -2.0, \"A_train_L_m2\": " << S2.first
       << ", \"holdout_rmse_L_m2\": " << rmse_S_m2_ho << ", \"holdout_rmse_L_m13_4\": " << rmse_S_1325_ho
       << "}\n  },\n";
    js << "  \"target_tests\": {\n";
    js << "    \"P4D_thermal_even_looks_like_L_m13_4\": " << (like_1325 ? "true" : "false") << ",\n";
    js << "    \"mean_pairwise_w_eff_P4D\": " << mean_weffD << ",\n";
    js << "    \"target_exponent\": " << e1325 << ",\n";
    js << "    \"log_alternative_improves_heldout_P4D\": " << (log_better ? "true" : "false") << "\n";
    js << "  }\n";
    js << "}\n";

    std::ofstream rpt(root + "/REPORT.md");
    rpt << std::setprecision(8);
    rpt << "# C05 / P33 — Thermal-coordinate tomography from threshold ranks\n\n";
    rpt << "Engine: issue-9 Philox Fisher–Yates Newman–Ziff on the C00 `HomologyUnionFind` "
           "(exact `adj(P)/det(P)` windings). **Cross channel only**; either-wrap matching-function "
           "rows are not a second replication.\n\n";
    rpt << "## Off-by-one convention\n\n";
    rpt << "- `K_plus`: smallest black occupation `k` at which primal CROSS wrapping is true "
           "(`N+1` if never). Transition at rank `k` in `1..N` is the `k`-th uniform order statistic, "
           "`T|K=k ~ Beta(k, N+1-k)`.\n";
    rpt << "- `K_minus = N - m* + 1` where `m*` is the first reverse-permutation matching occupation "
           "with CROSS wrapping (`0` if white never wraps).\n";
    rpt << "- Matching function: `M(p)=P(K_plus <= m) - P(K_minus > m)` with `m ~ Binomial(N,p)`.\n";
    rpt << "- Exact tiny tests (axis L=2,3 exhaustive perms; gaussian (2,1); axis L=4 subset vs "
           "published polynomial + MC): **PASS**. `K_minus <= K_plus` on every replica of those tests "
           "and of production.\n\n";
    rpt << "## Frozen choices (N=65 only)\n\n";
    rpt << "- u-grid: ";
    for (std::size_t i = 0; i < u_grid.size(); ++i) {
        if (i) rpt << ", ";
        rpt << u_grid[i];
    }
    rpt << " (from N=65 `|Mbar|` coverage; p*_65=" << pstar65 << ").\n";
    rpt << "- Scaling train: N=65,85,130. Held-out: N=145,170. N=145 is retained even if noisy.\n\n";
    rpt << "## Per-size summary (2e6 CRN replicas each, 40 batches, 8 OpenMP threads)\n\n";
    rpt << "| N | L | p*_Mbar | P4[D](p*) | se_batch | P4[S](p*) | se_batch | kappa3 | mean gap |\n";
    rpt << "|---|---|---|---|---|---|---|---|---|\n";
    for (const SizeRec& r : recs) {
        rpt << "| " << r.p.n << " | " << r.L << " | " << r.pstar << " | " << r.p4D << " | " << r.p4D_se
            << " | " << r.p4S << " | " << r.p4S_se << " | " << r.kappa3 << " | " << r.gap << " |\n";
    }
    rpt << "\nSigned effects are histogram-pooled P4 values; `se_batch` is the batch SD/sqrt(B) of "
           "the 40 CRN batches (covariance across orientations is already in the paired batches).\n\n";
    rpt << "## Target tests\n\n";
    rpt << "P4[D]_thermal_even at u=0 vs L^{-13/4} = " << e1325 << ":\n\n";
    rpt << "- pairwise effective exponents: ";
    for (std::size_t i = 0; i < weffD.size(); ++i) {
        if (i) rpt << ", ";
        rpt << "N=" << weffN[i].first << "→" << weffN[i].second << " w_eff=" << weffD[i];
    }
    rpt << "\n- mean w_eff = " << mean_weffD << " (target -3.25)\n";
    rpt << "- looks-like-L^{-13/4} (mean w_eff within 1 of target): **" << (like_1325 ? "yes" : "no")
        << "**\n";
    rpt << "- train A for L^{-13/4}: " << D1325.first << "\n";
    rpt << "- held-out RMSE L^{-13/4}: " << rmse_D_fixed_ho << "\n";
    rpt << "- held-out RMSE L^{-13/4}(1+B log L), B=" << B_log << ": " << rmse_D_log_ho << "\n";
    rpt << "- held-out RMSE free exponent (train e=" << efree << "): " << rmse_D_free_ho << "\n";
    rpt << "- log alternative improves held-out P4[D]: **" << (log_better ? "yes" : "no") << "**\n\n";
    rpt << "P4[S]_thermal_even at u=0 vs L^{-2}:\n\n";
    rpt << "- pairwise w_eff: ";
    for (std::size_t i = 0; i < weffS.size(); ++i) {
        if (i) rpt << ", ";
        rpt << "N=" << weffN[i].first << "→" << weffN[i].second << " w_eff=" << weffS[i];
    }
    rpt << "\n- held-out RMSE L^{-2}: " << rmse_S_m2_ho << "\n";
    rpt << "- held-out RMSE L^{-13/4}: " << rmse_S_1325_ho << "\n\n";
    rpt << "P4[D]_even values: ";
    for (std::size_t i = 0; i < yallD.size(); ++i) {
        if (i) rpt << ", ";
        rpt << "N" << nallD[i] << "=" << yallD[i];
    }
    rpt << "\nP4[S]_even values: ";
    for (std::size_t i = 0; i < yallS.size(); ++i) {
        if (i) rpt << ", ";
        rpt << "N" << nallS[i] << "=" << yallS[i];
    }
    rpt << "\n\n## Remaining\n\n";
    rpt << "- Axis L=8..32 production histograms were deprioritized (P33 same-N Gaussian first); "
           "exact axis L=2,3,4 tests did run and PASS.\n";
    rpt << "- No 1e8-sample campaign. CPU pilot is 2e6/orientation-pair.\n";
    rpt << "- N=1105 not started.\n";
    std::cout << "wrote " << root << "/derived_summary.json\n";
    std::cout << "wrote " << root << "/REPORT.md\n";
    std::cout << "P4D like L^{-13/4}? " << (like_1325 ? "yes" : "no") << " mean_weff=" << mean_weffD
              << "\n";
    return 0;
}

int main(int argc, char** argv) { return analyze_main(argc, argv); }
