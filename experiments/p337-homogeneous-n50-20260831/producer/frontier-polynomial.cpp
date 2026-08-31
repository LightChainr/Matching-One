// Exact black-NN lifted frontier, with K carried by dynamic value polynomials.
// No scientific root/source score.  Target resource probes must stop before N.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sys/resource.h>
#include <time.h>

using I128 = __int128_t;
constexpr int MAX_FRONT = 128;
volatile std::sig_atomic_t cpu_signal = 0;

void on_cpu(int) { cpu_signal = 1; }

double cpu_seconds() {
    timespec t{};
    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t)) throw std::runtime_error("CPU clock failed");
    return double(t.tv_sec) + double(t.tv_nsec)*1e-9;
}

double peak_rss_mib() {
    rusage r{};
    if (getrusage(RUSAGE_SELF, &r)) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return double(r.ru_maxrss)/(1024.0*1024.0);
#else
    return double(r.ru_maxrss)/1024.0;
#endif
}

std::string decimal(I128 v) {
    if (!v) return "0";
    bool negative = v < 0;
    if (negative) v = -v;
    std::string s;
    while (v) { s.push_back(char('0' + v%10)); v /= 10; }
    if (negative) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

struct Edge { int position, dx, dy; };
struct Layer { int before, after; std::vector<int> retained; std::vector<Edge> edges; };
struct Graph { int a, b, n; std::vector<Layer> layers; };

Graph read_graph(const std::string& path) {
    std::ifstream in(path);
    std::string magic;
    Graph g{};
    if (!(in >> magic >> g.a >> g.b >> g.n) || magic != "P337_BLACK_FRONTIER_V1")
        throw std::runtime_error("invalid geometry header");
    // |partial S| <= 6N+1, so all count/source coefficients fit signed64 at N<=50.
    // Arithmetic below still uses checked signed128 intermediates.
    if (g.n < 1 || g.n > 50) throw std::runtime_error("geometry outside N<=50 scope");
    int expected_before = 0;
    for (int v=0; v<g.n; ++v) {
        Layer x{}; int edge_count;
        if (!(in >> x.before >> x.after >> edge_count)) throw std::runtime_error("truncated layer");
        if (x.before != expected_before || x.before+1 > MAX_FRONT || x.after < 0 || x.after > x.before+1)
            throw std::runtime_error("inconsistent frontier lengths");
        std::array<bool, MAX_FRONT> seen{};
        for (int j=0; j<x.after; ++j) {
            int p; in >> p;
            if (!in || p<0 || p>x.before || seen[p]) throw std::runtime_error("invalid retained position");
            seen[p]=true; x.retained.push_back(p);
        }
        if (edge_count<0 || edge_count>4) throw std::runtime_error("not the fixed square NN graph");
        for (int j=0; j<edge_count; ++j) {
            Edge e{}; in >> e.position >> e.dx >> e.dy;
            if (!in || e.position<0 || e.position>=x.before || std::abs(e.dx)+std::abs(e.dy)!=1)
                throw std::runtime_error("invalid physical NN edge displacement");
            x.edges.push_back(e);
        }
        expected_before=x.after;
        g.layers.push_back(std::move(x));
    }
    if (expected_before != 0) throw std::runtime_error("terminal frontier must be empty");
    std::string extra;
    if (in >> extra) throw std::runtime_error("unexpected geometry tokens");
    return g;
}

int checked_int32(std::int64_t x) {
    if (x < std::numeric_limits<std::int32_t>::min() || x > std::numeric_limits<std::int32_t>::max())
        throw std::runtime_error("gain exceeds int32 working range");
    return static_cast<int>(x);
}

void put16(std::string& key, int value) {
    if (value < -32768 || value > 32767) throw std::runtime_error("gain exceeds signed int16 packed range");
    auto u=static_cast<std::uint16_t>(static_cast<std::int16_t>(value));
    key.push_back(static_cast<char>(u&255));
    key.push_back(static_cast<char>((u>>8)&255));
}

int get16(const std::string& key, std::size_t& pos) {
    if (pos+2>key.size()) throw std::runtime_error("truncated packed gain");
    auto lo=static_cast<unsigned char>(key[pos++]);
    auto hi=static_cast<unsigned char>(key[pos++]);
    int value=int(lo) + (int(hi)<<8);
    return value>=32768 ? value-65536 : value;
}

struct State {
    int r=0, hx=0, hy=0;
    std::array<int, MAX_FRONT> label{}, px{}, py{};
};

State decode(const std::string& key, int width) {
    if (key.empty()) throw std::runtime_error("truncated state header");
    State s;
    std::size_t pos=0;
    s.r=static_cast<unsigned char>(key[pos++]);
    if (s.r>2) throw std::runtime_error("invalid ambient rank");
    if (s.r==1) { s.hx=get16(key,pos); s.hy=get16(key,pos); }
    for (int i=0; i<width; ++i) {
        if (pos>=key.size()) throw std::runtime_error("truncated packed label");
        s.label[i]=static_cast<unsigned char>(key[pos++]);
        if (s.label[i]>=MAX_FRONT) throw std::runtime_error("packed label outside range");
        if (s.r<2) s.px[i]=get16(key,pos);
        if (s.r==0) s.py[i]=get16(key,pos);
        if (s.label[i]==0 && (s.px[i] || s.py[i])) throw std::runtime_error("vacant site carries a gain");
    }
    if (pos!=key.size()) throw std::runtime_error("packed state length mismatch");
    return s;
}

std::string encode(const State& s, const std::vector<int>& retained) {
    if (s.r<0 || s.r>2) throw std::runtime_error("invalid state header");
    if (s.r==1 && (std::gcd(std::abs(s.hx),std::abs(s.hy))!=1 || s.hx<0 || (s.hx==0 && s.hy<=0)))
        throw std::runtime_error("rank-one direction must be primitive with canonical sign");
    if (s.r!=1 && (s.hx || s.hy)) throw std::runtime_error("spurious global direction");
    std::string key;
    key.reserve(1+(s.r==1 ? 4:0)+retained.size()*(s.r==0 ? 5 : s.r==1 ? 3:1));
    key.push_back(static_cast<char>(s.r));
    if (s.r==1) { put16(key,s.hx); put16(key,s.hy); }
    std::array<int, MAX_FRONT> map{}, ax{}, ay{};
    int next=0;
    for (int i:retained) {
        int label=s.label[i], x=0, y=0, canon=0;
        if (label<0 || label>=MAX_FRONT) throw std::runtime_error("working label outside range");
        if (label) {
            if (!map[label]) { map[label]=++next; ax[label]=s.px[i]; ay[label]=s.py[i]; }
            canon=map[label];
            x=checked_int32(std::int64_t(s.px[i])-ax[label]);
            y=checked_int32(std::int64_t(s.py[i])-ay[label]);
        }
        key.push_back(static_cast<char>(canon));
        if (s.r<2) put16(key,x);
        if (s.r==0) put16(key,y);
        if (s.r==2 && (x || y)) throw std::runtime_error("rank2 should have no gains");
    }
    return key;
}

int advance(State& s, const Layer& layer, bool occupied) {
    const int v=layer.before;
    s.label[v]=0; s.px[v]=s.py[v]=0;
    if (!occupied) return 0;
    int largest=0;
    for (int i=0;i<v;++i) largest=std::max(largest,s.label[i]);
    s.label[v]=largest+1;
    int ds=-3;
    for (const Edge& edge:layer.edges) {
        int u=edge.position;
        if (!s.label[u]) continue;
        int ex=edge.dx, ey=edge.dy;
        if (s.r==1) { ex=checked_int32(std::int64_t(s.hx)*ey-std::int64_t(s.hy)*ex); ey=0; }
        if (s.r==2) ex=ey=0;
        int dx=checked_int32(std::int64_t(s.px[u])+ex-s.px[v]);
        int dy=checked_int32(std::int64_t(s.py[u])+ey-s.py[v]);
        int cu=s.label[u], cv=s.label[v];
        if (cu!=cv) {
            for (int j=0;j<=v;++j) if (s.label[j]==cv) {
                s.label[j]=cu;
                s.px[j]=checked_int32(std::int64_t(s.px[j])+dx);
                s.py[j]=checked_int32(std::int64_t(s.py[j])+dy);
            }
        } else {
            ds+=2; // every redundant edge increments ordinary graph cycle rank
            if (s.r==0 && (dx || dy)) {
                int divisor=std::gcd(std::abs(dx),std::abs(dy));
                s.hx=dx/divisor; s.hy=dy/divisor;
                if (s.hx<0 || (s.hx==0 && s.hy<0)) { s.hx=-s.hx; s.hy=-s.hy; }
                for (int j=0;j<=v;++j) {
                    s.px[j]=checked_int32(std::int64_t(s.hx)*s.py[j]-std::int64_t(s.hy)*s.px[j]);
                    s.py[j]=0;
                }
                s.r=1; --ds;
            } else if (s.r==1 && dx) {
                s.r=2; s.hx=s.hy=0;
                for (int j=0;j<=v;++j) s.px[j]=s.py[j]=0;
                --ds;
            }
        }
    }
    return ds;
}

struct Coefficient { std::int64_t count=0, sum_s=0; };
struct Polynomial { int low=0; std::vector<Coefficient> coefficients; };
struct HistogramValue { I128 count=0, sum_s=0; };
using Table=std::unordered_map<std::string,Polynomial>;

std::int64_t checked64(I128 v) {
    if (v<std::numeric_limits<std::int64_t>::min() || v>std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("signed64 polynomial coefficient overflow");
    return static_cast<std::int64_t>(v);
}

void add_shifted(Polynomial& target, const Polynomial& source, int occupied, int ds) {
    if (source.coefficients.empty()) throw std::runtime_error("empty source polynomial");
    int incoming_low=source.low+occupied;
    int incoming_high=incoming_low+static_cast<int>(source.coefficients.size())-1;
    if (incoming_low<0 || incoming_high>50) throw std::runtime_error("K support outside frozen N<=50 scope");
    if (target.coefficients.empty()) {
        target.low=incoming_low;
        target.coefficients.resize(source.coefficients.size());
    } else {
        int low=std::min(target.low,incoming_low);
        int high=std::max(target.low+static_cast<int>(target.coefficients.size())-1,incoming_high);
        if (low!=target.low || high!=target.low+static_cast<int>(target.coefficients.size())-1) {
            std::vector<Coefficient> expanded(static_cast<std::size_t>(high-low+1));
            std::copy(target.coefficients.begin(),target.coefficients.end(),expanded.begin()+(target.low-low));
            target.coefficients.swap(expanded);
            target.low=low;
        }
    }
    std::size_t offset=static_cast<std::size_t>(incoming_low-target.low);
    for (std::size_t i=0;i<source.coefficients.size();++i) {
        const Coefficient& from=source.coefficients[i];
        if (!from.count) {
            if (from.sum_s) throw std::runtime_error("zero-mass cell has a nonzero source sum");
            continue;
        }
        if (from.count<0) throw std::runtime_error("negative count coefficient");
        Coefficient& to=target.coefficients[offset+i];
        to.count=checked64(I128(to.count)+from.count);
        to.sum_s=checked64(I128(to.sum_s)+from.sum_s+I128(ds)*from.count);
    }
}

struct TableStats {
    std::size_t nonzero_cells=0, interval_cells=0, capacity_cells=0;
    I128 mass=0;
};

TableStats summarize(const Table& table) {
    TableStats stats;
    for (const auto& item:table) {
        stats.interval_cells+=item.second.coefficients.size();
        stats.capacity_cells+=item.second.coefficients.capacity();
        for (const Coefficient& cell:item.second.coefficients) {
            stats.nonzero_cells+=bool(cell.count);
            stats.mass+=cell.count;
        }
    }
    return stats;
}
struct Limit : std::runtime_error { using std::runtime_error::runtime_error; };
struct Config {
    std::string graph, prefix, authorization_commit;
    double seconds=10, rss=512;
    std::size_t cap=1000000;
    int max_layers=-1;
    bool resource_only=false;
};

Config arguments(int argc,char**argv) {
    Config c;
    for (int i=1;i<argc;++i) {
        std::string a=argv[i];
        if (a=="--resource-only") {c.resource_only=true;continue;}
        if (i+1>=argc) throw std::runtime_error("missing argument value");
        std::string value=argv[++i];
        if (a=="--graph") c.graph=value;
        else if (a=="--output-prefix") c.prefix=value;
        else if (a=="--authorization-commit") c.authorization_commit=value;
        else if (a=="--cpu-seconds") c.seconds=std::stod(value);
        else if (a=="--rss-mib") c.rss=std::stod(value);
        else if (a=="--state-cap") c.cap=std::stoull(value);
        else if (a=="--max-layers") c.max_layers=std::stoi(value);
        else throw std::runtime_error("unknown option "+a);
    }
    if (c.graph.empty() || c.prefix.empty() || c.seconds<=0 || c.seconds>600 || c.rss<=0 || c.rss>18432 || !c.cap || c.cap>5000000)
        throw std::runtime_error("explicit graph/output and CPU<=600s, RSS<=18432MiB, states<=5e6 required");
    return c;
}

int run(const Config& config) {
    Graph graph=read_graph(config.graph);
    int max_layers=config.max_layers<0 ? graph.n : config.max_layers;
    if (max_layers<1 || max_layers>graph.n) throw std::runtime_error("invalid layer cap");
    bool authorized=config.authorization_commit.size()==40 && std::all_of(
        config.authorization_commit.begin(),config.authorization_commit.end(),
        [](char ch){return (ch>='0' && ch<='9') || (ch>='a' && ch<='f');});
    if (!config.authorization_commit.empty() && !authorized)
        throw std::runtime_error("authorization commit must be a full lowercase SHA");
    // The invoking freeze-verification driver must verify this commit and source
    // hashes.  The C++ engine records it, but does not claim to inspect Git itself.
    if (graph.n==50 && !authorized && (max_layers>=graph.n || !config.resource_only))
        throw std::runtime_error("N50 is pre-freeze resource-only and must stop before final layer");
    for (const char* suffix:{".json",".jsonl",".csv"}) {
        std::ifstream existing(config.prefix+suffix);
        if (existing.good()) throw std::runtime_error("refusing to overwrite an output");
    }
    std::ofstream trace(config.prefix+".jsonl");
    if (!trace) throw std::runtime_error("cannot open layer trace");
    std::signal(SIGXCPU,on_cpu);
    rlimit limit{};
    limit.rlim_max=static_cast<rlim_t>(std::min(600.0,std::ceil(config.seconds)+5));
    limit.rlim_cur=static_cast<rlim_t>(std::min(double(limit.rlim_max),std::ceil(config.seconds)+2));
    if (setrlimit(RLIMIT_CPU,&limit)) throw std::runtime_error("cannot set hard CPU limit");
    const double start=cpu_seconds();
    const double actual_cpu_gate=std::min(config.seconds,double(limit.rlim_max)-3);
    const auto wall_start=std::chrono::steady_clock::now();
    auto guard=[&]() {
        if (cpu_signal || cpu_seconds()-start>actual_cpu_gate) throw Limit("cpu_gate");
        if (peak_rss_mib()>config.rss) throw Limit("rss_gate");
    };
    Table states, following;
    states.max_load_factor(0.8f); following.max_load_factor(0.8f);
    states.emplace(encode(State{},{}),Polynomial{0,{{1,2*graph.n+1}}});
    std::size_t max_states=1, max_key_bytes=1, transitions=0, max_cells=1, max_value_bytes=sizeof(Coefficient);
    int completed=0, building=0;
    std::string stop;
    try {
        for (int v=0;v<max_layers;++v) {
            building=v+1;
            guard();
            following.clear();
            // Reserve no more than a single layer's plausible first branching;
            // RSS is checked both before and after potentially large rehashes.
            std::size_t reserve=std::min(config.cap,std::max<std::size_t>(16,states.size()*2));
            following.reserve(reserve);
            guard();
            const Layer& layer=graph.layers[v];
            for (const auto& item:states) {
                State original=decode(item.first,layer.before);
                for (int occupied=0;occupied<=1;++occupied) {
                    State next=original;
                    int ds=advance(next,layer,bool(occupied));
                    std::string key=encode(next,layer.retained);
                    max_key_bytes=std::max(max_key_bytes,key.size());
                    auto inserted=following.try_emplace(std::move(key),Polynomial{});
                    add_shifted(inserted.first->second,item.second,occupied,ds);
                    ++transitions;
                    if (following.size()>config.cap) throw Limit("state_gate");
                    if ((transitions&4095)==0) guard();
                }
            }
            TableStats stats=summarize(following);
            if (stats.mass!=(I128(1)<<(v+1))) throw std::runtime_error("prefix mass != 2^layer");
            states.swap(following);
            completed=v+1;
            max_states=std::max(max_states,states.size());
            max_cells=std::max(max_cells,stats.nonzero_cells);
            max_value_bytes=std::max(max_value_bytes,stats.capacity_cells*sizeof(Coefficient));
            trace<<std::setprecision(12)<<"{\"layer\":"<<completed<<",\"frontier\":"<<layer.after
                 <<",\"states\":"<<states.size()<<",\"boundarykeys\":"<<states.size()
                 <<",\"nonzero_K_cells\":"<<stats.nonzero_cells<<",\"interval_K_cells\":"<<stats.interval_cells
                 <<",\"valuebytes\":"<<stats.interval_cells*sizeof(Coefficient)
                 <<",\"capacity_valuebytes\":"<<stats.capacity_cells*sizeof(Coefficient)
                 <<",\"prefix_count\":\""<<decimal(stats.mass)
                 <<"\",\"cpu_seconds\":"<<cpu_seconds()-start<<",\"peak_rss_mib\":"<<peak_rss_mib()<<"}\n";
            trace.flush();
            guard();
        }
        if (completed<graph.n) stop="pre_final_layer_gate";
    } catch (const Limit& error) {stop=error.what();}
      catch (const std::bad_alloc&) {stop="allocation_failed";}
    bool complete=(completed==graph.n && stop.empty());
    bool binomial=false;
    I128 total=0;
    if (complete) {
        std::array<std::array<HistogramValue,3>,65> hist{};
        std::array<I128,65> totals{};
        for (const auto& item:states) {
            State s=decode(item.first,0);
            for (std::size_t i=0;i<item.second.coefficients.size();++i) {
                int k=item.second.low+static_cast<int>(i);
                const Coefficient& cell=item.second.coefficients[i];
                hist[k][s.r].count+=cell.count;
                hist[k][s.r].sum_s+=cell.sum_s;
                totals[k]+=cell.count;
                total+=cell.count;
            }
        }
        I128 choose=1;
        for (int k=0;k<=graph.n;++k) {
            if (totals[k]!=choose) throw std::runtime_error("K marginal != binomial(N,K)");
            if (k<graph.n) choose=choose*(graph.n-k)/(k+1);
        }
        if (total!=(I128(1)<<graph.n)) throw std::runtime_error("terminal mass != 2^N");
        binomial=true;
        if (!config.resource_only) {
            std::ofstream csv(config.prefix+".csv");
            csv<<"K,q,count,sum_S\n";
            for (int k=0;k<=graph.n;++k) for (int r=0;r<3;++r) if(hist[k][r].count)
                csv<<k<<','<<r-1<<','<<decimal(hist[k][r].count)<<','<<decimal(hist[k][r].sum_s)<<'\n';
            if (!csv) throw std::runtime_error("histogram write failed");
        }
    }
    std::ofstream out(config.prefix+".json");
    TableStats last_stats=summarize(states), following_stats=summarize(following);
    out<<std::setprecision(12)<<"{\n\"geometry\":["<<graph.a<<','<<graph.b<<"],\n\"N\":"<<graph.n
       <<",\n\"complete\":"<<(complete?"true":"false")<<",\n\"stop\":\""<<stop
       <<"\",\n\"completed_layers\":"<<completed<<",\n\"building_layer\":"<<building
       <<",\n\"last_complete_states\":"<<states.size()<<",\n\"partial_or_previous_map_states\":"<<following.size()
       <<",\n\"maximum_complete_layer_states\":"<<max_states<<",\n\"maximum_packed_key_bytes\":"<<max_key_bytes
       <<",\n\"last_complete_nonzero_K_cells\":"<<last_stats.nonzero_cells
       <<",\n\"last_complete_valuebytes\":"<<last_stats.interval_cells*sizeof(Coefficient)
       <<",\n\"maximum_complete_nonzero_K_cells\":"<<max_cells
       <<",\n\"maximum_complete_capacity_valuebytes\":"<<max_value_bytes
       <<",\n\"partial_or_previous_nonzero_K_cells\":"<<following_stats.nonzero_cells
       <<",\n\"partial_or_previous_valuebytes\":"<<following_stats.interval_cells*sizeof(Coefficient)
       <<",\n\"transitions\":"<<transitions<<",\n\"cpu_seconds_to_receipt\":"<<cpu_seconds()-start
       <<",\n\"wall_seconds_to_receipt\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_start).count()
       <<",\n\"peak_rss_mib\":"<<peak_rss_mib()<<",\n\"binomial_and_total_checks\":"<<(binomial?"true":"false")
       <<",\n\"resource_only\":"<<(config.resource_only?"true":"false")
       <<",\n\"authorization_commit_from_verified_driver\":\""<<config.authorization_commit<<"\""
       <<",\n\"CPU_requested_gate\":"<<config.seconds<<",\n\"CPU_soft_gate\":"<<actual_cpu_gate<<",\n\"CPU_OS_hard_limit\":"<<limit.rlim_max
       <<",\n\"RSS_gate_mib\":"<<config.rss<<",\n\"state_gate\":"<<config.cap
       <<",\n\"max_layers_gate\":"<<max_layers
       <<",\n\"value_representation\":\"dynamic contiguous K support with signed64 count/sumS and checked signed128 arithmetic\",\n\"gain_representation\":\"range-checked signed16, rank0 xy, rank1 perpendicular, rank2 absent\"\n}\n";
    out.flush(); trace.flush();
    if (!out || !trace) throw std::runtime_error("receipt write failed");
    std::cout<<"geometry="<<graph.a<<','<<graph.b<<" completed_layers="<<completed<<" stop="<<stop
             <<" states="<<states.size()<<" peak_rss_mib="<<peak_rss_mib()<<'\n';
    return 0;
}

int main(int argc,char**argv) {
    try { return run(arguments(argc,argv)); }
    catch (const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
