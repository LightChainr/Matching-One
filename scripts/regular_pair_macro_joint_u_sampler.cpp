// Frozen-pilot collector for the canonical macro-window joint-U source.
//
// One invocation handles both same-N Gaussian geometries. It reuses one
// Bernoulli occupation vector and one list of 16 distinct anchor indices for
// the pair, then exhausts the complete frozen displacement window at every
// anchor. It emits batch-by-K integer moments; roots and U are scored later.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
using i128 = __int128_t;
using Kernel = std::vector<std::int16_t>;

constexpr std::uint64_t bernoulli_threshold = 10934234699625173385ULL;
constexpr std::uint64_t n100_occupation_seed = 2026090101001001ULL;
constexpr std::uint64_t n100_anchor_seed = 2026090101001002ULL;
constexpr std::uint64_t n400_occupation_seed = 2026090101004001ULL;
constexpr std::uint64_t n400_anchor_seed = 2026090101004002ULL;
constexpr std::uint64_t frozen_batches = 100;
constexpr std::uint64_t frozen_configurations_per_batch = 500;
constexpr std::uint64_t frozen_anchors = 16;

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    __uint128_t magnitude = negative ? static_cast<__uint128_t>(-(value+1))+1
                                     : static_cast<__uint128_t>(value);
    std::string result;
    while (magnitude) {
        result.push_back(static_cast<char>('0'+magnitude%10));
        magnitude /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::uint64_t unsigned_integer(const std::string& text) {
    if (text.empty() || text.front() == '-')
        throw std::invalid_argument("invalid unsigned integer: "+text);
    std::size_t end = 0;
    const auto result = std::stoull(text, &end);
    if (end != text.size()) throw std::invalid_argument("invalid unsigned integer: "+text);
    return result;
}

struct Options {
    int n = 0;
    std::uint64_t occupation_seed = 0, anchor_seed = 0;
    std::uint64_t batches = 0, configurations_per_batch = 0, anchors = 0;
    std::string kernel_path, output_path;
};

Options parse_options(int argc, char** argv) {
    if (argc != 9)
        throw std::invalid_argument(
            "usage: regular_pair_macro_joint_u_sampler N occupation_seed anchor_seed "
            "batches configurations_per_batch anchors kernel.tsv output.csv");
    Options o;
    o.n = static_cast<int>(unsigned_integer(argv[1]));
    o.occupation_seed = unsigned_integer(argv[2]);
    o.anchor_seed = unsigned_integer(argv[3]);
    o.batches = unsigned_integer(argv[4]);
    o.configurations_per_batch = unsigned_integer(argv[5]);
    o.anchors = unsigned_integer(argv[6]);
    o.kernel_path = argv[7];
    o.output_path = argv[8];
    if (o.n != 100 && o.n != 400) throw std::invalid_argument("frozen pilot accepts N=100 or N=400");
    const auto expected_occupation = o.n == 100 ? n100_occupation_seed : n400_occupation_seed;
    const auto expected_anchor = o.n == 100 ? n100_anchor_seed : n400_anchor_seed;
    if (o.occupation_seed != expected_occupation || o.anchor_seed != expected_anchor ||
        o.batches != frozen_batches || o.configurations_per_batch != frozen_configurations_per_batch ||
        o.anchors != frozen_anchors)
        throw std::invalid_argument("arguments do not match frozen pilot seeds/counts");
    return o;
}

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

Kernel read_kernel(const std::string& path, std::size_t& row_count) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read kernel: "+path);
    Kernel kernel(1U<<24, 0); // Valid omitted keys are exact zero.
    std::unordered_set<std::uint32_t> seen;
    std::string line;
    int key_column = -1, g_column = -1;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tsv(line);
        if (key_column < 0) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "key" || fields[i] == "packed_key") key_column = static_cast<int>(i);
                if (fields[i] == "g16") g_column = static_cast<int>(i);
            }
            if (key_column < 0 || g_column < 0)
                throw std::runtime_error("kernel TSV needs key (or packed_key) and g16 columns");
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max(key_column, g_column)))
            throw std::runtime_error("incomplete kernel row");
        const auto wide_key = unsigned_integer(fields[key_column]);
        std::size_t end = 0;
        const auto g = std::stoll(fields[g_column], &end);
        if (end != fields[g_column].size() || wide_key >= kernel.size() ||
            g < std::numeric_limits<std::int16_t>::min() ||
            g > std::numeric_limits<std::int16_t>::max())
            throw std::runtime_error("invalid kernel entry");
        const auto key = static_cast<std::uint32_t>(wide_key);
        int max_label = -1;
        for (int i = 0; i < 8; ++i) {
            const int label = (key>>(3*i))&7;
            if (label > max_label+1) throw std::runtime_error("noncanonical kernel key");
            max_label = std::max(max_label, label);
        }
        if (!seen.insert(key).second) throw std::runtime_error("duplicate kernel key");
        kernel[key] = static_cast<std::int16_t>(g);
    }
    if (!in.eof() || seen.empty()) throw std::runtime_error("empty or failed kernel read");
    row_count = seen.size();
    return kernel;
}

class Components {
    std::vector<int> parent, size;
public:
    explicit Components(int n) : parent(n), size(n) {}
    void reset() { std::iota(parent.begin(), parent.end(), 0); std::fill(size.begin(), size.end(), 1); }
    int root(int v) {
        while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
        return v;
    }
    void join(int a, int b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a, b);
        parent[b] = a; size[a] += size[b];
    }
};

std::uint64_t bounded_random(std::mt19937_64& rng, std::uint64_t bound) {
    if (!bound) throw std::logic_error("zero random bound");
    const std::uint64_t rejection = -bound%bound;
    std::uint64_t word;
    do word = rng(); while (word < rejection);
    return word%bound;
}

struct WindowEntry { int vertex = -1, r2 = 0, dx = 0, dy = 0; };
struct GroupValue { std::int64_t b16 = 0; std::uint64_t eligible = 0, nonzero = 0; };
struct ConfigurationValue {
    int q = 0, e = 0;
    std::array<GroupValue, 3> group; // total, exactly s=2, s>=3
    std::uint64_t s_le1_pairs = 0, s_le1_nonzero = 0;
};

class Geometry {
    int a, b, n;
    std::string name;
    const Kernel& kernel;
    std::vector<std::pair<int,int>> representatives;
    std::vector<int> coordinate_index;
    std::vector<std::array<int,4>> neighbors, port_edges; // N,E,S,W
    std::vector<WindowEntry> window;
    std::vector<unsigned char> occupied;
    std::vector<int> occupied_root;
    Components black, white;

    static int mod(int x, int modulus) { x%=modulus; return x<0 ? x+modulus : x; }
    int quotient_key(int x, int y) const { return n*mod(a*x+b*y,n)+mod(-b*x+a*y,n); }
    int translate(int x, int displacement) const {
        const auto [xx,xy] = representatives[x];
        const auto [dx,dy] = representatives[displacement];
        return coordinate_index[quotient_key(xx+dx,xy+dy)];
    }
    static std::int64_t floor_div(std::int64_t p, std::int64_t q) {
        std::int64_t result=p/q, remainder=p%q;
        if (remainder<0) --result;
        return result;
    }
    WindowEntry shortest_entry(int vertex) const {
        const auto [x,y]=representatives[vertex];
        const std::int64_t du=static_cast<std::int64_t>(a)*x+b*y;
        const std::int64_t dv=-static_cast<std::int64_t>(b)*x+a*y;
        const auto mf=floor_div(du,n), nf=floor_div(dv,n);
        WindowEntry best; best.vertex=vertex; best.r2=std::numeric_limits<int>::max();
        for (const auto m:std::array<std::int64_t,2>{mf,mf+1})
            for (const auto k:std::array<std::int64_t,2>{nf,nf+1}) {
                const auto dx=x-m*a+k*b, dy=y-m*b-k*a;
                const auto r2=dx*dx+dy*dy;
                if (r2<best.r2 || (r2==best.r2 &&
                    std::pair<std::int64_t,std::int64_t>{dx,dy}<
                    std::pair<std::int64_t,std::int64_t>{best.dx,best.dy})) {
                    best.r2=static_cast<int>(r2); best.dx=static_cast<int>(dx); best.dy=static_cast<int>(dy);
                }
            }
        return best;
    }
    int component_count(Components& components, bool black_sites) {
        int result=0;
        for (int v=0;v<n;++v)
            if ((occupied[v]!=0)==black_sites && components.root(v)==v) ++result;
        return result;
    }
    std::array<int,2> topology() {
        black.reset(); white.reset();
        int k=0,edges=0,faces=0;
        for (int v=0;v<n;++v) {
            k+=occupied[v]!=0;
            if (occupied[v]) {
                for (int d:{0,1}) if (occupied[neighbors[v][d]]) {
                    black.join(v,neighbors[v][d]); ++edges;
                }
                const int east=neighbors[v][1], north=neighbors[v][0];
                faces+=occupied[east] && occupied[north] && occupied[neighbors[east][0]];
            } else {
                for (int d:{0,1}) if (!occupied[neighbors[v][d]]) white.join(v,neighbors[v][d]);
                const int northeast=neighbors[neighbors[v][1]][0];
                const int northwest=neighbors[neighbors[v][3]][0];
                if (!occupied[northeast]) white.join(v,northeast);
                if (!occupied[northwest]) white.join(v,northwest);
            }
        }
        const int q=component_count(black,true)-component_count(white,false)-(k-edges+faces);
        if (q < -1 || q > 1) throw std::logic_error("digital-Alexander q outside {-1,0,1}");
        for (int v=0;v<n;++v) if (occupied[v]) occupied_root[v]=black.root(v);
        return {q,q*q};
    }
    int outside_id(int center,int direction) const {
        const int u=neighbors[center][direction];
        return occupied[u] ? occupied_root[u] : n+port_edges[center][direction];
    }
    std::pair<std::int64_t,int> pair_kernel(int x,int y) const {
        if (occupied[x] || occupied[y]) return {0,-1};
        std::array<int,8> outside{},labels{};
        std::uint32_t key=0; int next=0;
        for (int i=0;i<8;++i) {
            outside[i]=outside_id(i<4?x:y,i%4);
            int previous=0;
            while (previous<i && outside[previous]!=outside[i]) ++previous;
            labels[i]=previous<i ? labels[previous] : next++;
            key|=static_cast<std::uint32_t>(labels[i])<<(3*i);
        }
        std::array<bool,8> left{},right{};
        for (int i=0;i<4;++i) left[labels[i]]=true;
        for (int i=4;i<8;++i) right[labels[i]]=true;
        int shared=0;
        for (int label=0;label<next;++label) shared+=left[label]&&right[label];
        const auto value=kernel[key];
        if (shared<=1 && value!=0) throw std::logic_error("s<=1 nonzero kernel control failed");
        return {value,shared};
    }

public:
    Geometry(int aa,int bb,std::string label,const Kernel& table)
        : a(aa),b(bb),n(a*a+b*b),name(std::move(label)),kernel(table),
          coordinate_index(n*n,-1),neighbors(n),port_edges(n),occupied(n),occupied_root(n),
          black(n),white(n) {
        representatives.emplace_back(0,0); coordinate_index[quotient_key(0,0)]=0;
        for (std::size_t i=0;i<representatives.size();++i) {
            const auto [x,y]=representatives[i];
            for (const auto step:std::array<std::pair<int,int>,2>{{{1,0},{0,1}}}) {
                const int xx=x+step.first,yy=y+step.second,key=quotient_key(xx,yy);
                if (coordinate_index[key]<0) {
                    coordinate_index[key]=static_cast<int>(representatives.size());
                    representatives.emplace_back(xx,yy);
                }
            }
        }
        if (static_cast<int>(representatives.size())!=n) throw std::logic_error("wrong quotient area");
        for (int v=0;v<n;++v) {
            const auto [x,y]=representatives[v];
            neighbors[v]={coordinate_index[quotient_key(x,y+1)],coordinate_index[quotient_key(x+1,y)],
                          coordinate_index[quotient_key(x,y-1)],coordinate_index[quotient_key(x-1,y)]};
            port_edges[v]={2*v,2*v+1,2*neighbors[v][2],2*neighbors[v][3]+1};
        }
        for (int v=1;v<n;++v) {
            const auto entry=shortest_entry(v);
            if (16*entry.r2>=n && 25*entry.r2<=4*n) window.push_back(entry);
        }
        std::sort(window.begin(),window.end(),[](const auto& x,const auto& y){
            return std::tuple<int,int,int>{x.r2,x.dx,x.dy}<std::tuple<int,int,int>{y.r2,y.dx,y.dy};
        });
        if (window.empty()) throw std::logic_error("frozen window has no displacement");
        for (const auto& entry:window) for (int ex:port_edges[0]) for (int ey:port_edges[entry.vertex])
            if (ex==ey) throw std::logic_error("accepted displacement has overlapping edge ports");
    }
    const std::string& geometry_name() const { return name; }
    int geometry_a() const { return a; }
    int geometry_b() const { return b; }
    std::size_t displacement_count() const { return window.size(); }
    void write_window(const std::string& path) const {
        if (std::ifstream(path).good()) throw std::runtime_error("window table exists: "+path);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot create window table: "+path);
        out<<"ordinal,vertex,r2,canonical_dx,canonical_dy\n";
        for (std::size_t i=0;i<window.size();++i)
            out<<i<<','<<window[i].vertex<<','<<window[i].r2<<','<<window[i].dx<<','<<window[i].dy<<'\n';
        out.close(); if (!out) throw std::runtime_error("window table write failed");
    }
    ConfigurationValue evaluate(const std::vector<unsigned char>& bits,const std::vector<int>& anchors) {
        occupied=bits;
        ConfigurationValue result;
        const auto [q,e]=topology(); result.q=q; result.e=e;
        for (int x:anchors) for (const auto& displacement:window) {
            const int y=translate(x,displacement.vertex);
            const auto [g16,shared]=pair_kernel(x,y);
            if (shared<0) continue; // At least one endpoint occupied.
            ++result.group[0].eligible;
            result.group[0].nonzero+=g16!=0;
            result.group[0].b16+=g16;
            if (shared<=1) {
                ++result.s_le1_pairs; result.s_le1_nonzero+=g16!=0;
            } else {
                const int group=shared==2?1:2;
                ++result.group[group].eligible;
                result.group[group].nonzero+=g16!=0;
                result.group[group].b16+=g16;
            }
        }
        if (result.group[0].b16!=result.group[1].b16+result.group[2].b16 ||
            result.s_le1_nonzero!=0)
            throw std::logic_error("configuration support/additivity control failed");
        return result;
    }
};

struct GroupSums {
    i128 b16=0,qb16=0,eb16=0;
    std::uint64_t eligible=0,nonzero=0;
};
struct Cell {
    std::uint64_t count=0;
    i128 q=0,e=0;
    std::array<GroupSums,3> group;
    std::uint64_t s_le1_pairs=0,s_le1_nonzero=0;
    void add(const ConfigurationValue& value) {
        ++count; q+=value.q; e+=value.e;
        for (int index=0;index<3;++index) {
            const auto& source=value.group[index]; auto& target=group[index];
            target.b16+=source.b16; target.qb16+=static_cast<i128>(value.q)*source.b16;
            target.eb16+=static_cast<i128>(value.e)*source.b16;
            target.eligible+=source.eligible; target.nonzero+=source.nonzero;
        }
        s_le1_pairs+=value.s_le1_pairs; s_le1_nonzero+=value.s_le1_nonzero;
    }
};

std::vector<int> sample_anchors(std::mt19937_64& rng,int n,std::uint64_t count) {
    std::vector<int> pool(n); std::iota(pool.begin(),pool.end(),0);
    for (std::uint64_t i=0;i<count;++i) {
        const auto j=i+bounded_random(rng,n-i);
        std::swap(pool[i],pool[j]);
    }
    pool.resize(count);
    return pool;
}
void write_group(std::ostream& out,const GroupSums& group) {
    out<<','<<decimal(group.b16)<<','<<decimal(group.qb16)<<','<<decimal(group.eb16)
       <<','<<group.eligible<<','<<group.nonzero;
}
} // namespace

int main(int argc,char** argv) {
    try {
        const auto start=std::chrono::steady_clock::now();
        const auto o=parse_options(argc,argv);
        std::size_t kernel_rows=0;
        const auto kernel=read_kernel(o.kernel_path,kernel_rows);
        const std::array<int,2> axis=o.n==100?std::array<int,2>{10,0}:std::array<int,2>{20,0};
        const std::array<int,2> tilted=o.n==100?std::array<int,2>{8,6}:std::array<int,2>{16,12};
        Geometry axis_geometry(axis[0],axis[1],"axis",kernel);
        Geometry tilted_geometry(tilted[0],tilted[1],"tilted",kernel);
        const std::string axis_window=o.output_path+".axis.window.csv";
        const std::string tilted_window=o.output_path+".tilted.window.csv";
        if (std::ifstream(o.output_path).good() || std::ifstream(axis_window).good() ||
            std::ifstream(tilted_window).good())
            throw std::runtime_error("output or window sidecar exists");
        axis_geometry.write_window(axis_window);
        tilted_geometry.write_window(tilted_window);
        std::ofstream out(o.output_path);
        if (!out) throw std::runtime_error("cannot create output");
        out<<"N,batch,geometry,a,b,K,count,sum_q,sum_E,occupation_seed,anchor_seed,"
              "anchors_per_configuration,window_displacements,source_denominator";
        for (const std::string group:{"total","s2","sge3"})
            out<<",sum_B16_"<<group<<",sum_qB16_"<<group<<",sum_EB16_"<<group
               <<",eligible_pair_count_"<<group<<",nonzero_pair_count_"<<group;
        out<<",s_le1_pair_count,s_le1_nonzero_g16_count,total_minus_s2_minus_sge3_B16\n";

        std::mt19937_64 occupation_rng(o.occupation_seed),anchor_rng(o.anchor_seed);
        std::vector<unsigned char> bits(o.n);
        for (std::uint64_t batch=0;batch<o.batches;++batch) {
            std::array<std::vector<Cell>,2> histogram{
                std::vector<Cell>(o.n+1),std::vector<Cell>(o.n+1)};
            for (std::uint64_t configuration=0;configuration<o.configurations_per_batch;++configuration) {
                int k=0;
                for (int v=0;v<o.n;++v) { bits[v]=occupation_rng()<bernoulli_threshold; k+=bits[v]!=0; }
                const auto anchors=sample_anchors(anchor_rng,o.n,o.anchors);
                histogram[0][k].add(axis_geometry.evaluate(bits,anchors));
                histogram[1][k].add(tilted_geometry.evaluate(bits,anchors));
            }
            const std::array<Geometry*,2> geometries{&axis_geometry,&tilted_geometry};
            for (int geometry_index=0;geometry_index<2;++geometry_index) for (int k=0;k<=o.n;++k) {
                const auto& geometry=*geometries[geometry_index];
                const auto& cell=histogram[geometry_index][k];
                const i128 denominator=static_cast<i128>(16)*o.anchors*o.n;
                const i128 residual=cell.group[0].b16-cell.group[1].b16-cell.group[2].b16;
                if (residual!=0 || cell.s_le1_nonzero!=0)
                    throw std::logic_error("batch/K support/additivity control failed");
                out<<o.n<<','<<batch<<','<<geometry.geometry_name()<<','<<geometry.geometry_a()<<','
                   <<geometry.geometry_b()<<','<<k<<','<<cell.count<<','<<decimal(cell.q)<<','<<decimal(cell.e)
                   <<','<<o.occupation_seed<<','<<o.anchor_seed<<','<<o.anchors<<','
                   <<geometry.displacement_count()<<','<<decimal(denominator);
                for (const auto& group:cell.group) write_group(out,group);
                out<<','<<cell.s_le1_pairs<<','<<cell.s_le1_nonzero<<','<<decimal(residual)<<'\n';
            }
        }
        out.close(); if (!out) throw std::runtime_error("output write failed");
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<std::setprecision(17)<<"{\"status\":\"completed\",\"N\":"<<o.n
                 <<",\"kernel_rows\":"<<kernel_rows<<",\"axis_window\":"
                 <<axis_geometry.displacement_count()<<",\"tilted_window\":"
                 <<tilted_geometry.displacement_count()<<",\"configurations\":"
                 <<o.batches*o.configurations_per_batch<<",\"elapsed_seconds\":"<<seconds<<"}\n";
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
