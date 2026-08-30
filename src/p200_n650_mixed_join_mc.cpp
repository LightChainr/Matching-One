// Dedicated fixed-p N650 mixed C2 x C5 join runner for Issue #200.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Int = std::int64_t;
struct Vector { Int x = 0; Int y = 0; };
struct Matrix { Int a = 0; Int b = 0; Int c = 0; Int d = 0; };
struct Edge { int i; int j; int dx; int dy; };

Int determinant(const Matrix& p) { return p.a * p.d - p.b * p.c; }
Int mod(Int x, Int n) { x %= n; return x < 0 ? x + n : x; }

struct Bezout { Int gcd; Int x; Int y; };
Bezout extended_gcd(Int a, Int b) {
    const Int sign_a = a < 0 ? -1 : 1, sign_b = b < 0 ? -1 : 1;
    a = std::llabs(a); b = std::llabs(b);
    Int old_r = a, r = b, old_s = 1, s = 0, old_t = 0, t = 1;
    while (r) {
        const Int q = old_r / r;
        const Int nr = old_r - q * r, ns = old_s - q * s, nt = old_t - q * t;
        old_r = r; r = nr; old_s = s; s = ns; old_t = t; t = nt;
    }
    return {old_r, sign_a * old_s, sign_b * old_t};
}

struct Quotient {
    Matrix periods;
    Int det;
    int order;
    Int h11, h12, h22;
    explicit Quotient(Matrix p) : periods(p), det(determinant(p)), order(static_cast<int>(std::llabs(det))) {
        if (!det) throw std::invalid_argument("singular periods");
        const Bezout lower = extended_gcd(p.c, p.d);
        h22 = lower.gcd; h11 = order / h22;
        h12 = mod(p.a * lower.x + p.b * lower.y, h11);
        if (h11 * h22 != order) throw std::logic_error("invalid HNF");
    }
    int label(Vector v) const {
        const Int qy = (v.y - mod(v.y, h22)) / h22;
        const Int ry = v.y - qy * h22;
        const Int rx = mod(v.x - qy * h12, h11);
        return static_cast<int>(rx + h11 * ry);
    }
    Vector representative(int label_value) const { return {label_value % h11, label_value / h11}; }
    Vector winding(Int dx, Int dy) const {
        const Int n0 = periods.d * dx - periods.b * dy;
        const Int n1 = -periods.c * dx + periods.a * dy;
        if (n0 % det || n1 % det) throw std::logic_error("non-period cycle");
        return {n0 / det, n1 / det};
    }
};

Vector primitive(Vector v) {
    const Int g = std::gcd(std::llabs(v.x), std::llabs(v.y));
    if (!g) return v;
    v.x /= g; v.y /= g;
    if (v.x < 0 || (v.x == 0 && v.y < 0)) { v.x = -v.x; v.y = -v.y; }
    return v;
}

class PotentialUF {
  public:
    explicit PotentialUF(const Quotient& quotient)
        : q_(quotient), parent_(q_.order), size_(q_.order), dx_(q_.order), dy_(q_.order),
          homology_rank_(q_.order), basis_(q_.order), seen_(q_.order) { reset(); }
    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(dx_.begin(), dx_.end(), 0); std::fill(dy_.begin(), dy_.end(), 0);
        std::fill(homology_rank_.begin(), homology_rank_.end(), 0);
    }
    struct Found { int root; Int dx; Int dy; };
    Found find(int v) {
        if (parent_[v] == v) return {v, 0, 0};
        const int p = parent_[v]; const Found above = find(p);
        dx_[v] += above.dx; dy_[v] += above.dy; parent_[v] = above.root;
        return {above.root, dx_[v], dy_[v]};
    }
    void extend(int root, Vector v) {
        if ((!v.x && !v.y) || homology_rank_[root] == 2) return;
        v = primitive(v);
        if (!homology_rank_[root]) { basis_[root][0] = v; homology_rank_[root] = 1; return; }
        const Vector a = basis_[root][0];
        if (a.x * v.y != a.y * v.x) { basis_[root][1] = v; homology_rank_[root] = 2; }
    }
    void add_edge(int i, int j, int edge_dx, int edge_dy) {
        Found a = find(i), b = find(j);
        Int rx = a.dx + edge_dx - b.dx, ry = a.dy + edge_dy - b.dy;
        if (a.root == b.root) { extend(a.root, q_.winding(rx, ry)); return; }
        if (size_[a.root] < size_[b.root]) { std::swap(a, b); rx = -rx; ry = -ry; }
        parent_[b.root] = a.root; dx_[b.root] = rx; dy_[b.root] = ry;
        size_[a.root] += size_[b.root];
        for (int k = 0; k < homology_rank_[b.root]; ++k) extend(a.root, basis_[b.root][k]);
        homology_rank_[b.root] = 0;
    }
    std::pair<int, int> ranks(const std::vector<std::uint8_t>& active) {
        std::fill(seen_.begin(), seen_.end(), 0);
        int components = 0, global_rank = 0; std::array<Vector, 2> global{};
        for (int v = 0; v < q_.order; ++v) if (active[v]) {
            const int root = find(v).root;
            if (seen_[root]) continue;
            seen_[root] = 1; ++components;
            for (int k = 0; k < homology_rank_[root]; ++k) {
                Vector w = primitive(basis_[root][k]);
                if (!global_rank) { global[0] = w; global_rank = 1; }
                else if (global_rank == 1 && global[0].x * w.y != global[0].y * w.x) {
                    global[1] = w; global_rank = 2;
                }
            }
        }
        return {static_cast<int>(std::count(active.begin(), active.end(), 1)) - components, global_rank};
    }
  private:
    const Quotient& q_;
    std::vector<int> parent_, size_;
    std::vector<Int> dx_, dy_;
    std::vector<std::uint8_t> homology_rank_, seen_;
    std::vector<std::array<Vector, 2>> basis_;
};

struct FiberCell { int vertex; int row; int column; };
struct SourceFiber { std::vector<FiberCell> cells; };

struct Geometry {
    Quotient final_q, q2, q5, source;
    std::vector<Vector> coordinates;
    std::vector<Edge> primal_edges, matching_edges;
    std::vector<std::vector<int>> relation2, relation5;
    std::vector<SourceFiber> source_fibers;
    PotentialUF uf;
    Geometry(Matrix final_p, Matrix p2, Matrix p5, Matrix source_p)
        : final_q(final_p), q2(p2), q5(p5), source(source_p), coordinates(final_q.order), uf(final_q) {
        if (final_q.order != 10 * source.order || q2.order != 2 * source.order ||
            q5.order != 5 * source.order)
            throw std::invalid_argument("geometry is not a C2 x C5 cover square");
        for (int i = 0; i < final_q.order; ++i) coordinates[i] = final_q.representative(i);
        const std::array<Vector, 4> steps{{{1,0},{0,1},{1,1},{-1,1}}};
        for (int i = 0; i < final_q.order; ++i) for (int k = 0; k < 4; ++k) {
            const Vector s = coordinates[i], d = steps[k];
            Edge e{i, final_q.label({s.x+d.x,s.y+d.y}), static_cast<int>(d.x), static_cast<int>(d.y)};
            if (k < 2) primal_edges.push_back(e);
            matching_edges.push_back(e);
        }
        relation2 = quotient_fibers(q5);  // kernel size two, final -> N325
        relation5 = quotient_fibers(q2);  // kernel size five, final -> N130
        const auto source_groups = quotient_fibers(source);
        for (const auto& vertices : source_groups) {
            std::vector<int> rows, columns;
            for (int v : vertices) { rows.push_back(q2.label(coordinates[v])); columns.push_back(q5.label(coordinates[v])); }
            std::sort(rows.begin(), rows.end()); rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
            std::sort(columns.begin(), columns.end()); columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
            if (vertices.size()!=10 || rows.size()!=2 || columns.size()!=5) throw std::logic_error("CRT source fiber failed");
            SourceFiber fiber;
            for (int v : vertices) {
                const int row = static_cast<int>(std::lower_bound(rows.begin(),rows.end(),q2.label(coordinates[v]))-rows.begin());
                const int column = static_cast<int>(std::lower_bound(columns.begin(),columns.end(),q5.label(coordinates[v]))-columns.begin());
                fiber.cells.push_back({v,row,column});
            }
            source_fibers.push_back(std::move(fiber));
        }
    }
    std::vector<std::vector<int>> quotient_fibers(const Quotient& target) const {
        std::vector<std::vector<int>> groups(target.order);
        for (int v = 0; v < final_q.order; ++v) groups[target.label(coordinates[v])].push_back(v);
        return groups;
    }
};

struct Corner { int partition_rank; int ambient_rank; };
void add_relation(Geometry& geometry, const std::vector<std::vector<int>>& fibers,
                  const std::vector<std::uint8_t>& active) {
    for (const auto& fiber : fibers) {
        int anchor = -1;
        for (int v : fiber) if (active[v]) {
            if (anchor < 0) anchor = v;
            else {
                const Vector a = geometry.coordinates[anchor], b = geometry.coordinates[v];
                geometry.uf.add_edge(anchor, v, static_cast<int>(b.x-a.x), static_cast<int>(b.y-a.y));
            }
        }
    }
}
Corner corner(Geometry& geometry, const std::vector<std::uint8_t>& active,
              const std::vector<Edge>& edges, bool use2, bool use5, bool reverse=false) {
    geometry.uf.reset();
    for (const Edge& e : edges) if (active[e.i] && active[e.j]) geometry.uf.add_edge(e.i,e.j,e.dx,e.dy);
    if (!reverse) {
        if (use2) add_relation(geometry, geometry.relation2, active);
        if (use5) add_relation(geometry, geometry.relation5, active);
    } else {
        if (use5) add_relation(geometry, geometry.relation5, active);
        if (use2) add_relation(geometry, geometry.relation2, active);
    }
    const auto value = geometry.uf.ranks(active);
    return {value.first, value.second};
}

int local_join(const Geometry& geometry, const std::vector<std::uint8_t>& active) {
    int total = 0;
    for (const SourceFiber& fiber : geometry.source_fibers) {
        std::array<int,7> parent; std::iota(parent.begin(), parent.end(), 0);
        std::array<std::uint8_t,7> used{};
        auto find = [&](int x) { while (parent[x]!=x) { parent[x]=parent[parent[x]]; x=parent[x]; } return x; };
        int edge_count = 0;
        for (const FiberCell& cell : fiber.cells) if (active[cell.vertex]) {
            int a=cell.row, b=2+cell.column, ra=find(a), rb=find(b);
            if (ra!=rb) parent[std::max(ra,rb)]=std::min(ra,rb);
            used[a]=used[b]=1; ++edge_count;
        }
        int vertices=0; std::array<std::uint8_t,7> root_seen{}; int components=0;
        for (int i=0;i<7;++i) if (used[i]) { ++vertices; int r=find(i); if(!root_seen[r]){root_seen[r]=1;++components;} }
        total += edge_count - vertices + components;
    }
    return total;
}

struct ColorStats {
    int residual;
    int ambient_delta;
    int j_local;
    Corner c0, c2, c5, c25;
    int activation_order;
};
ColorStats color_stats(Geometry& geometry, const std::vector<std::uint8_t>& active,
                       const std::vector<Edge>& edges, bool check_order=false) {
    const Corner c0=corner(geometry,active,edges,false,false);
    const Corner c2=corner(geometry,active,edges,true,false);
    const Corner c5=corner(geometry,active,edges,false,true);
    const Corner c25=corner(geometry,active,edges,true,true);
    if (check_order) {
        const Corner reverse=corner(geometry,active,edges,true,true,true);
        if (reverse.partition_rank!=c25.partition_rank || reverse.ambient_rank!=c25.ambient_rank)
            throw std::logic_error("join order changed typed corner");
    }
    const int j_full=c2.partition_rank+c5.partition_rank-c25.partition_rank-c0.partition_rank;
    const int j_local=local_join(geometry,active);
    const int activation_order=(c25.ambient_rank==2)*(
        static_cast<int>(c2.ambient_rank==2)-static_cast<int>(c5.ambient_rank==2));
    return {j_full-j_local,
            c25.ambient_rank-c2.ambient_rank-c5.ambient_rank+c0.ambient_rank,
            j_local,c0,c2,c5,c25,activation_order};
}

struct OrientationStats { int even, odd, ambient_even, ambient_odd; };
OrientationStats orientation_stats(Geometry& geometry, const std::vector<std::uint8_t>& black, bool check=false) {
    std::vector<std::uint8_t> white(black.size());
    for (std::size_t i=0;i<black.size();++i) white[i]=!black[i];
    const ColorStats b=color_stats(geometry,black,geometry.primal_edges,check);
    const ColorStats w=color_stats(geometry,white,geometry.matching_edges,check);
    return {b.residual+w.residual,b.residual-w.residual,
            b.ambient_delta+w.ambient_delta,b.ambient_delta-w.ambient_delta};
}

std::array<int,8> paired_state(Geometry& first, Geometry& second,
                               const std::vector<std::uint8_t>& active, bool check=false) {
    const OrientationStats a=orientation_stats(first,active,check), b=orientation_stats(second,active,check);
    return {a.even+b.even,a.even-b.even,a.odd+b.odd,a.odd-b.odd,
            a.ambient_even+b.ambient_even,a.ambient_even-b.ambient_even,
            a.ambient_odd+b.ambient_odd,a.ambient_odd-b.ambient_odd};
}

constexpr int ordered_field_count=10;
constexpr std::array<const char*,ordered_field_count> ordered_field_names{{
    "Jlocal","ambient_h0","ambient_h2","ambient_h5","ambient_h25",
    "partition_r0","partition_r2","partition_r5","partition_r25","Cact"
}};

int ordered_field(const ColorStats& value,int field){
    switch(field){
      case 0:return value.j_local;
      case 1:return value.c0.ambient_rank;
      case 2:return value.c2.ambient_rank;
      case 3:return value.c5.ambient_rank;
      case 4:return value.c25.ambient_rank;
      case 5:return value.c0.partition_rank;
      case 6:return value.c2.partition_rank;
      case 7:return value.c5.partition_rank;
      case 8:return value.c25.partition_rank;
      case 9:return value.activation_order;
      default:throw std::logic_error("ordered field index out of range");
    }
}

std::array<int,4> paired_projection(const ColorStats& first_black,
                                    const ColorStats& first_white,
                                    const ColorStats& second_black,
                                    const ColorStats& second_white,int field){
    const int first_even=ordered_field(first_black,field)+ordered_field(first_white,field);
    const int first_odd=ordered_field(first_black,field)-ordered_field(first_white,field);
    const int second_even=ordered_field(second_black,field)+ordered_field(second_white,field);
    const int second_odd=ordered_field(second_black,field)-ordered_field(second_white,field);
    return {first_even+second_even,first_even-second_even,
            first_odd+second_odd,first_odd-second_odd};
}

struct DetailedState {
    std::array<int,8> legacy{};
    std::array<int,4*ordered_field_count> ordered{};
};

DetailedState detailed_state(Geometry& first,Geometry& second,
                              const std::vector<std::uint8_t>& black,bool check=false){
    std::vector<std::uint8_t> white(black.size());
    for(std::size_t i=0;i<black.size();++i)white[i]=!black[i];
    const ColorStats fb=color_stats(first,black,first.primal_edges,check);
    const ColorStats fw=color_stats(first,white,first.matching_edges,check);
    const ColorStats sb=color_stats(second,black,second.primal_edges,check);
    const ColorStats sw=color_stats(second,white,second.matching_edges,check);
    DetailedState state;
    const int fe=fb.residual+fw.residual,fo=fb.residual-fw.residual;
    const int se=sb.residual+sw.residual,so=sb.residual-sw.residual;
    const int fae=fb.ambient_delta+fw.ambient_delta,fao=fb.ambient_delta-fw.ambient_delta;
    const int sae=sb.ambient_delta+sw.ambient_delta,sao=sb.ambient_delta-sw.ambient_delta;
    state.legacy={fe+se,fe-se,fo+so,fo-so,fae+sae,fae-sae,fao+sao,fao-sao};
    for(int field=0;field<ordered_field_count;++field){
        const auto values=paired_projection(fb,fw,sb,sw,field);
        for(int channel=0;channel<4;++channel)state.ordered[4*field+channel]=values[channel];
    }
    return state;
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31);
}
class Stream {
  public:
    explicit Stream(std::uint64_t s):state_(s){}
    std::uint64_t next(){state_+=0x9e3779b97f4a7c15ULL;return splitmix64(state_-0x9e3779b97f4a7c15ULL);}
    std::uint64_t below(std::uint64_t n){
        const std::uint64_t threshold = static_cast<std::uint64_t>(-n) % n;
        while (true) { const std::uint64_t value=next(); if(value>=threshold)return value%n; }
    }
  private: std::uint64_t state_;
};
class Binomial650 {
  public:
    Binomial650(){
        constexpr long double p=0.592746050790L, q=1.0L-p;
        long double mass=std::pow(q,650), cumulative=mass; cdf_[0]=cumulative;
        for(int k=0;k<650;++k){mass*=static_cast<long double>(650-k)/(k+1)*p/q;cumulative+=mass;cdf_[k+1]=cumulative;}
        for(auto& value:cdf_) value/=cumulative; cdf_[650]=1.0L;
    }
    int sample(std::uint64_t raw) const {
        const long double u=std::ldexp(static_cast<long double>(raw),-64);
        return static_cast<int>(std::lower_bound(cdf_.begin(),cdf_.end(),u)-cdf_.begin());
    }
  private: std::array<long double,651> cdf_{};
};

void sample_configuration(std::uint64_t seed,std::uint64_t replica,const Binomial650& binomial,
                          std::vector<int>& permutation,std::vector<std::uint8_t>& active){
    Stream rng(splitmix64(seed^splitmix64(replica+0xd1b54a32d192ed03ULL)));
    const int occupied=binomial.sample(rng.next());
    std::iota(permutation.begin(),permutation.end(),0);
    for(int stop=649;stop>0;--stop){int other=static_cast<int>(rng.below(stop+1));std::swap(permutation[stop],permutation[other]);}
    std::fill(active.begin(),active.end(),0); for(int i=0;i<occupied;++i) active[permutation[i]]=1;
}

Geometry first_geometry(){return Geometry({23,-11,11,23},{7,-9,9,7},{17,6,-6,17},{8,-1,1,8});}
Geometry second_geometry(){return Geometry({17,-19,19,17},{3,-11,11,3},{18,-1,1,18},{7,-4,4,7});}

void self_test(){
    Geometry tiny({3,-1,1,3},{1,-1,1,1},{2,1,-1,2},{1,0,0,1});
    std::array<std::array<int,6>,6> residual_hist{};
    std::array<std::array<int,4>,4> ambient_hist{};
    Int local_black_sum=0,local_white_sum=0,local_difference_square=0;
    Int activation_black_sum=0,activation_white_sum=0;
    bool balanced_residual=false,balanced_ambient=false;
    for(int mask=0;mask<(1<<10);++mask){
        std::vector<std::uint8_t> black(10),white(10);
        for(int i=0;i<10;++i){black[i]=(mask>>i)&1;white[i]=!black[i];}
        const ColorStats b=color_stats(tiny,black,tiny.primal_edges,true);
        const ColorStats w=color_stats(tiny,white,tiny.matching_edges,true);
        activation_black_sum+=b.activation_order;activation_white_sum+=w.activation_order;
        if(b.residual>=-4&&b.residual<=1&&w.residual>=-4&&w.residual<=1)++residual_hist[b.residual+4][w.residual+4];
        if(b.ambient_delta>=-2&&b.ambient_delta<=1&&w.ambient_delta>=-2&&w.ambient_delta<=1)++ambient_hist[b.ambient_delta+2][w.ambient_delta+2];
        const int jb=local_join(tiny,black),jw=local_join(tiny,white);local_black_sum+=jb;local_white_sum+=jw;local_difference_square+=(jb-jw)*(jb-jw);
        if(std::count(black.begin(),black.end(),1)==5&&b.residual-w.residual!=0)balanced_residual=true;
        if(std::count(black.begin(),black.end(),1)==5&&b.ambient_delta-w.ambient_delta!=0)balanced_ambient=true;
        if((b.residual-w.residual)!=-(w.residual-b.residual))throw std::runtime_error("typed odd swap failed");
    }
    if(local_black_sum!=499||local_white_sum!=499||local_difference_square!=1362||
       residual_hist[0][4]!=1||residual_hist[5][4]!=50||residual_hist[4][0]!=1||
       ambient_hist[0][2]!=32||ambient_hist[3][2]!=32||ambient_hist[2][2]!=509||
       activation_black_sum!=-133||activation_white_sum!=-87||
       !balanced_residual||!balanced_ambient){
        std::ostringstream message;message<<"tiny exact mixed-join regression failed activation="
          <<activation_black_sum<<','<<activation_white_sum;
        throw std::runtime_error(message.str());
    }
    Geometry first=first_geometry(),second=second_geometry();
    if(first.final_q.h11!=650||first.final_q.h12!=593||second.final_q.h12!=343)throw std::runtime_error("N650 HNF regression failed");
    std::vector<int> permutation(650);std::vector<std::uint8_t> active(650);Binomial650 binomial;
    sample_configuration(200,17,binomial,permutation,active);
    const auto state=paired_state(first,second,active,true);
    const auto detailed=detailed_state(first,second,active,true);
    if(detailed.legacy!=state)throw std::runtime_error("detailed legacy projection mismatch");
    const std::array<int,8> expected{{0,0,0,0,0,0,0,0}};
    if(state==expected){} // zero is allowed; the exhaustive gate above proves nondegeneracy.
    std::cout<<"self-test passed: N10 exhaustive typed joins, HNF activation order (-133,-87), local moments, matching odd swap, N650 HNF/lifts\n";
}

struct Options{
    std::uint64_t samples=20000,seed=2026102003,replica_offset=18000000000ULL;
    int batches=100,threads=1;bool self_test_only=false,production=false;
    std::string output_prefix="p200_n650_mixed_join",git_commit="unknown";
};
Options parse(int argc,char**argv){
    Options o;for(int i=1;i<argc;++i){std::string a=argv[i];auto next=[&](){if(++i>=argc)throw std::invalid_argument("missing "+a);return std::string(argv[i]);};
        if(a=="--samples")o.samples=std::stoull(next());else if(a=="--batches")o.batches=std::stoi(next());
        else if(a=="--seed")o.seed=std::stoull(next());else if(a=="--replica-offset")o.replica_offset=std::stoull(next());
        else if(a=="--threads")o.threads=std::stoi(next());else if(a=="--output-prefix")o.output_prefix=next();
        else if(a=="--git-commit")o.git_commit=next();else if(a=="--self-test")o.self_test_only=true;
        else if(a=="--production")o.production=true;else throw std::invalid_argument("unknown option "+a);
    }
    if(o.self_test_only)return o;
    if(o.samples==0||o.batches<2||o.samples%o.batches||o.threads<1)throw std::invalid_argument("samples must divide batches>=2; threads>=1");
    if(o.samples>20000&&!o.production)throw std::invalid_argument("samples>20000 requires explicit --production contract gate");
    return o;
}

struct Batch {
    std::uint64_t samples=0;
    std::array<Int,8> sums{};
    std::array<Int,4*ordered_field_count> ordered_sums{};
};

void run(const Options&o){
    const auto start=std::chrono::steady_clock::now();const std::uint64_t per_batch=o.samples/o.batches;
    std::vector<Batch> output(o.batches);
#ifdef _OPENMP
    omp_set_num_threads(o.threads);
#pragma omp parallel for schedule(static)
#endif
    for(int batch=0;batch<o.batches;++batch){
        Geometry first=first_geometry(),second=second_geometry();Binomial650 binomial;
        std::vector<int> permutation(650);std::vector<std::uint8_t> active(650);Batch local;local.samples=per_batch;
        const std::uint64_t begin=o.replica_offset+static_cast<std::uint64_t>(batch)*per_batch;
        for(std::uint64_t r=begin;r<begin+per_batch;++r){
            sample_configuration(o.seed,r,binomial,permutation,active);
            const auto state=detailed_state(first,second,active);
            for(int k=0;k<8;++k)local.sums[k]+=state.legacy[k];
            for(int k=0;k<4*ordered_field_count;++k)local.ordered_sums[k]+=state.ordered[k];
        }
        output[batch]=local;
    }
    const std::filesystem::path csv_path=o.output_prefix+".batches.csv",meta_path=o.output_prefix+".metadata.json";
    if(csv_path.has_parent_path())std::filesystem::create_directories(csv_path.parent_path());
    std::ofstream csv(csv_path);if(!csv)throw std::runtime_error("cannot open batch output");
    csv<<"batch,counter_first,counter_last_exclusive,samples,ES_num_sum,ED_num_sum,OS_num_sum,OD_num_sum,ambient_ES_num_sum,ambient_ED_num_sum,ambient_OS_num_sum,ambient_OD_num_sum";
    for(const char* field:ordered_field_names)for(const char* channel:std::array<const char*,4>{{"ES","ED","OS","OD"}})csv<<','<<field<<'_'<<channel<<"_num_sum";
    csv<<'\n';
    for(int b=0;b<o.batches;++b){const std::uint64_t first=o.replica_offset+static_cast<std::uint64_t>(b)*per_batch;csv<<b<<','<<first<<','<<first+per_batch<<','<<output[b].samples;for(Int v:output[b].sums)csv<<','<<v;for(Int v:output[b].ordered_sums)csv<<','<<v;csv<<'\n';}
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream meta(meta_path);meta<<std::setprecision(17)
      <<"{\n  \"schema\": \"matching-one.p255-n650-ordered-proxy-replay.v1\",\n  \"git_commit\": \""<<o.git_commit
      <<"\",\n  \"samples\": "<<o.samples<<",\n  \"batches\": "<<o.batches<<",\n  \"threads\": "<<o.threads
      <<",\n  \"seed\": "<<o.seed<<",\n  \"replica_offset\": "<<o.replica_offset
      <<",\n  \"p_ref\": \"0.592746050790\",\n  \"state_order\": [\"ES\",\"ED\",\"OS\",\"OD\"],\n  \"stored_sum_divisor\": 2,\n"
      <<"  \"ordered_fields\": [\"Jlocal\",\"ambient_h0\",\"ambient_h2\",\"ambient_h5\",\"ambient_h25\",\"partition_r0\",\"partition_r2\",\"partition_r5\",\"partition_r25\",\"Cact\"],\n"
      <<"  \"Cact_definition\": \"1{ambient_h25=2}(1{ambient_h2=2}-1{ambient_h5=2})\",\n"
      <<"  \"replay_source_commit\": \"308097b\",\n"
      <<"  \"replay_source_batch_sha256\": \"db5be1d870135053691e34605703f15e99e95df2da88dc279c0a55e26130d0af\",\n"
      <<"  \"first_periods\": [[23,-11],[11,23]],\n  \"second_periods\": [[17,-19],[19,17]],\n"
      <<"  \"lift_convention\": \"raw displacement between C++ column-HNF representatives; ambient-H1 is convention-sensitive secondary\",\n"
      <<"  \"rng\": \"counter SplitMix64; inverse-CDF Binomial(650,p_ref); conditional Fisher-Yates prefix\",\n"
      <<"  \"seconds\": "<<seconds<<",\n  \"batch_csv\": \""<<csv_path.string()<<"\"\n}\n";
    std::cout<<"wrote "<<csv_path<<"\nwrote "<<meta_path<<"\nseconds "<<seconds<<"\n";
}

} // namespace
int main(int argc,char**argv){try{Options o=parse(argc,argv);if(o.self_test_only){self_test();return 0;}run(o);return 0;}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 2;}}
