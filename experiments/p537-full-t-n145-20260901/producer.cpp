// N145 Monte Carlo for the complete canonical Kreg pair source.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
constexpr int N=145, BATCHES=100, KEY_SPACE=1<<24;
using Point=std::pair<int,int>;
int mod(int x){x%=N;return x<0?x+N:x;}
int qkey(int a,int b,int x,int y){return N*mod(a*x+b*y)+mod(-b*x+a*y);}
std::uint64_t mix(std::uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
double unit(std::uint64_t x){return double(mix(x)>>11)*0x1.0p-53;}

std::vector<std::int16_t> read_kernel(const std::string& path){
  std::ifstream f(path);if(!f)throw std::runtime_error("kernel open");
  std::vector<std::int16_t> out(KEY_SPACE);std::string s;int kc=-1,gc=-1;
  while(std::getline(f,s)){if(s.empty()||s[0]=='#')continue;std::stringstream z(s);std::vector<std::string> v;std::string x;while(std::getline(z,x,'\t'))v.push_back(x);
    if(kc<0){for(int i=0;i<int(v.size());++i){if(v[i]=="key"||v[i]=="packed_key")kc=i;if(v[i]=="g16")gc=i;}continue;}
    auto key=std::stoul(v.at(kc));auto value=std::stol(v.at(gc));if(key>=out.size()||value<std::numeric_limits<std::int16_t>::min()||value>std::numeric_limits<std::int16_t>::max())throw std::runtime_error("kernel range");out[key]=value;
  }if(kc<0||gc<0)throw std::runtime_error("kernel schema");return out;
}

struct DSU{std::array<int,N> p,s;DSU(){for(int i=0;i<N;++i)p[i]=i,s[i]=1;}int root(int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}void join(int x,int y){x=root(x);y=root(y);if(x==y)return;if(s[x]<s[y])std::swap(x,y);p[y]=x;s[x]+=s[y];}};
struct State{int q;std::array<int,N> black;};
struct Sum{std::uint64_t count=0;std::int64_t a0=0,qa0=0,ea0=0,a1=0,qa1=0,ea1=0,n0=0,qn0=0,en0=0,n1=0,qn1=0,en1=0;};

std::vector<Point> catalogue(){
  std::vector<std::tuple<int,int,int>> c;for(int x=-40;x<=40;++x)for(int y=-40;y<=40;++y)c.emplace_back(x*x+y*y,x,y);std::sort(c.begin(),c.end());
  std::array<bool,N*N> one{},two{};std::vector<Point> out;
  for(auto [r,x,y]:c){int u=qkey(12,1,x,y),v=qkey(9,8,x,y);if(!one[u]&&!two[v]){one[u]=two[v]=true;out.push_back({x,y});if(out.size()==N)break;}}
  if(out.size()!=N)throw std::logic_error("common catalogue");return out;
}

struct Geometry{
  int a,b;std::array<std::array<int,4>,N> nn,diag,edge;std::array<int,N> at;int x,z;
  Geometry(int aa,int bb,const std::vector<Point>& d):a(aa),b(bb){
    for(auto&v:edge)v.fill(-1);std::array<int,N*N> index;index.fill(-1);std::vector<Point> rep{{0,0}};index[qkey(a,b,0,0)]=0;
    for(size_t i=0;i<rep.size();++i){auto [u,v]=rep[i];for(auto [dx,dy]:std::array<Point,2>{{{1,0},{0,1}}}){int k=qkey(a,b,u+dx,v+dy);if(index[k]<0){index[k]=rep.size();rep.push_back({u+dx,v+dy});}}}if(rep.size()!=N)throw std::logic_error("quotient");
    for(int i=0;i<N;++i){auto [u,v]=rep[i];nn[i]={index[qkey(a,b,u,v+1)],index[qkey(a,b,u+1,v)],index[qkey(a,b,u,v-1)],index[qkey(a,b,u-1,v)]};diag[i]={index[qkey(a,b,u+1,v+1)],index[qkey(a,b,u-1,v+1)],index[qkey(a,b,u-1,v-1)],index[qkey(a,b,u+1,v-1)]};}
    int ne=0;for(int i=0;i<N;++i)for(int j:{0,1}){int u=nn[i][j];edge[i][j]=edge[u][j+2]=ne++;}if(ne!=2*N)throw std::logic_error("edges");
    for(int i=0;i<N;++i)at[i]=index[qkey(a,b,d[i].first,d[i].second)];std::array<bool,N> seen{};for(int v:at){if(seen[v])throw std::logic_error("catalogue collision");seen[v]=true;}
    auto find=[&](Point p){auto it=std::find(d.begin(),d.end(),p);if(it==d.end())throw std::logic_error("missing direction");return int(it-d.begin());};x=at[find({0,0})];z=at[find({1,0})];
  }
  State eval(const std::array<unsigned char,N>&o)const{
    DSU bd,wd;int k=0,e=0,f=0;for(int v=0;v<N;++v)k+=o[v];for(int v=0;v<N;++v){if(o[v]){for(int j:{0,1})if(o[nn[v][j]]){++e;bd.join(v,nn[v][j]);}}else{for(int j:{0,1})if(!o[nn[v][j]])wd.join(v,nn[v][j]);for(int j:{0,1})if(!o[diag[v][j]])wd.join(v,diag[v][j]);}if(o[v]&&o[nn[v][1]]&&o[nn[v][0]]&&o[diag[v][0]])++f;}
    State s;s.black.fill(-1);std::array<bool,N> bs{},ws{};int bc=0,wc=0;for(int v=0;v<N;++v){if(o[v]){int r=bd.root(v);s.black[v]=r;if(!bs[r])bs[r]=1,++bc;}else{int r=wd.root(v);if(!ws[r])ws[r]=1,++wc;}}s.q=bc-wc-(k-e+f);if(std::abs(s.q)>1)throw std::logic_error("rank");return s;
  }
  std::uint32_t bell(int y,const State&s,const std::array<unsigned char,N>&o)const{
    std::array<int,3*N> label;label.fill(-1);int next=0;std::uint32_t key=0;for(int side=0;side<2;++side)for(int j=0;j<4;++j){int c=side?y:x,u=nn[c][j],id=o[u]?s.black[u]:N+edge[c][j];if(label[id]<0)label[id]=next++;key|=std::uint32_t(label[id])<<(3*(4*side+j));}return key;
  }
};
}

int main(int argc,char**argv){try{
  if(argc!=10)throw std::invalid_argument("KERNEL OUTPUT SAMPLES SHARD_INDEX SHARD_COUNT SEED PROPOSAL_P BATCHES TOKEN");auto K=read_kernel(argv[1]);std::string output=argv[2];std::uint64_t samples=std::stoull(argv[3]);int shard=std::stoi(argv[4]),shards=std::stoi(argv[5]);std::uint64_t seed=std::stoull(argv[6]);double p=std::stod(argv[7]);if(std::stoi(argv[8])!=BATCHES||std::string(argv[9])!="p537-N145-full-T")throw std::invalid_argument("contract");if(shard<0||shard>=shards||!samples||p<=.5||p>=.7||std::ifstream(output).good())throw std::invalid_argument("arguments");
  auto disp=catalogue();Geometry geo[2]={{12,1,disp},{9,8,disp}};auto index_of=[&](Point p){return int(std::find(disp.begin(),disp.end(),p)-disp.begin());};int x=index_of({0,0}),z=index_of({1,0});std::array<int,3> nn3{index_of({-1,0}),index_of({0,1}),index_of({0,-1})};std::vector<int> source;for(int d=0;d<N;++d)if(d!=x&&d!=z)source.push_back(d);std::vector<Sum> sums(BATCHES*2*N);std::uint64_t begin=samples*shard/shards,end=samples*(shard+1)/shards;
  for(std::uint64_t s=begin;s<end;++s){int batch=s%BATCHES;std::array<unsigned char,N> bits{};for(int d:source)bits[d]=unit(seed^(s*0xd6e8feb86659fd93ULL)^(std::uint64_t(d)*0xa0761d6478bd642fULL))<p;int k=0;for(int d:source)k+=bits[d];
    for(int g=0;g<2;++g){std::array<unsigned char,N> o{};for(int d=0;d<N;++d)o[geo[g].at[d]]=bits[d];State state[2];std::int64_t total[2]{},near[2]{};for(int state_z=0;state_z<2;++state_z){o[geo[g].z]=state_z;state[state_z]=geo[g].eval(o);for(int d:source)if(!bits[d]){auto value=K[geo[g].bell(geo[g].at[d],state[state_z],o)];total[state_z]+=value;if(std::find(nn3.begin(),nn3.end(),d)!=nn3.end())near[state_z]+=value;}}auto&v=sums[(batch*2+g)*N+k];++v.count;v.a0+=total[0];v.qa0+=state[0].q*total[0];v.ea0+=state[0].q*state[0].q*total[0];v.a1+=total[1];v.qa1+=state[1].q*total[1];v.ea1+=state[1].q*state[1].q*total[1];v.n0+=near[0];v.qn0+=state[0].q*near[0];v.en0+=state[0].q*state[0].q*near[0];v.n1+=near[1];v.qn1+=state[1].q*near[1];v.en1+=state[1].q*state[1].q*near[1];}
  }
  std::ofstream out(output);if(!out)throw std::runtime_error("output");out<<"# schema=matching-one/p537-full-t-n145/v1\n# samples="<<samples<<"\n# shard_index="<<shard<<"\n# shard_count="<<shards<<"\n# seed="<<seed<<"\n# proposal_p="<<std::setprecision(17)<<p<<"\n# begin="<<begin<<"\n# end="<<end<<"\n";out<<"batch\tgeometry\tk\tcount\ta0\tqa0\tea0\ta1\tqa1\tea1\tnn0\tqnn0\tenn0\tnn1\tqnn1\tenn1\n";for(int b=0;b<BATCHES;++b)for(int g=0;g<2;++g)for(int k=0;k<N;++k){auto&v=sums[(b*2+g)*N+k];if(v.count)out<<b<<'\t'<<(g?"tilted":"axis")<<'\t'<<k<<'\t'<<v.count<<'\t'<<v.a0<<'\t'<<v.qa0<<'\t'<<v.ea0<<'\t'<<v.a1<<'\t'<<v.qa1<<'\t'<<v.ea1<<'\t'<<v.n0<<'\t'<<v.qn0<<'\t'<<v.en0<<'\t'<<v.n1<<'\t'<<v.qn1<<'\t'<<v.en1<<'\n';}
  std::cout<<"samples="<<end-begin<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
