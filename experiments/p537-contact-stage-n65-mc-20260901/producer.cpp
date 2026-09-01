// Frozen N65 MC for the P537 kernel-changing contact x birth-stage carrier.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
constexpr int N=65, BATCHES=100, KEY_SPACE=1<<24;
using Point=std::pair<int,int>;
int mod(int x){x%=N;return x<0?x+N:x;}
int qkey(int a,int b,int x,int y){return N*mod(a*x+b*y)+mod(-b*x+a*y);}
std::uint64_t mix(std::uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
double unit(std::uint64_t x){return double(mix(x)>>11)*0x1.0p-53;}

std::vector<std::int32_t> kernel(const std::string& path){
  std::ifstream f(path);if(!f)throw std::runtime_error("kernel open");
  std::vector<std::int32_t> out(KEY_SPACE);std::string s;int kc=-1,gc=-1;
  while(std::getline(f,s)){if(s.empty()||s[0]=='#')continue;std::stringstream z(s);std::vector<std::string> v;std::string x;while(std::getline(z,x,'\t'))v.push_back(x);
    if(kc<0){for(int i=0;i<int(v.size());++i){if(v[i]=="key"||v[i]=="packed_key")kc=i;if(v[i]=="g16")gc=i;}continue;}
    out.at(std::stoul(v.at(kc)))=std::stoi(v.at(gc));
  }if(kc<0||gc<0)throw std::runtime_error("kernel schema");return out;
}

struct DSU{std::array<int,N> p,s;DSU(){for(int i=0;i<N;++i)p[i]=i,s[i]=1;}int root(int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}void join(int x,int y){x=root(x);y=root(y);if(x==y)return;if(s[x]<s[y])std::swap(x,y);p[y]=x;s[x]+=s[y];}};
struct State{int q;std::array<int,N> black,white;};
struct Sum{std::uint64_t n=0;std::int64_t q0=0,e0=0,a0=0,qa0=0,ea0=0,q1=0,e1=0,a1=0,qa1=0,ea1=0;void add(const State&x,const State&z,int u,int v){++n;q0+=x.q;e0+=x.q*x.q;a0+=u;qa0+=x.q*u;ea0+=x.q*x.q*u;q1+=z.q;e1+=z.q*z.q;a1+=v;qa1+=z.q*v;ea1+=z.q*z.q*v;}};

std::vector<Point> catalogue(){
  std::vector<std::tuple<int,int,int>> candidates;for(int x=-24;x<=24;++x)for(int y=-24;y<=24;++y)candidates.emplace_back(x*x+y*y,x,y);std::sort(candidates.begin(),candidates.end());
  std::array<bool,N*N> one{},two{};std::vector<Point> out;
  for(auto [r,x,y]:candidates){int u=qkey(8,1,x,y),v=qkey(7,4,x,y);if(!one[u]&&!two[v]){one[u]=two[v]=true;out.push_back({x,y});if(out.size()==N)break;}}
  if(out.size()!=N)throw std::logic_error("common displacement catalogue");return out;
}

struct Geometry{
  int a,b;std::array<std::array<int,4>,N> nn,diag,edge;std::array<int,N> at;int x,z;std::array<int,4> card,corner;
  Geometry(int aa,int bb,const std::vector<Point>& d):a(aa),b(bb){
    for(auto&v:edge)v.fill(-1);std::array<int,N*N> index;index.fill(-1);std::vector<Point> rep{{0,0}};index[qkey(a,b,0,0)]=0;
    for(size_t i=0;i<rep.size();++i){auto [u,v]=rep[i];for(auto [dx,dy]:std::array<Point,2>{{{1,0},{0,1}}}){int k=qkey(a,b,u+dx,v+dy);if(index[k]<0){index[k]=rep.size();rep.push_back({u+dx,v+dy});}}}if(rep.size()!=N)throw std::logic_error("quotient");
    for(int i=0;i<N;++i){auto [u,v]=rep[i];nn[i]={index[qkey(a,b,u,v+1)],index[qkey(a,b,u+1,v)],index[qkey(a,b,u,v-1)],index[qkey(a,b,u-1,v)]};diag[i]={index[qkey(a,b,u+1,v+1)],index[qkey(a,b,u-1,v+1)],index[qkey(a,b,u-1,v-1)],index[qkey(a,b,u+1,v-1)]};}
    int ne=0;for(int i=0;i<N;++i)for(int j:{0,1}){int u=nn[i][j];edge[i][j]=edge[u][j+2]=ne++;}if(ne!=2*N)throw std::logic_error("edges");
    for(int i=0;i<N;++i){at[i]=index[qkey(a,b,d[i].first,d[i].second)];}std::array<bool,N> seen{};for(int v:at){if(seen[v])throw std::logic_error("catalogue collision");seen[v]=true;}
    auto find=[&](Point p){return int(std::find(d.begin(),d.end(),p)-d.begin());};x=at[find({0,0})];z=at[find({1,0})];card=nn[z];corner={diag[z][0],diag[z][3],diag[z][2],diag[z][1]};std::array<bool,N> c{};c[z]=1;for(int v:card)c[v]=1;for(int v:corner)c[v]=1;if(std::count(c.begin(),c.end(),true)!=9)throw std::logic_error("collar not injective");
  }
  State eval(const std::array<unsigned char,N>&o)const{
    DSU bds,wds;int k=0,e=0,f=0;for(int v=0;v<N;++v)k+=o[v];for(int v=0;v<N;++v){if(o[v]){for(int j:{0,1})if(o[nn[v][j]]){++e;bds.join(v,nn[v][j]);}}else{for(int j:{0,1})if(!o[nn[v][j]])wds.join(v,nn[v][j]);for(int j:{0,1})if(!o[diag[v][j]])wds.join(v,diag[v][j]);}if(o[v]&&o[nn[v][1]]&&o[nn[v][0]]&&o[diag[v][0]])++f;}
    State s;s.black.fill(-1);s.white.fill(-1);std::array<bool,N> bs{},ws{};int bc=0,wc=0;for(int v=0;v<N;++v){if(o[v]){int r=bds.root(v);s.black[v]=r;if(!bs[r])bs[r]=1,++bc;}else{int r=wds.root(v);s.white[v]=r;if(!ws[r])ws[r]=1,++wc;}}s.q=bc-wc-(k-e+f);if(std::abs(s.q)>1)throw std::logic_error("rank");return s;
  }
  int outside(int c,int j,const State&s,const std::array<unsigned char,N>&o)const{int u=nn[c][j];return o[u]?s.black[u]:N+edge[c][j];}
  std::uint32_t bell(int y,const State&s,const std::array<unsigned char,N>&o)const{std::array<int,3*N> label;label.fill(-1);int next=0;std::uint32_t key=0;for(int side=0;side<2;++side)for(int j=0;j<4;++j){int id=outside(side?y:x,j,s,o);if(label[id]<0)label[id]=next++;key|=std::uint32_t(label[id])<<(3*(4*side+j));}return key;}
  int arm(int v,const std::array<unsigned char,N>&o)const{if(!o[v])return 0;if(v==card[0]||v==corner[0]||v==corner[3])return 1;if(v==card[2]||v==corner[1]||v==corner[2])return 2;return 0;}
  int contact(int y,const std::array<unsigned char,N>&o)const{int m=0;for(int c:{x,y})for(int j=0;j<4;++j)m|=arm(nn[c][j],o);return m;}
};

std::uint64_t ckey(int b,int g,int d,int k,int stage,int mask){return (((((std::uint64_t(b)*2+g)*N+d)*N+k)*2+stage)*4+mask);}
void emit(std::ofstream&o,int b,const char*kind,const char*g,Point d,const char*stage,int mask,int k,const Sum&s){o<<b<<'\t'<<kind<<'\t'<<g<<'\t'<<d.first<<'\t'<<d.second<<'\t'<<stage<<'\t'<<mask<<'\t'<<k<<'\t'<<s.n<<'\t'<<s.q0<<'\t'<<s.e0<<'\t'<<s.a0<<'\t'<<s.qa0<<'\t'<<s.ea0<<'\t'<<s.q1<<'\t'<<s.e1<<'\t'<<s.a1<<'\t'<<s.qa1<<'\t'<<s.ea1<<'\n';}
}

int main(int argc,char**argv){try{
  if(argc!=10)throw std::invalid_argument("KERNEL OUTPUT SAMPLES SHARD_INDEX SHARD_COUNT SEED PROPOSAL_P BATCHES TOKEN");auto K=kernel(argv[1]);std::string output=argv[2];std::uint64_t samples=std::stoull(argv[3]);int shard=std::stoi(argv[4]),shards=std::stoi(argv[5]);std::uint64_t seed=std::stoull(argv[6]);double p=std::stod(argv[7]);if(std::stoi(argv[8])!=BATCHES||std::string(argv[9])!="frozen-N65-contact-stage")throw std::invalid_argument("frozen contract");if(shard<0||shard>=shards||!samples||p<=.5||p>=.7||std::ifstream(output).good())throw std::invalid_argument("arguments");
  auto disp=catalogue();Geometry geo[2]={{8,1,disp},{7,4,disp}};std::vector<int> source;for(int d=0;d<N;++d)if(disp[d]!=Point{0,0}&&disp[d]!=Point{1,0})source.push_back(d);std::vector<Sum> global(BATCHES*2*N*N);std::unordered_map<std::uint64_t,Sum> carrier;std::uint64_t begin=samples*shard/shards,end=samples*(shard+1)/shards;
  for(std::uint64_t s=begin;s<end;++s){int batch=s%BATCHES;std::array<unsigned char,N> bits{};for(int d:source)bits[d]=unit(seed^(s*0xd6e8feb86659fd93ULL)^(std::uint64_t(d)*0xa0761d6478bd642fULL))<p;int k=0;for(int d:source)k+=bits[d];
    for(int g=0;g<2;++g){std::array<unsigned char,N> o{};for(int d=0;d<N;++d)o[geo[g].at[d]]=bits[d];State s0=geo[g].eval(o);int arms=0;for(int j=0;j<4;++j)arms|=int(o[geo[g].card[j]])<<j;o[geo[g].z]=1;State s1=geo[g].eval(o);o[geo[g].z]=0;int stage=s0.q==-1&&s1.q==0?0:(s0.q==0&&s1.q==1?1:-1);
      for(int d:source){int y=geo[g].at[d],a0=0,a1=0;std::uint32_t b0=KEY_SPACE,b1=KEY_SPACE;if(!o[y]){b0=geo[g].bell(y,s0,o);o[geo[g].z]=1;b1=geo[g].bell(y,s1,o);o[geo[g].z]=0;a0=K[b0];a1=K[b1];}global[(((batch*2+g)*N+d)*N+k)].add(s0,s1,a0,a1);int m=!o[y]?geo[g].contact(y,o):0;if((arms==5||arms==10)&&stage>=0&&!o[y]&&b0!=b1&&a0!=a1&&m>=1&&m<=3)carrier[ckey(batch,g,d,k,stage,m)].add(s0,s1,a0,a1);}
    }}
  std::ofstream out(output);if(!out)throw std::runtime_error("output");out<<"# schema=matching-one/p537-contact-stage-n65-mc/v1\n# samples="<<samples<<"\n# shard_index="<<shard<<"\n# shard_count="<<shards<<"\n# seed="<<seed<<"\n# proposal_p="<<std::setprecision(17)<<p<<"\n# begin="<<begin<<"\n# end="<<end<<"\n";out<<"batch\tkind\tgeometry\tdx\tdy\tstage\tcontact_mask\tk\tcount\tsum_q0\tsum_E0\tsum_a16_0\tsum_q0_a16_0\tsum_E0_a16_0\tsum_q1\tsum_E1\tsum_a16_1\tsum_q1_a16_1\tsum_E1_a16_1\n";for(int b=0;b<BATCHES;++b)for(int g=0;g<2;++g)for(int d:source)for(int k=0;k<N;++k){auto&s=global[(((b*2+g)*N+d)*N+k)];if(s.n)emit(out,b,"global",g?"tilted":"axis",disp[d],"-",0,k,s);for(int st=0;st<2;++st)for(int m=1;m<=3;++m){auto it=carrier.find(ckey(b,g,d,k,st,m));if(it!=carrier.end())emit(out,b,"carrier",g?"tilted":"axis",disp[d],st?"12":"01",m,k,it->second);}}std::cout<<"samples="<<end-begin<<" carrier_cells="<<carrier.size()<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
