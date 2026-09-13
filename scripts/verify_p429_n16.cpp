// Independent N16 verification: weighted union-find on row-major square tori,
// direct subset survival, and exhaustive selected-prefix permutations.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using U=std::uint32_t;
struct PotentialUF {
 std::array<int,16> parent{},sz{},dx{},dy{};
 int firstx=0,firsty=0,rank=0;
 PotentialUF(){std::iota(parent.begin(),parent.end(),0);sz.fill(1);}
 struct Root{int r,x,y;};
 Root find(int v){if(parent[v]==v)return{v,0,0}; auto z=find(parent[v]);dx[v]+=z.x;dy[v]+=z.y;parent[v]=z.r;return{z.r,dx[v],dy[v]};}
 void edge(int i,int j,int ex,int ey){
  auto a=find(i),b=find(j); int x=a.x+ex-b.x,y=a.y+ey-b.y;
  if(a.r==b.r){
   if(x%4||y%4)throw std::runtime_error("invalid closed displacement");
   if(x||y){if(!rank){firstx=x;firsty=y;rank=1;}else if(firstx*y-firsty*x)rank=2;}
   return;
  }
  if(sz[a.r]<sz[b.r]){std::swap(a,b);x=-x;y=-y;}
  parent[b.r]=a.r;dx[b.r]=x;dy[b.r]=y;sz[a.r]+=sz[b.r];
 }
};
int main(int argc,char**argv){
 std::array<int,65536> rr{},lx{},ly{};
 for(U s=0;s<65536;++s){
  PotentialUF f;
  for(int y=0;y<4;++y)for(int x=0;x<4;++x){int u=x+4*y; if(!(s&(1u<<u)))continue;
   int v=(x+1)%4+4*y; if(s&(1u<<v))f.edge(u,v,1,0);
   v=x+4*((y+1)%4);if(s&(1u<<v))f.edge(u,v,0,1);
  }
  rr[s]=f.rank;
  if(f.rank==1){int x=f.firstx/4,y=f.firsty/4,g=std::gcd(std::abs(x),std::abs(y));x/=g;y/=g;if(x<0||(x==0&&y<0)){x=-x;y=-y;}lx[s]=x;ly[s]=y;}
 }
 if(argc==2){std::ofstream o(argv[1]);for(U s=0;s<65536;++s)o<<s<<' '<<rr[s]<<' '<<lx[s]<<' '<<ly[s]<<'\n';}
 std::array<int,9> target{{1,7,18,20,8,0,0,0,0}};
 std::vector<U> stratum;
 std::array<std::uint64_t,9> prefixes{},next_three{},fork_success{};
 auto survival=[&](U s){std::array<int,9> c{};U rest=65535^s;for(U a=rest;;a=(a-1)&rest){if(rr[s|a]==1)++c[__builtin_popcount(a)];if(!a)break;}return c;};
 auto vacancy_statistics=[&](U s){std::array<int,2> z{};U rem=65535^s;
  for(int v=0;v<16;++v)if(rem&(1u<<v)){U t=s|(1u<<v);if(rr[t]!=1)continue;int safe=0;
   for(int w=0;w<16;++w)if(!(t&(1u<<w))&&rr[t|(1u<<w)]==1)++safe;
   if(safe==4){++z[0];}
   z[1]+=safe*safe;
  }return z;};
 for(U s=0;s<65536;++s)if(__builtin_popcount(s)==8 && rr[s]==1 && lx[s]==1 && ly[s]==0 && survival(s)==target){
  stratum.push_back(s);std::vector<int> order;
  for(int v=0;v<16;++v)if(s&(1u<<v))order.push_back(v);
  auto st=vacancy_statistics(s);
  do{U m=0;int birth=0;for(int j=0;j<8;++j){m|=1u<<order[j];if(!birth && rr[m]>=1)birth=j+1;}
   if(!birth)throw std::runtime_error("rank-one final state has no birth");
   ++prefixes[birth];next_three[birth]+=st[0];fork_success[birth]+=st[1];
  }while(std::next_permutation(order.begin(),order.end()));
 }
 std::cout<<"{\"rank_one_states\":"<<std::count(rr.begin(),rr.end(),1)<<",\"stratum_states\":"<<stratum.size()<<",\"witnesses\":[";
 bool comma=false; for(U s: {12463u,4343u}){if(comma)std::cout<<',';comma=true;auto st=vacancy_statistics(s);auto surv=survival(s);
 std::cout<<"{\"rowmajor_mask\":"<<s<<",\"signature\":[";for(int i=0;i<9;++i){if(i)std::cout<<',';std::cout<<surv[i];}std::cout<<"],\"next_three_choices\":"<<st[0]<<",\"fork_success_of_392\":"<<st[1]<<'}';}
 std::cout<<"],\"cohorts\":[";comma=false; for(int j=1;j<=8;++j)if(prefixes[j]){if(comma)std::cout<<',';comma=true;std::cout<<"{\"K1\":"<<j<<",\"prefixes\":"<<prefixes[j]<<",\"next_three_weight\":"<<next_three[j]<<",\"fork_success_weight\":"<<fork_success[j]<<'}';}
 std::cout<<"]}\n";
}
