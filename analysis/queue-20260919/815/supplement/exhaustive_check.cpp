// Issue 815 only. Integer lift-DFS adaptation of PR773's torus_rank,
// not the missing rank775 C kernel. White paths use a planar face-hub graph.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct Data { int comp[25], x[25], y[25]; };
struct Grid {
    int L,N,nb[25][4];
    const int dx[4]={1,-1,0,0}, dy[4]={0,0,1,-1};
    explicit Grid(int l):L(l),N(l*l) {
        for(int a=0;a<N;a++) for(int d=0;d<4;d++)
            nb[a][d]=index(a%L+dx[d],a/L+dy[d]);
    }
    int index(int x,int y) const { return ((y%L+L)%L)*L+(x%L+L)%L; }
    bool black(uint32_t mask,int x,int y) const {return (mask>>index(x,y))&1u;}
    int rank(uint32_t mask,Data& s) const {
        fill(s.comp,s.comp+N,-1); int bx=0,by=0,r=0;
        for(int root=0;root<N;root++) {
            if(!((mask>>root)&1u)||s.comp[root]>=0) continue;
            int stack[25],n=0; stack[n++]=root;s.comp[root]=root;s.x[root]=s.y[root]=0;
            while(n) {
                int a=stack[--n];
                for(int d=0;d<4;d++) {
                    int b=nb[a][d]; if(!((mask>>b)&1u)) continue;
                    int x=s.x[a]+dx[d], y=s.y[a]+dy[d];
                    if(s.comp[b]<0) {
                        s.comp[b]=root;s.x[b]=x;s.y[b]=y;stack[n++]=b;
                    } else {
                        int vx=x-s.x[b],vy=y-s.y[b];
                        if(!vx&&!vy) continue;
                        if(vx%L||vy%L) throw runtime_error("nonintegral deck offset");
                        vx/=L;vy/=L;
                        if(!r){bx=vx;by=vy;r=1;}
                        else if(bx*vy-by*vx) return 2;
                    }
                }
            }
        }
        return r;
    }
};
struct Flow {
    struct E {int v,rev,cap;}; vector<vector<E>> a;
    explicit Flow(int n):a(n){}
    void edge(int u,int v,int cap) {
        E x{v,(int)a[v].size(),cap},y{u,(int)a[u].size(),0};a[u].push_back(x);a[v].push_back(y);
    }
    int solve(int source,int sink) {
        int ans=0,n=a.size();
        while(true) {
            vector<int> pv(n,-1),pe(n,-1);queue<int> q;q.push(source);pv[source]=source;
            while(!q.empty()&&pv[sink]<0) {
                int u=q.front();q.pop();
                for(int j=0;j<(int)a[u].size();j++) if(a[u][j].cap&&pv[a[u][j].v]<0) {
                    int v=a[u][j].v;pv[v]=u;pe[v]=j;q.push(v);
                }
            }
            if(pv[sink]<0) break;
            for(int v=sink;v!=source;v=pv[v]) {
                E &e=a[pv[v]][pe[v]];--e.cap;++a[v][e.rev].cap;
            }
            ++ans;
        }
        return ans;
    }
};
const int ringx[8]={1,1,0,-1,-1,-1,0,1};
const int ringy[8]={0,1,1,1,0,-1,-1,-1};
const int dirring[4]={0,4,2,6};

bool arm_check(const Grid& g,uint32_t mask,const vector<int>& dirs,int R) {
    auto inside=[&](int x,int y){int r=max(abs(x),abs(y));return r>=1&&r<=R;};
    auto ix=[&](int x,int y){return (y+R)*(2*R+1)+x+R;};
    int width=2*R+1;
    // Each selected black lift component must reach the outer lattice boundary.
    // Distinctness of the lift components is checked by the caller's integer offsets.
    for(int d:dirs) {
        vector<int> seen(width*width,0);queue<pair<int,int>> q;
        q.push({g.dx[d],g.dy[d]});seen[ix(g.dx[d],g.dy[d])]=1;bool hit=false;
        while(!q.empty()) {
            auto [x,y]=q.front();q.pop();
            if(max(abs(x),abs(y))==R){hit=true;break;}
            for(int j=0;j<4;j++) {
                int xx=x+g.dx[j],yy=y+g.dy[j];
                if(inside(xx,yy)&&g.black(mask,xx,yy)&&!seen[ix(xx,yy)]) {
                    seen[ix(xx,yy)]=1;q.push({xx,yy});
                }
            }
        }
        if(!hit) return false;
    }
    // White lattice vertices have capacity 1. Face hubs also have capacity 1,
    // prohibiting two crossing diagonal paths through one face.
    int ids[25];fill(ids,ids+25,-1);int vertices=0;
    for(int y=-R;y<=R;y++) for(int x=-R;x<=R;x++)
        if(inside(x,y)&&!g.black(mask,x,y))ids[ix(x,y)]=vertices++;
    struct Face {int id;vector<int> corners;};vector<Face> faces;
    for(int y=-R;y<R;y++) for(int x=-R;x<R;x++) {
        bool ok=true;vector<int> white;
        for(int b=0;b<=1;b++)for(int a=0;a<=1;a++) {
            if(!inside(x+a,y+b)){ok=false;continue;}
            int z=ids[ix(x+a,y+b)];if(z>=0)white.push_back(z);
        }
        if(ok&&white.size()>=2)faces.push_back({vertices++,white});
    }
    int k=dirs.size(),source=2*vertices+k,sink=source+1;Flow f(sink+1);
    for(int z=0;z<vertices;z++) f.edge(2*z,2*z+1,1);
    for(int y=-R;y<=R;y++)for(int x=-R;x<=R;x++)if(inside(x,y)) {
        int u=ids[ix(x,y)];if(u<0)continue;
        if(max(abs(x),abs(y))==R)f.edge(2*u+1,sink,1);
        for(int d=0;d<4;d++) {
            int xx=x+g.dx[d],yy=y+g.dy[d];
            if(inside(xx,yy)) {int v=ids[ix(xx,yy)];if(v>=0)f.edge(2*u+1,2*v,k);}
        }
    }
    for(const auto& face:faces)for(int v:face.corners) {
        f.edge(2*v+1,2*face.id,k);f.edge(2*face.id+1,2*v,k);
    }
    vector<int> starts;for(int d:dirs)starts.push_back(dirring[d]);sort(starts.begin(),starts.end());
    for(int j=0;j<k;j++) {
        int from=starts[j],to=starts[(j+1)%k],group=2*vertices+j;f.edge(source,group,1);
        for(int h=(from+1)%8;h!=to;h=(h+1)%8) {
            int z=ids[ix(ringx[h],ringy[h])];if(z>=0)f.edge(group,2*z,1);
        }
    }
    return f.solve(source,sink)==k;
}

int main(int argc,char** argv) {
    try {
        if(argc!=3)throw runtime_error("usage: exhaustive_check --ranks|--census L (1<=L<=5)");
        int L=stoi(argv[2]);if(L<1||L>5)throw runtime_error("L outside authorized range");
        Grid g(L);Data before,after;string mode=argv[1];
        if(mode=="--ranks") {
            uint32_t mask;while(cin>>mask) {
                if(mask>=(1u<<g.N))throw runtime_error("mask out of range");
                cout<<g.rank(mask,before)<<'\n';
            }
            return 0;
        }
        if(mode!="--census")throw runtime_error("unknown mode");
        auto begin=chrono::steady_clock::now();uint64_t total=1ull<<(g.N-1),zero=0,jump=0,t3=0,t4=0,rs=0,bad=0;
        int R=max(1,(L-1)/2);vector<uint32_t> failures;
        for(uint64_t rest=0;rest<total;rest++) {
            uint32_t mask=rest<<1;
            if(g.rank(mask,before)!=0)continue;
            ++zero;if(g.rank(mask|1u,after)!=2)continue;++jump;
            // L=1 has periodic self-loops and no outside-site arm event.
            if(L==1)continue;
            vector<int> ds;int cx[4],cy[4],cc[4];
            for(int d=0;d<4;d++) {
                int u=g.nb[0][d];if(!((mask>>u)&1u))continue;
                ds.push_back(d);cc[d]=before.comp[u];cx[d]=g.dx[d]-before.x[u];cy[d]=g.dy[d]-before.y[u];
            }
            bool isT=false,pass=false;int contacts=0;
            for(int a=0;a<(int)ds.size();a++)for(int b=a+1;b<(int)ds.size();b++)for(int c=b+1;c<(int)ds.size();c++) {
                int i=ds[a],j=ds[b],k=ds[c];
                if(cc[i]!=cc[j]||cc[i]!=cc[k])continue;
                if((cx[j]-cx[i])*(cy[k]-cy[i])-(cy[j]-cy[i])*(cx[k]-cx[i])==0)continue;
                isT=true;contacts=0;for(int d:ds)if(cc[d]==cc[i])++contacts;
                if(arm_check(g,mask,{i,j,k},R))pass=true;
            }
            if(isT){if(contacts==3)++t3;else ++t4;}
            else {
                bool form=false;
                if(ds.size()==4)for(int a=1;a<4;a++) {
                    int i=ds[0],j=ds[a];if(cc[i]!=cc[j])continue;
                    vector<int> other;for(int b=1;b<4;b++)if(b!=a)other.push_back(ds[b]);
                    int k=other[0],l=other[1];if(cc[k]!=cc[l]||cc[i]==cc[k])continue;
                    if((cx[j]-cx[i])*(cy[l]-cy[k])-(cy[j]-cy[i])*(cx[l]-cx[k]))form=true;
                }
                if(form){++rs;pass=arm_check(g,mask,ds,R);}
            }
            if(!pass){++bad;if(failures.size()<8)failures.push_back(mask);}
        }
        double seconds=chrono::duration<double>(chrono::steady_clock::now()-begin).count();
        cout<<"{\"L\":"<<L<<",\"rest_masks\":"<<total<<",\"rank_calls\":"<<total+zero
            <<",\"closed_rank_zero\":"<<zero<<",\"jump2\":"<<jump<<",\"T3\":"<<t3
            <<",\"T4plus\":"<<t4<<",\"Rsplit\":"<<rs<<",\"arm_failures\":"<<bad
            <<",\"inner_radius\":1,\"outer_radius\":"<<R<<",\"periodic_self_loop_exclusion\":"<<(L==1?"true":"false")
            <<",\"injective_box\":"<<(L>=3?"true":"false")<<",\"wall_seconds\":"<<seconds<<",\"first_failure_masks\":[";
        for(size_t j=0;j<failures.size();j++){if(j)cout<<',';cout<<failures[j];}
        cout<<"]}\n";return bad?1:0;
    } catch(const exception& e){cerr<<e.what()<<'\n';return 2;}
}
