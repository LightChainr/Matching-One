// Independent re-derivation of the ambient H1 rank of the occupied site subgraph on an LxL torus.
// Consumer: docs/research-status/control-20260919-rev2.md section 1.
// Usage: ./rank_census <lattice 0=square-NN 1=triangular> <L>
// lattice 0 = square NN (G4), 1 = triangular (6 nbrs: (+-1,0),(0,+-1),(1,-1),(-1,1))
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int L, N, LT;
static int nnb;
static int off[6][2];

static int di_[6], dj_[6];

// returns ambient rank 0/1/2 of configuration bits `occ`
static inline int rank_of(unsigned int occ){
    static int dispI[64], dispJ[64], seen[64], stack[64];
    int vx=0, vy=0, have1=0, rank=0;
    for(int s=0;s<N;s++) seen[s]=0;
    for(int s0=0;s0<N;s0++){
        if(!((occ>>s0)&1u) || seen[s0]) continue;
        int sp=0; stack[sp++]=s0; seen[s0]=1; dispI[s0]=0; dispJ[s0]=0;
        while(sp){
            int cur=stack[--sp];
            int ci=cur/L, cj=cur%L;
            for(int k=0;k<nnb;k++){
                int ri=ci+off[k][0], rj=cj+off[k][1];
                int wi=0,wj=0;
                int ni=ri, nj=rj;
                if(ni<0){ni+=L; wi=-1;} else if(ni>=L){ni-=L; wi=1;}
                if(nj<0){nj+=L; wj=-1;} else if(nj>=L){nj-=L; wj=1;}
                int nb=ni*L+nj;
                if(!((occ>>nb)&1u)) continue;
                int pi=dispI[cur]+wi, pj=dispJ[cur]+wj;
                if(!seen[nb]){ seen[nb]=1; dispI[nb]=pi; dispJ[nb]=pj; stack[sp++]=nb; }
                else {
                    int cx=pi-dispI[nb], cy=pj-dispJ[nb];
                    if(cx||cy){
                        if(rank==2) continue;
                        if(!have1){ vx=cx; vy=cy; have1=1; rank=1; }
                        else if(vx*cy-vy*cx != 0){ rank=2; }
                    }
                }
            }
        }
        if(rank==2) { /* keep scanning is unnecessary for rank but cheap to stop */ break; }
    }
    return rank;
}

int main(int argc,char**argv){
    int lat=atoi(argv[1]); L=atoi(argv[2]); N=L*L;
    if(lat==0){ nnb=4; int o[4][2]={{1,0},{-1,0},{0,1},{0,-1}}; memcpy(off,o,sizeof(o)); }
    else { nnb=6; int o[6][2]={{1,0},{-1,0},{0,1},{0,-1},{1,-1},{-1,1}}; memcpy(off,o,sizeof(o)); }

    long long rk[3]={0,0,0};
    unsigned long long total=1ULL<<N;
    for(unsigned long long c=0;c<total;c++) rk[rank_of((unsigned int)c)]++;
    printf("lattice=%d L=%d  full rank counts over 2^%d: r0=%lld r1=%lld r2=%lld\n",lat,L,N,rk[0],rk[1],rk[2]);

    // fixed site v = index 0; enumerate external configs on the other N-1 sites
    static long long H[3][40]; memset(H,0,sizeof(H)); long long d[3]={0,0,0}; long long bad=0;
    unsigned long long ext=1ULL<<(N-1);
    for(unsigned long long e=0;e<ext;e++){
        unsigned int base=(unsigned int)(e<<1); // bit0 reserved for v
        int r0=rank_of(base);
        int r1=rank_of(base|1u);
        int dd=r1-r0;
        if(dd<0||dd>2){bad++; continue;}
        d[dd]++; H[dd][__builtin_popcount(base)]++;
    }
    printf("           fixed-site delta counts over 2^%d: d0=%lld d1=%lld d2=%lld  (out-of-range=%lld)\n",N-1,d[0],d[1],d[2],bad);
    for(int dd=0;dd<3;dd++){ printf("H%d:",dd); for(int k=0;k<N;k++) printf(" %lld",H[dd][k]); printf("\n"); }
    return 0;
}
