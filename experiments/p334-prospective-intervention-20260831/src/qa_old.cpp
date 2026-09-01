// Reconstruct only named old archived prefixes; no fresh counter and no tail.
#define P334_PROSPECTIVE_LIBRARY_ONLY
#include "prospective.cpp"
int main(int argc,char**argv) {
    if(argc!=3)return 2;
    int n=std::stoi(argv[1]);if(n!=325&&n!=425)return 2;
    std::ifstream input(argv[2]);std::string line;std::getline(input,line);
    auto g0=make_geometry({n,n==325?57:132,0,1}),g1=make_geometry({n,n==325?18:268,0,1});
    Checkpoint a(g0),b(g1);std::vector<int> perm,vacant;
    std::uint64_t seed=n==325?20260831430325ULL:20260831430425ULL;
    std::cout<<"index,counter,first_mass,first_energy,first_degree,first_loop,second_mass,second_energy,second_degree,second_loop,W_first,W_second\n";
    int index=-1;
    while(std::getline(input,line)&&++index<1024) {
        auto x=parse_csv(line);if(x.size()!=10||x[0]!=static_cast<unsigned>(n))return 3;
        if(x[2]!=(n==325?43032500000ULL:43042500000ULL)+index)return 4;
        counter_permutation(n,seed,x[2],perm);int k0=x[3];a.load(perm,k0);b.load(perm,k0);
        if(prospective::checkpoint_rank(a)!=static_cast<int>(x[6])||prospective::checkpoint_rank(b)!=static_cast<int>(x[9]))return 5;
        if(x[6]||x[9])continue;
        prospective::Row row{};vacant.assign(perm.begin()+k0,perm.end());std::sort(vacant.begin(),vacant.end());
        prospective::census(a,b,vacant,row);
        std::cout<<index<<','<<x[2]<<std::setprecision(17);
        for(int o=0;o<2;++o)for(int j=0;j<4;++j)std::cout<<','<<row.feature[o][j];
        for(int o=0;o<2;++o)std::cout<<','<<row.weight[o];
        std::cout<<'\n';
    }
    return 0;
}
