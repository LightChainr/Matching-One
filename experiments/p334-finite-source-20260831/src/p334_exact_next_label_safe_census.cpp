// Exact finite vacant-label census on the original checkpoints, never a tail.
#define P334_CONTACT_LIBRARY_ONLY
#include "p334_next_label_contact_coordinates.cpp"

namespace {
Vector essential_line(Checkpoint& checkpoint, int old_rank) {
    if (old_rank != 1) return {0,0};
    for (int v = 0; v < checkpoint.geometry.n; ++v) {
        if (!checkpoint.active[v]) continue;
        const auto mark = checkpoint.uf.component_mark(v);
        if (mark.rank == 1) return mark.line;
    }
    throw std::runtime_error("R1 checkpoint has no existing essential component");
}

int after_rank(Checkpoint& checkpoint, int vertex, int old_rank, Vector line) {
    if (old_rank == 2) return 2;
    std::map<int,Vector> anchors;
    Vector basis = line;
    int rank = old_rank;
    for (const int edge_index: checkpoint.geometry.primal_incident[vertex]) {
        const auto& edge = checkpoint.geometry.primal_edges[edge_index];
        const int other = edge.i == vertex ? edge.j : edge.i;
        if (!checkpoint.active[other]) continue;
        const auto found = checkpoint.uf.find(other);
        const Vector step = edge.i == vertex ? Vector{edge.dx,edge.dy} : Vector{-edge.dx,-edge.dy};
        const Vector alpha{step.x-found.dx,step.y-found.dy};
        const auto inserted = anchors.emplace(found.root,alpha);
        if (inserted.second) continue;
        const auto anchor = inserted.first->second;
        const auto w = checkpoint.geometry.quotient.winding(alpha.x-anchor.x,alpha.y-anchor.y);
        if (w.x == 0 && w.y == 0) continue;
        if (rank == 0) { rank=1; basis=w; }
        else if (static_cast<__int128>(basis.x)*w.y != static_cast<__int128>(basis.y)*w.x)
            return 2;
    }
    return rank;
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: safe_census N output_directory code_commit");
        const int n = std::stoi(argv[1]);
        if (n != 325 && n != 425) throw std::runtime_error("only originalN325/N425");
        const std::filesystem::path out(argv[2]);
        if (std::filesystem::exists(out)) throw std::runtime_error("refuse to overwrite census");
        const auto prefixes = original_prefixes(n);
        const auto geometry0 = make_geometry({n,n==325?57:132,0,1});
        const auto geometry1 = make_geometry({n,n==325?18:268,0,1});
        Checkpoint first(geometry0),second(geometry1);
        const std::uint64_t seed = n == 325 ? 20260831430325ULL : 20260831430425ULL;
        std::filesystem::create_directories(out);
        gzFile output = gzopen((out/("N"+std::to_string(n)+".csv.gz")).c_str(),"wb1");
        if (!output) throw std::runtime_error("cannot write census");
        gzprintf(output,"N,batch,counter,k0,d,first_oldrank,second_oldrank,first_safe_count,second_safe_count,joint_safe_count\n");
        std::vector<int> permutation;
        std::uint64_t positions = 0, safe_sum = 0, allsafe = 0, zerosafe = 0;
        const auto start = std::chrono::steady_clock::now();
        for (const auto& item: prefixes) {
            const auto counter = item.first;
            const auto prefix = item.second;
            const int batch=prefix[0],k0=prefix[1],r0=prefix[2],r1=prefix[3],d=n-k0;
            counter_permutation(n,seed,counter,permutation);
            first.load(permutation,k0); second.load(permutation,k0);
            const auto line0=essential_line(first,r0),line1=essential_line(second,r1);
            int safe0=0,safe1=0,joint=0;
            for (int k=k0;k<n;++k) {
                const int vertex=permutation[k];
                const bool a=after_rank(first,vertex,r0,line0)==r0;
                const bool b=after_rank(second,vertex,r1,line1)==r1;
                safe0+=a; safe1+=b; joint+=a&&b; ++positions;
            }
            gzprintf(output,"%d,%d,%llu,%d,%d,%d,%d,%d,%d,%d\n",n,batch,
                     static_cast<unsigned long long>(counter),k0,d,r0,r1,safe0,safe1,joint);
            safe_sum+=joint; allsafe+=joint==d; zerosafe+=joint==0;
        }
        if (gzclose(output)!=Z_OK) throw std::runtime_error("census gzip failed");
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream meta(out/"metadata.json");
        meta<<std::setprecision(15)<<"{\n\"N\":"<<n<<",\n\"code_commit\":\""<<argv[3]
            <<"\",\n\"prefix_source_commit\":\"9c495ab13e65f2bc93dc0849ee3b73f88724c4b1\",\n"
            <<"\"prefixes\":"<<prefixes.size()<<",\n\"vacant_positions_enumerated\":"<<positions
            <<",\n\"joint_safe_sum\":"<<safe_sum<<",\n\"all_joint_safe_prefixes\":"<<allsafe
            <<",\n\"zero_joint_safe_prefixes\":"<<zerosafe<<",\n\"elapsed_seconds\":"<<seconds
            <<",\n\"new_samples\":0,\n\"tail_replays\":0,\n\"DP_calls\":0,\n"
            <<"\"safe_definition\":\"Both orientations retain their own original checkpoint rank; R2 is always rank-preserving; pi_joint_safe=joint_safe_count/d exactly\"\n}\n";
        meta.close();
        std::cout<<"N"<<n<<" CENSUS COMPLETE "<<positions<<" vacant positions in "<<seconds<<"s\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
