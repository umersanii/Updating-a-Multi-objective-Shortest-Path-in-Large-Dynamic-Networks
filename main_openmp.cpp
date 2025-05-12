#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_map>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

const double INF = numeric_limits<double>::infinity();

struct Edge {
    int u, v;
    vector<double> weights;
    Edge(int _u, int _v, const vector<double>& _w) : u(_u), v(_v), weights(_w) {}
};

struct CombinedEdge {
    int u, v;
    double weight;
    vector<double> original_weights;
    CombinedEdge(int _u, int _v, double _w, const vector<double>& _ow)
      : u(_u), v(_v), weight(_w), original_weights(_ow) {}
};

void dijkstra(
    const vector<vector<pair<int, double>>>& graph,
    int source,
    vector<int>& parent,
    vector<double>& dist,
    int N)
{
    dist.assign(N, numeric_limits<double>::infinity());
    parent.assign(N, -1);
    dist[source] = 0.0;

    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    pq.emplace(0.0, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        for (const auto& [v, w] : graph[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.emplace(dist[v], v);
            }
        }
        // cout << ".";
    }
}

void sosp_update(
    vector<vector<pair<int,double>>>& graph,
    int source,
    vector<int>& parent,
    vector<double>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices,
    int obj_idx)
{
    const double INF = 1e9;

    // Reverse graph construction
    vector<vector<pair<int, double>>> rev(num_vertices);
    for (int u = 0; u < num_vertices; u++) {
        for (int i = 0; i < graph[u].size(); i++) {
            int v = graph[u][i].first;
            double w = graph[u][i].second;
            if (w == INF) continue;
            rev[v].push_back(make_pair(u, w));
        }
    }

    // byV[v] = list of incoming edges to v from inserted_edges
    vector<vector<pair<int, double>>> byV(num_vertices);
    for (int i = 0; i < inserted_edges.size(); i++) {
        const Edge& e = inserted_edges[i];
        if (e.weights[obj_idx] != INF) {
            byV[e.v].push_back(make_pair(e.u, e.weights[obj_idx]));
        }
    }

    vector<int> affected;
    vector<char> marked(num_vertices, 0);

    if (inserted_edges.empty()) {
        distances.assign(num_vertices, INF);
        parent.assign(num_vertices, -1);
        distances[source] = 0;
        marked[source] = 1;
        affected.push_back(source);

        for (int i = 0; i < graph[source].size(); i++) {
            int v = graph[source][i].first;
            double w = graph[source][i].second;
            if (w != INF && !marked[v]) {
                marked[v] = 1;
                affected.push_back(v);
            }
        }
    } else {
        for (int v = 0; v < num_vertices; v++) {
            for (int i = 0; i < byV[v].size(); i++) {
                int u = byV[v][i].first;
                double w = byV[v][i].second;
                if (w == INF) continue;
                double nd = distances[u] + w;
                if (nd < distances[v]) {
                    distances[v] = nd;
                    parent[v] = u;
                    if (!marked[v]) {
                        marked[v] = 1;
                        affected.push_back(v);
                    }
                }
            }
        }
    }

    while (!affected.empty()) {
        vector<char> inN(num_vertices, 0);

        for (int i = 0; i < affected.size(); i++) {
            int v = affected[i];
            for (int j = 0; j < graph[v].size(); j++) {
                int nbr = graph[v][j].first;
                inN[nbr] = 1;
            }
        }

        vector<int> N;
        for (int v = 0; v < num_vertices; v++) {
            if (inN[v]) {
                N.push_back(v);
            }
        }

        vector<int> new_aff;
        for (int i = 0; i < N.size(); i++) {
            int v = N[i];

            // Backward edges (rev)
            for (int j = 0; j < rev[v].size(); j++) {
                int u = rev[v][j].first;
                double w = rev[v][j].second;
                if (w == INF) continue;
                if (marked[u]) {
                    double nd = distances[u] + w;
                    if (nd < distances[v]) {
                        distances[v] = nd;
                        parent[v] = u;
                        if (!marked[v]) {
                            marked[v] = 1;
                            new_aff.push_back(v);
                        }
                    }
                }
            }

            // Forward edges
            for (int j = 0; j < graph[v].size(); j++) {
                int w_node = graph[v][j].first;
                double weight = graph[v][j].second;
                if (weight == INF) continue;
                double nd = distances[v] + weight;
                if (nd < distances[w_node]) {
                    distances[w_node] = nd;
                    parent[w_node] = v;
                    if (!marked[w_node]) {
                        marked[w_node] = 1;
                        new_aff.push_back(w_node);
                    }
                }
            }
        }

        affected.swap(new_aff);
    }
}


struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

void create_combined_graph(const vector<vector<int>>& parents, int num_vertices, int num_obj, const vector<Edge>& orig, const vector<Edge>& ins,vector<vector<pair<int, double>>>& comb, vector<vector<CombinedEdge>>& combE)
{
    comb.assign(num_vertices, {});
    combE.assign(num_vertices, {});

    //////////////////////////////////////////////////////////////////debug
    // auto start_step1 = high_resolution_clock::now();
    ////////////////////////////////////////////////////////////////////////
    // Step 1: Build set S in parallel (collect and merge later)
    vector<set<pair<int, int>>> local_Sets(omp_get_max_threads());

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < num_obj; i++) {
        for (int v = 0; v < num_vertices; v++) {
            if (parents[i][v] >= 0) {
                int tid = omp_get_thread_num();
                local_Sets[tid].insert({parents[i][v], v});
            }
        }
    }

    // Merge local sets
    set<pair<int, int>> S;
    for (auto& ls : local_Sets)
        S.insert(ls.begin(), ls.end());

    ////////////////////////////////////////////////////////// debug
    // auto end_step1 = high_resolution_clock::now();
    // cout << "Step 1 completed in " 
    //      << duration_cast<milliseconds>(end_step1 - start_step1).count() << " ms." << endl;

    // auto start_step2 = high_resolution_clock::now();
    //////////////////////////////////////////////////////////////////////////
    // Step 2: Build edge_weights
    unordered_map<pair<int, int>, vector<double>, pair_hash> edge_weights;

    #pragma omp parallel
    {
        unordered_map<pair<int, int>, vector<double>, pair_hash> local_map;

        #pragma omp for nowait
        for (int i = 0; i < (int)orig.size(); i++) {
            const Edge& e = orig[i];
            local_map[make_pair(e.u, e.v)] = e.weights;
        }

        #pragma omp for nowait
        for (int i = 0; i < (int)ins.size(); i++) {
            const Edge& e = ins[i];
            local_map[make_pair(e.u, e.v)] = e.weights;
        }

        #pragma omp critical
        {
            for (auto& kv : local_map)
                edge_weights[kv.first] = kv.second;
        }
    }

    // auto end_step2 = high_resolution_clock::now();
    // cout << "Step 2 completed in " 
    //      << duration_cast<milliseconds>(end_step2 - start_step2).count() << " ms." << endl;

    // auto start_step3 = high_resolution_clock::now();
    // Step 3: Process S in parallel to fill comb and combE
    #pragma omp parallel for
    for (int idx = 0; idx < (int)S.size(); idx++) {
        auto it = std::next(S.begin(), idx);
        int u = it->first;
        int v = it->second;

        int count = 0;
        for (int i = 0; i < num_obj; i++) {
            if (parents[i][v] == u)
                count++;
        }

        vector<double> ow(num_obj, 1e9);
        auto found = edge_weights.find({u, v});
        if (found != edge_weights.end()) {
            ow = found->second;
        }

        bool hasAny = false;
        for (double w : ow) {
            if (w != 1e9) {
                hasAny = true;
                break;
            }
        }
        if (!hasAny) continue;

        double cw = (num_obj - count + 1);

        #pragma omp critical
        {
            comb[u].emplace_back(v, cw);
            combE[u].emplace_back(u, v, cw, ow);
        }
    }

    // auto end_step3 = high_resolution_clock::now();
    // cout << "Step 3 completed in " 
    //      << duration_cast<milliseconds>(end_step3 - start_step3).count() << " ms." << endl;
}

void mosp_update(
    vector<vector<vector<pair<int,double>>>>& graph,
    int source,
    vector<vector<int>>& parents,
    vector<vector<double>>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices, int num_obj,
    const vector<Edge>& orig_edges)
{
    // cout << "In mosp_update" << endl;
    #pragma omp parallel for
    for (auto &e : inserted_edges) {
        for (int i = 0; i < num_obj; i++) {
            if (e.weights[i] != INF)
                graph[i][e.u].push_back({e.v, e.weights[i]});
        }
    }

    // cout << "Running SOSP updates..." << endl;
    #pragma omp parallel for
    for (int i = 0; i < num_obj; i++) {
        sosp_update(graph[i], source,
                    parents[i], distances[i],
                    inserted_edges, num_vertices, i);
    }

    // cout << "Creating combined graph..." << endl;
    vector<vector<pair<int,double>>> comb(num_vertices);
    vector<vector<CombinedEdge>> combE(num_vertices);
    create_combined_graph(parents, num_vertices, num_obj, orig_edges, inserted_edges, comb, combE);

    // cout << "Running SOSP update on combined graph..." << endl;
    vector<int> cpar(num_vertices, -1);
    vector<double> cd(num_vertices, INF);
    sosp_update(comb, source, cpar, cd, {}, num_vertices, 0);

}

void initialize_graph(
    const string& filepath,
    vector<vector<vector<pair<int, double>>>>& graph,
    vector<Edge>& ins,
    int& N,
    int num_obj,
    vector<Edge>& orig,
    int num_changes)
{
    unordered_map<int, int> id_to_index;
    vector<pair<int, int>> edge_list;
    int index = 0;

    // Read from file
    ifstream file(filepath);
    string line;

    // ignore first 4 rows
    for (int i = 0; i < 4; ++i) {
        getline(file, line);
    }

    while (getline(file, line)) {
        stringstream ss(line);
        int u_id, v_id;
        ss >> u_id >> v_id;
        if (!id_to_index.count(u_id)) id_to_index[u_id] = index++;
        if (!id_to_index.count(v_id)) id_to_index[v_id] = index++;
        edge_list.emplace_back(id_to_index[u_id], id_to_index[v_id]);
    }

    N = index;
    graph.assign(num_obj, vector<vector<pair<int, double>>>(N));

    // Random number setup
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(1.0, 10.0);

    // Original edges
    for (auto [u, v] : edge_list) {
        vector<double> weights(num_obj);
        for (int i = 0; i < num_obj; ++i)
            weights[i] = dist(gen);
        orig.emplace_back(u, v, weights);
        for (int i = 0; i < num_obj; ++i)
            graph[i][u].emplace_back(v, weights[i]);
    }

    // Random insertions
    uniform_int_distribution<> node_dist(0, N - 1);
    for (int i = 0; i < num_changes; ++i) {
        int u = node_dist(gen);
        int v = node_dist(gen);
        while (u == v) v = node_dist(gen);
        vector<double> weights(num_obj);
        for (int j = 0; j < num_obj; ++j)
            weights[j] = dist(gen);
        ins.emplace_back(u, v, weights);
    }
}


// int main() {
//     // string dataset_path = "datasets/road_usa/road_usa.mtx";
//     // string dataset_path = "datasets/bio-CE-CX/bio-CE-CX.edges";
//     // string dataset_path = "datasets/GL7d16/GL7d16.mtx";
//     // string dataset_path = "datasets/3Dspectralwave2/3Dspectralwave2.mtx";
//     // string dataset_path = "datasets/test.data";
//     string dataset_path = "datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx";
//     int num_obj = 20, num_changes = 50, N, source = 1;
int main(int argc, char* argv[]) {
    if (argc < 5) {
        cout << "Usage: " << argv[0] << " <dataset_path> <num_obj> <num_changes> <source_node>" << endl;
        return 1;
    }

    string dataset_path = argv[1];
    int num_obj = stoi(argv[2]);
    int num_changes = stoi(argv[3]);
    int source = stoi(argv[4]);
    int N;

    cout<<num_obj<<","<<num_changes<<",";

    vector<Edge> inserted, original;
    vector<vector<vector<pair<int, double>>>> graph;

    // cout << "Initializing graph..." << endl;
    auto start_graph = high_resolution_clock::now();
    initialize_graph(dataset_path, graph, inserted, N, num_obj, original, num_changes);
    auto end_graph = high_resolution_clock::now();
    cout << duration_cast<milliseconds>(end_graph - start_graph).count() << ",";

    // cout << "Number of vertices: " << N << endl;
    // cout << "Number of edges: " << original.size() + inserted.size() << endl;
    // cout << "Number of objects: " << num_obj << endl;
    // cout << "Number of changes: " << inserted.size() << endl;
    // cout << "Number of inserted edges: " << inserted.size() << endl;
    // cout << "Number of original edges: " << original.size() << endl;

    vector<vector<int>> parents(num_obj, vector<int>(N));
    vector<vector<double>> distances(num_obj, vector<double>(N));

    // cout << "Running Dijkstra's algorithm..." << endl;
    auto start_dijkstra = high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < num_obj; i++) {
        dijkstra(graph[i], source, parents[i], distances[i], N);
    }
    auto end_dijkstra = high_resolution_clock::now();
    cout << duration_cast<milliseconds>(end_dijkstra - start_dijkstra).count() << ",";


    // cout << "Running MOSP update..." << endl;
    auto start_mosp = high_resolution_clock::now();
    mosp_update(graph, source, parents, distances, inserted, N, num_obj, original);
    auto end_mosp = high_resolution_clock::now();
    cout << duration_cast<milliseconds>(end_mosp - start_mosp).count() << ", Parallel" << endl;

    return 0;
}