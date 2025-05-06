#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>
#include <omp.h>
#include <limits>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_map>
#include <chrono>
#include <mpi.h>
#include <mpi.h>

using namespace std;
using namespace std::chrono;

const double INF = numeric_limits<double>::infinity();

// Struct definitions
// Struct definitions
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
        : u(_u), v(_v), weight(_w), original_weights(_ow) {}
};

// Function definitions
// Function definitions
void dijkstra(const vector<vector<pair<int, double>>>& graph, int source,
              vector<int>& parent, vector<double>& distances, int num_vertices) {
              vector<int>& parent, vector<double>& distances, int num_vertices) {
    if (source >= num_vertices || source < 0) {
        cerr << "Error: source node " << source << " out of bounds." << endl;
        return;
    }

    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq;
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq;
    distances.assign(num_vertices, INF);
    parent.assign(num_vertices, -1);
    distances[source] = 0;
    pq.push({0, source});

    cout << "[Dijkstra] Started from node " << source << endl;

    int nodes_processed = 0, edges_checked = 0, updates = 0;

    while (!pq.empty()) {
        auto [dist, u] = pq.top();
        pq.pop();
        auto [dist, u] = pq.top();
        pq.pop();
        if (dist > distances[u]) continue;

        nodes_processed++;
        for (auto &e : graph[u]) {
            int v = e.first;
            double w = e.second;
            int v = e.first;
            double w = e.second;
            edges_checked++;
            if (w == INF) continue;
            if (distances[u] + w < distances[v]) {
                distances[v] = distances[u] + w;
                parent[v] = u;
                pq.push({distances[v], v});
                updates++;
            }
        }
    }
        }
    }

    cout << "[Dijkstra] Done. Nodes: " << nodes_processed
         << ", Edges: " << edges_checked
         << ", Updates: " << updates << endl;
    cout << "[Dijkstra] Done. Nodes: " << nodes_processed
         << ", Edges: " << edges_checked
         << ", Updates: " << updates << endl;
}

void sosp_update(
    vector<vector<pair<int, double>>>& graph,
    vector<vector<pair<int, double>>>& graph,
    int source,
    vector<int>& parent,
    vector<double>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices,
    int obj_idx)
{
    vector<vector<pair<int, double>>> rev(num_vertices);
    vector<vector<pair<int, double>>> rev(num_vertices);
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < num_vertices; u++) {
        for (auto &pr : graph[u]) {
            int v = pr.first;
            double w = pr.second;
    for (int u = 0; u < num_vertices; u++) {
        for (auto &pr : graph[u]) {
            int v = pr.first;
            double w = pr.second;
            if (w == INF) continue;
            #pragma omp critical
            rev[v].push_back({u, w});
            rev[v].push_back({u, w});
        }
    }

    vector<vector<pair<int, double>>> byV(num_vertices);
    for (auto &e : inserted_edges) {
    vector<vector<pair<int, double>>> byV(num_vertices);
    for (auto &e : inserted_edges) {
        if (e.weights[obj_idx] != INF)
            byV[e.v].push_back({e.u, e.weights[obj_idx]});
    }

    vector<int> affected;
    vector<char> marked(num_vertices, 0);

    if (inserted_edges.empty()) {
        distances.assign(num_vertices, INF);
        parent.assign(num_vertices, -1);
        distances[source] = 0;
        marked[source] = 1;
        affected.push_back(source);
        for (auto &pr : graph[source]) {
            int v = pr.first;
            if (pr.second != INF && !marked[v]) {
                marked[v] = 1;
                affected.push_back(v);
            }
        }
    } else {
        #pragma omp parallel for schedule(dynamic)
        for (int v = 0; v < num_vertices; v++) {
            for (auto &ue : byV[v]) {
                int u = ue.first;
                double w = ue.second;
        for (int v = 0; v < num_vertices; v++) {
            for (auto &ue : byV[v]) {
                int u = ue.first;
                double w = ue.second;
                if (w == INF) continue;
                double nd = distances[u] + w;
                if (nd < distances[v]) {
                if (nd < distances[v]) {
                    #pragma omp critical
                    {
                        distances[v] = nd;
                        parent[v] = u;
                        if (!marked[v]) {
                        if (!marked[v]) {
                            marked[v] = 1;
                            affected.push_back(v);
                        }
                    }
                }
            }
        }
    }

    while (!affected.empty()) {
    while (!affected.empty()) {
        vector<char> inN(num_vertices, 0);
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)affected.size(); i++) {
        for (int i = 0; i < (int)affected.size(); i++) {
            int v = affected[i];
            for (auto &vw : graph[v]) {
            for (auto &vw : graph[v]) {
                int nbr = vw.first;
                inN[nbr] = 1;
            }
        }

        vector<int> N;
        for (int v = 0; v < num_vertices; v++) {
            if (inN[v]) N.push_back(v);
        for (int v = 0; v < num_vertices; v++) {
            if (inN[v]) N.push_back(v);
        }

        vector<int> new_aff;
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)N.size(); i++) {
        for (int i = 0; i < (int)N.size(); i++) {
            int v = N[i];
            for (auto &uw : rev[v]) {
                int u = uw.first;
                double w = uw.second;
            for (auto &uw : rev[v]) {
                int u = uw.first;
                double w = uw.second;
                if (w == INF) continue;
                if (marked[u]) {
                if (marked[u]) {
                    double nd = distances[u] + w;
                    if (nd < distances[v]) {
                    if (nd < distances[v]) {
                        #pragma omp critical
                        {
                            distances[v] = nd;
                            parent[v] = u;
                            if (!marked[v]) {
                            if (!marked[v]) {
                                marked[v] = 1;
                                new_aff.push_back(v);
                            }
                        }
                    }
                }
            }
            for (auto &vw : graph[v]) {
                int w = vw.first;
                double weight = vw.second;
            for (auto &vw : graph[v]) {
                int w = vw.first;
                double weight = vw.second;
                if (weight == INF) continue;
                double nd = distances[v] + weight;
                if (nd < distances[w]) {
                if (nd < distances[w]) {
                    #pragma omp critical
                    {
                        distances[w] = nd;
                        parent[w] = v;
                        if (!marked[w]) {
                        if (!marked[w]) {
                            marked[w] = 1;
                            new_aff.push_back(w);
                        }
                    }
                }
            }
        }
        affected.swap(new_aff);
    }
}

void create_combined_graph(
    const vector<vector<int>>& parents,
    int num_vertices, int num_obj,
    const vector<Edge>& orig, const vector<Edge>& ins,
    vector<vector<pair<int, double>>>& comb,
    vector<vector<pair<int, double>>>& comb,
    vector<vector<CombinedEdge>>& combE)
{
    comb.assign(num_vertices, {});
    combE.assign(num_vertices, {});
    set<pair<int, int>> S;
    set<pair<int, int>> S;

    for (int i = 0; i < num_obj; i++) {
        for (int v = 0; v < num_vertices; v++) {
            if (parents[i][v] >= 0)
    for (int i = 0; i < num_obj; i++) {
        for (int v = 0; v < num_vertices; v++) {
            if (parents[i][v] >= 0)
                S.insert({parents[i][v], v});
        }
    }

    for (auto &uv : S) {
        int u = uv.first, v = uv.second;
        int count = 0;
        for (int i = 0; i < num_obj; i++) { // Fixed typo: 'pans' to 'i < num_obj'
        for (int i = 0; i < num_obj; i++) { // Fixed typo: 'pans' to 'i < num_obj'
            if (parents[i][v] == u) count++;
        }

        vector<double> ow(num_obj, INF);
        for (auto &e : orig) {
            if (e.u == u && e.v == v) {
                for (int i = 0; i < num_obj; i++) {
                    if (e.weights[i] != INF)
                        ow[i] = e.weights[i];
                }
                break;
            }
        }
        for (auto &e : ins) {
            if (e.u == u && e.v == v) {
                for (int i = 0; i < num_obj; i++) {
                    if (e.weights[i] != INF)
                        ow[i] = e.weights[i];
                }
                break;
            }
        }

        bool hasAny = false;
        for (double w : ow) {
            if (w != INF) { hasAny = true; break; }
        }
        if (!hasAny) continue;

        double cw = (num_obj - count + 1);
        comb[u].push_back({v, cw});
        combE[u].push_back(CombinedEdge(u, v, cw, ow));
    }
}

void mosp_update(
    vector<vector<vector<pair<int, double>>>>& graph,
    vector<vector<vector<pair<int, double>>>>& graph,
    int source,
    vector<vector<int>>& parents,
    vector<vector<double>>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices, int num_obj,
    const vector<Edge>& orig_edges)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Distribute objectives across MPI processes
    int objs_per_process = num_obj / size;
    int start_obj = rank * objs_per_process;
    int end_obj = (rank == size - 1) ? num_obj : start_obj + objs_per_process;

    if (rank == 0) {
        cout << "[MOSP] Distributing " << num_obj << " objectives across " << size << " processes." << endl;
        cout << "[MOSP] Process " << rank << " handles objectives " << start_obj << " to " << end_obj - 1 << endl;
    }

    // Update graph with inserted edges
    for (const auto &e : inserted_edges) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Distribute objectives across MPI processes
    int objs_per_process = num_obj / size;
    int start_obj = rank * objs_per_process;
    int end_obj = (rank == size - 1) ? num_obj : start_obj + objs_per_process;

    if (rank == 0) {
        cout << "[MOSP] Distributing " << num_obj << " objectives across " << size << " processes." << endl;
        cout << "[MOSP] Process " << rank << " handles objectives " << start_obj << " to " << end_obj - 1 << endl;
    }

    // Update graph with inserted edges
    for (const auto &e : inserted_edges) {
        for (int i = 0; i < num_obj; i++) {
            if (e.weights[i] != INF)
                graph[i][e.u].push_back({e.v, e.weights[i]});
        }
    }

    // Process assigned objectives
    for (int i = start_obj; i < end_obj; i++) {
        sosp_update(graph[i], source, parents[i], distances[i], inserted_edges, num_vertices, i);
    }

    // Gather distances and parents
    vector<double> all_distances(num_obj * num_vertices, INF);
    vector<int> all_parents(num_obj * num_vertices, -1);

    vector<double> local_distances((end_obj - start_obj) * num_vertices);
    vector<int> local_parents((end_obj - start_obj) * num_vertices);

    for (int i = start_obj; i < end_obj; i++) {
        int offset = (i - start_obj) * num_vertices;
        for (int v = 0; v < num_vertices; v++) {
            local_distances[offset + v] = distances[i][v];
            local_parents[offset + v] = parents[i][v];
        }
    }

    MPI_Gather(
        local_distances.data(), (end_obj - start_obj) * num_vertices, MPI_DOUBLE,
        all_distances.data(), (end_obj - start_obj) * num_vertices, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    MPI_Gather(
        local_parents.data(), (end_obj - start_obj) * num_vertices, MPI_INT,
        all_parents.data(), (end_obj - start_obj) * num_vertices, MPI_INT,
        0, MPI_COMM_WORLD);

    // Broadcast combined results
    MPI_Bcast(all_distances.data(), num_obj * num_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(all_parents.data(), num_obj * num_vertices, MPI_INT, 0, MPI_COMM_WORLD);

    // Process assigned objectives
    for (int i = start_obj; i < end_obj; i++) {
        sosp_update(graph[i], source, parents[i], distances[i], inserted_edges, num_vertices, i);
    }

    // Gather distances and parents
    vector<double> all_distances(num_obj * num_vertices, INF);
    vector<int> all_parents(num_obj * num_vertices, -1);

    vector<double> local_distances((end_obj - start_obj) * num_vertices);
    vector<int> local_parents((end_obj - start_obj) * num_vertices);

    for (int i = start_obj; i < end_obj; i++) {
        int offset = (i - start_obj) * num_vertices;
        for (int v = 0; v < num_vertices; v++) {
            local_distances[offset + v] = distances[i][v];
            local_parents[offset + v] = parents[i][v];
        }
    }

    MPI_Gather(
        local_distances.data(), (end_obj - start_obj) * num_vertices, MPI_DOUBLE,
        all_distances.data(), (end_obj - start_obj) * num_vertices, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    MPI_Gather(
        local_parents.data(), (end_obj - start_obj) * num_vertices, MPI_INT,
        all_parents.data(), (end_obj - start_obj) * num_vertices, MPI_INT,
        0, MPI_COMM_WORLD);

    // Broadcast combined results
    MPI_Bcast(all_distances.data(), num_obj * num_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(all_parents.data(), num_obj * num_vertices, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < num_obj; i++) {
        int offset = i * num_vertices;
        for (int v = 0; v < num_vertices; v++) {
            distances[i][v] = all_distances[offset + v];
            parents[i][v] = all_parents[offset + v];
        }
    }

    // Combined graph processing
    vector<vector<pair<int, double>>> comb(num_vertices);
        int offset = i * num_vertices;
        for (int v = 0; v < num_vertices; v++) {
            distances[i][v] = all_distances[offset + v];
            parents[i][v] = all_parents[offset + v];
        }
    }

    // Combined graph processing
    vector<vector<pair<int, double>>> comb(num_vertices);
    vector<vector<CombinedEdge>> combE(num_vertices);
    create_combined_graph(parents, num_vertices, num_obj, orig_edges, inserted_edges, comb, combE);

    vector<int> cpar(num_vertices, -1);
    vector<double> cd(num_vertices, INF);
    sosp_update(comb, source, cpar, cd, {}, num_vertices, 0);

    if (rank == 0) {
        cout << "[MOSP] Update completed." << endl;
    }
}
    if (rank == 0) {
        cout << "[MOSP] Update completed." << endl;
    }
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
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << "Initializing graph..." << endl;

        ifstream file(filepath);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filepath << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        unordered_map<long long, int> id_to_index;
        int index = 0;
        vector<pair<int, int>> edge_list;
        string line;

        while (getline(file, line)) {
            if (line.empty() || line[0] == '%') continue;
            stringstream ss(line);
            long long u_id, v_id;
            if (ss >> u_id >> v_id) {
                if (!id_to_index.count(u_id)) id_to_index[u_id] = index++;
                if (!id_to_index.count(v_id)) id_to_index[v_id] = index++;
                edge_list.emplace_back(id_to_index[u_id], id_to_index[v_id]);
            }
        }
        file.close();

        N = index;
        cout << "Number of vertices: " << N << ", Number of edges: " << edge_list.size() << endl;

        graph.assign(num_obj, vector<vector<pair<int, double>>>(N));
        orig.reserve(edge_list.size());

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dist(1.0, 10.0);

        for (const auto& [u, v] : edge_list) {
            vector<double> weights(num_obj);
            for (int i = 0; i < num_obj; ++i) {
                weights[i] = dist(gen);
            }
            orig.emplace_back(u, v, weights);
            for (int i = 0; i < num_obj; ++i) {
                graph[i][u].emplace_back(v, weights[i]);
            }
        }

        uniform_int_distribution<> node_dist(0, N - 1);
        ins.reserve(num_changes);
        for (int i = 0; i < num_changes; ++i) {
            int u = node_dist(gen);
            int v = node_dist(gen);
            while (u == v) v = node_dist(gen);
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; ++j) {
                weights[j] = dist(gen);
            }
            ins.emplace_back(u, v, weights);
        }

        cout << "Graph initialized successfully." << endl;
    }

    // Broadcast N
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast original edges
    int orig_size = (rank == 0) ? orig.size() : 0;
    MPI_Bcast(&orig_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> orig_u(orig_size), orig_v(orig_size);
    vector<double> orig_weights(orig_size * num_obj);
    if (rank == 0) {
        for (int i = 0; i < orig_size; i++) {
            orig_u[i] = orig[i].u;
            orig_v[i] = orig[i].v;
            for (int j = 0; j < num_obj; j++) {
                orig_weights[i * num_obj + j] = orig[i].weights[j];
            }
        }
    }
    MPI_Bcast(orig_u.data(), orig_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(orig_v.data(), orig_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(orig_weights.data(), orig_size * num_obj, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        orig.clear();
        orig.reserve(orig_size);
        for (int i = 0; i < orig_size; i++) {
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; j++) {
                weights[j] = orig_weights[i * num_obj + j];
            }
            orig.emplace_back(orig_u[i], orig_v[i], weights);
        }
    }

    // Broadcast inserted edges
    int ins_size = (rank == 0) ? ins.size() : 0;
    MPI_Bcast(&ins_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> ins_u(ins_size), ins_v(ins_size);
    vector<double> ins_weights(ins_size * num_obj);
    if (rank == 0) {
        for (int i = 0; i < ins_size; i++) {
            ins_u[i] = ins[i].u;
            ins_v[i] = ins[i].v;
            for (int j = 0; j < num_obj; j++) {
                ins_weights[i * num_obj + j] = ins[i].weights[j];
            }
        }
    }
    MPI_Bcast(ins_u.data(), ins_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ins_v.data(), ins_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ins_weights.data(), ins_size * num_obj, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        ins.clear();
        ins.reserve(ins_size);
        for (int i = 0; i < ins_size; i++) {
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; j++) {
                weights[j] = ins_weights[i * num_obj + j];
            }
            ins.emplace_back(ins_u[i], ins_v[i], weights);
        }
    }

    // Initialize graph on non-root processes
    if (rank != 0) {
        graph.assign(num_obj, vector<vector<pair<int, double>>>(N));
        for (const auto& e : orig) {
            for (int i = 0; i < num_obj; i++) {
                graph[i][e.u].emplace_back(e.v, e.weights[i]);
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << "Initializing graph..." << endl;

        ifstream file(filepath);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filepath << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        unordered_map<long long, int> id_to_index;
        int index = 0;
        vector<pair<int, int>> edge_list;
        string line;

        while (getline(file, line)) {
            if (line.empty() || line[0] == '%') continue;
            stringstream ss(line);
            long long u_id, v_id;
            if (ss >> u_id >> v_id) {
                if (!id_to_index.count(u_id)) id_to_index[u_id] = index++;
                if (!id_to_index.count(v_id)) id_to_index[v_id] = index++;
                edge_list.emplace_back(id_to_index[u_id], id_to_index[v_id]);
            }
        }
        file.close();

        N = index;
        cout << "Number of vertices: " << N << ", Number of edges: " << edge_list.size() << endl;

        graph.assign(num_obj, vector<vector<pair<int, double>>>(N));
        orig.reserve(edge_list.size());

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dist(1.0, 10.0);

        for (const auto& [u, v] : edge_list) {
            vector<double> weights(num_obj);
            for (int i = 0; i < num_obj; ++i) {
                weights[i] = dist(gen);
            }
            orig.emplace_back(u, v, weights);
            for (int i = 0; i < num_obj; ++i) {
                graph[i][u].emplace_back(v, weights[i]);
            }
        }

        uniform_int_distribution<> node_dist(0, N - 1);
        ins.reserve(num_changes);
        for (int i = 0; i < num_changes; ++i) {
            int u = node_dist(gen);
            int v = node_dist(gen);
            while (u == v) v = node_dist(gen);
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; ++j) {
                weights[j] = dist(gen);
            }
            ins.emplace_back(u, v, weights);
        }

        cout << "Graph initialized successfully." << endl;
    }

    // Broadcast N
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast original edges
    int orig_size = (rank == 0) ? orig.size() : 0;
    MPI_Bcast(&orig_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> orig_u(orig_size), orig_v(orig_size);
    vector<double> orig_weights(orig_size * num_obj);
    if (rank == 0) {
        for (int i = 0; i < orig_size; i++) {
            orig_u[i] = orig[i].u;
            orig_v[i] = orig[i].v;
            for (int j = 0; j < num_obj; j++) {
                orig_weights[i * num_obj + j] = orig[i].weights[j];
            }
        }
    }
    MPI_Bcast(orig_u.data(), orig_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(orig_v.data(), orig_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(orig_weights.data(), orig_size * num_obj, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        orig.clear();
        orig.reserve(orig_size);
        for (int i = 0; i < orig_size; i++) {
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; j++) {
                weights[j] = orig_weights[i * num_obj + j];
            }
            orig.emplace_back(orig_u[i], orig_v[i], weights);
        }
    }

    // Broadcast inserted edges
    int ins_size = (rank == 0) ? ins.size() : 0;
    MPI_Bcast(&ins_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> ins_u(ins_size), ins_v(ins_size);
    vector<double> ins_weights(ins_size * num_obj);
    if (rank == 0) {
        for (int i = 0; i < ins_size; i++) {
            ins_u[i] = ins[i].u;
            ins_v[i] = ins[i].v;
            for (int j = 0; j < num_obj; j++) {
                ins_weights[i * num_obj + j] = ins[i].weights[j];
            }
        }
    }
    MPI_Bcast(ins_u.data(), ins_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ins_v.data(), ins_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ins_weights.data(), ins_size * num_obj, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        ins.clear();
        ins.reserve(ins_size);
        for (int i = 0; i < ins_size; i++) {
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; j++) {
                weights[j] = ins_weights[i * num_obj + j];
            }
            ins.emplace_back(ins_u[i], ins_v[i], weights);
        }
    }

    // Initialize graph on non-root processes
    if (rank != 0) {
        graph.assign(num_obj, vector<vector<pair<int, double>>>(N));
        for (const auto& e : orig) {
            for (int i = 0; i < num_obj; i++) {
                graph[i][e.u].emplace_back(e.v, e.weights[i]);
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    string dataset_path = "datasets/road_usa/road_usa.mtx";
    int num_obj = 2, num_changes = 5, N, source = 1;
    vector<Edge> inserted, original;
    vector<vector<vector<pair<int, double>>>> graph;

    if (rank == 0) {
        cout << "Initializing graph..." << endl;
    }
    if (rank == 0) {
        cout << "Initializing graph..." << endl;
    }
    auto start_graph = high_resolution_clock::now();
    initialize_graph(dataset_path, graph, inserted, N, num_obj, original, num_changes);
    auto end_graph = high_resolution_clock::now();
    if (rank == 0) {
        cout << "Graph initialized in "
             << duration_cast<milliseconds>(end_graph - start_graph).count() << " ms." << endl;
        cout << "Number of vertices: " << N << endl;
        cout << "Number of edges: " << (original.size() + inserted.size()) << endl;
        cout << "Number of objects: " << num_obj << endl;
        cout << "Number of changes: " << inserted.size() << endl;
    }
    if (rank == 0) {
        cout << "Graph initialized in "
             << duration_cast<milliseconds>(end_graph - start_graph).count() << " ms." << endl;
        cout << "Number of vertices: " << N << endl;
        cout << "Number of edges: " << (original.size() + inserted.size()) << endl;
        cout << "Number of objects: " << num_obj << endl;
        cout << "Number of changes: " << inserted.size() << endl;
    }

    vector<vector<int>> parents(num_obj, vector<int>(N));
    vector<vector<double>> distances(num_obj, vector<double>(N));

    if (rank == 0) {
        cout << "Running Dijkstra's algorithm..." << endl;
    }
    if (rank == 0) {
        cout << "Running Dijkstra's algorithm..." << endl;
    }
    auto start_dijkstra = high_resolution_clock::now();
    for (int i = 0; i < num_obj; i++) {
        dijkstra(graph[i], source, parents[i], distances[i], N);
    }
    auto end_dijkstra = high_resolution_clock::now();
    if (rank == 0) {
        cout << "Dijkstra's algorithm completed in "
             << duration_cast<milliseconds>(end_dijkstra - start_dijkstra).count() << " ms." << endl;
    }
    if (rank == 0) {
        cout << "Dijkstra's algorithm completed in "
             << duration_cast<milliseconds>(end_dijkstra - start_dijkstra).count() << " ms." << endl;
    }

    if (rank == 0) {
        cout << "Running MOSP update..." << endl;
    }
    if (rank == 0) {
        cout << "Running MOSP update..." << endl;
    }
    auto start_mosp = high_resolution_clock::now();
    mosp_update(graph, source, parents, distances, inserted, N, num_obj, original);
    auto end_mosp = high_resolution_clock::now();
    if (rank == 0) {
        cout << "MOSP update completed in "
             << duration_cast<milliseconds>(end_mosp - start_mosp).count() << " ms." << endl;
    }
    if (rank == 0) {
        cout << "MOSP update completed in "
             << duration_cast<milliseconds>(end_mosp - start_mosp).count() << " ms." << endl;
    }

    MPI_Finalize();
    MPI_Finalize();
    return 0;
}