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

void dijkstra(const vector<vector<pair<int, double>>>& graph, int source,
    vector<int>& parent, vector<double>& distances, int num_vertices) {
    if (source >= num_vertices || source < 0) {
        cerr << "Error: source node " << source << " out of bounds." << endl;
        return;
    }

    priority_queue<pair<double,int>, vector<pair<double,int>>, greater<>> pq;
    distances.assign(num_vertices, INF);
    parent.assign(num_vertices, -1);
    distances[source] = 0;
    pq.push({0, source});

    cout << "[Dijkstra] Started from node " << source << endl;

    int nodes_processed = 0, edges_checked = 0, updates = 0;

    while (!pq.empty()) {
        auto [dist,u] = pq.top(); pq.pop();
        if (dist > distances[u]) continue;

        nodes_processed++;
        for (auto &e : graph[u]) {
            int v = e.first; double w = e.second;
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

cout << "[Dijkstra] Done. Nodes: " << nodes_processed 
<< ", Edges: " << edges_checked 
<< ", Updates: " << updates << endl;
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
    vector<vector<pair<int,double>>> rev(num_vertices);
    #pragma omp parallel for schedule(dynamic)
    for(int u=0; u<num_vertices; u++){
        for(auto &pr: graph[u]){
            int v = pr.first; double w = pr.second;
            if (w == INF) continue;
            #pragma omp critical
            rev[v].push_back({u,w});
        }
    }

    vector<vector<pair<int,double>>> byV(num_vertices);
    for(auto &e: inserted_edges){
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
        for(int v=0; v<num_vertices; v++){
            for(auto &ue : byV[v]){
                int u = ue.first; double w = ue.second;
                if (w == INF) continue;
                double nd = distances[u] + w;
                if(nd < distances[v]){
                    #pragma omp critical
                    {
                        distances[v] = nd;
                        parent[v] = u;
                        if(!marked[v]){
                            marked[v] = 1;
                            affected.push_back(v);
                        }
                    }
                }
            }
        }
    }

    while(!affected.empty()){
        vector<char> inN(num_vertices, 0);
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<(int)affected.size(); i++){
            int v = affected[i];
            for(auto &vw: graph[v]){
                int nbr = vw.first;
                inN[nbr] = 1;
            }
        }

        vector<int> N;
        for(int v=0; v<num_vertices; v++){
            if(inN[v]) N.push_back(v);
        }

        vector<int> new_aff;
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<(int)N.size(); i++){
            int v = N[i];
            for(auto &uw: rev[v]){
                int u = uw.first; double w = uw.second;
                if (w == INF) continue;
                if(marked[u]){
                    double nd = distances[u] + w;
                    if(nd < distances[v]){
                        #pragma omp critical
                        {
                            distances[v] = nd;
                            parent[v] = u;
                            if(!marked[v]){
                                marked[v] = 1;
                                new_aff.push_back(v);
                            }
                        }
                    }
                }
            }
            for(auto &vw: graph[v]){
                int w = vw.first; double weight = vw.second;
                if (weight == INF) continue;
                double nd = distances[v] + weight;
                if(nd < distances[w]){
                    #pragma omp critical
                    {
                        distances[w] = nd;
                        parent[w] = v;
                        if(!marked[w]){
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
    vector<vector<pair<int,double>>>& comb,
    vector<vector<CombinedEdge>>& combE)
{
    comb.assign(num_vertices, {});
    combE.assign(num_vertices, {});
    set<pair<int,int>> S;

    for(int i = 0; i < num_obj; i++) {
        for(int v = 0; v < num_vertices; v++) {
            if (parents[i][v] >= 0) 
                S.insert({parents[i][v], v});
        }
    }

    for (auto &uv : S) {
        int u = uv.first, v = uv.second;
        int count = 0;
        for (int i = 0; i < num_obj; i++) {
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
    vector<vector<vector<pair<int,double>>>>& graph,
    int source,
    vector<vector<int>>& parents,
    vector<vector<double>>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices, int num_obj,
    const vector<Edge>& orig_edges)
{
    for (auto &e : inserted_edges) {
        for (int i = 0; i < num_obj; i++) {
            if (e.weights[i] != INF)
                graph[i][e.u].push_back({e.v, e.weights[i]});
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_obj; i++) {
        sosp_update(graph[i], source,
                    parents[i], distances[i],
                    inserted_edges, num_vertices, i);
    }

    vector<vector<pair<int,double>>> comb(num_vertices);
    vector<vector<CombinedEdge>> combE(num_vertices);
    create_combined_graph(parents, num_vertices, num_obj, orig_edges, inserted_edges, comb, combE);

    vector<int> cpar(num_vertices, -1);
    vector<double> cd(num_vertices, INF);
    sosp_update(comb, source, cpar, cd, {}, num_vertices, 0);

}


// void initialize_graph(
//     const string& filepath,
//     vector<vector<vector<pair<int, double>>>>& graph,
//     vector<Edge>& ins,
//     int& N,
//     int num_obj,
//     vector<Edge>& orig,
//     int num_changes)
// {
//     unordered_map<int, int> id_to_index;
//     int index = 0;
//     vector<pair<int, int>> edge_list;
//     ifstream file(filepath);
//     string line;
//     while (getline(file, line)) {
//         stringstream ss(line);
//         int u_id, v_id;
//         ss >> u_id >> v_id;
//         if (!id_to_index.count(u_id)) id_to_index[u_id] = index++;
//         if (!id_to_index.count(v_id)) id_to_index[v_id] = index++;
//         edge_list.emplace_back(id_to_index[u_id], id_to_index[v_id]);
//     }
//     N = index;
//     graph.assign(num_obj, vector<vector<pair<int, double>>>(N));
//     random_device rd;
//     mt19937 gen(rd());
//     uniform_real_distribution<> dist(1.0, 10.0);
//     for (auto [u, v] : edge_list) {
//         vector<double> weights(num_obj);
//         for (int i = 0; i < num_obj; ++i)
//             weights[i] = dist(gen);
//         orig.emplace_back(u, v, weights);
//         for (int i = 0; i < num_obj; ++i)
//             graph[i][u].emplace_back(v, weights[i]);
//     }
//     // Generate random inserted edges
//     uniform_int_distribution<> node_dist(0, N - 1);
//     for (int i = 0; i < num_changes; ++i) {
//         int u = node_dist(gen);
//         int v = node_dist(gen);
//         while (u == v) v = node_dist(gen); // avoid self-loop
//         vector<double> weights(num_obj);
//         for (int j = 0; j < num_obj; ++j)
//             weights[j] = dist(gen);
//         ins.emplace_back(u, v, weights);
//     }
// }

void initialize_graph(
    const string& filepath,
    vector<vector<vector<pair<int, double>>>>& graph,
    vector<Edge>& ins,
    int& N,
    int num_obj,
    vector<Edge>& orig,
    int num_changes)
{
    ifstream file(filepath);
    string line;
    // Skip header line
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }
    stringstream dim_ss(line);
    int num_edges;
    dim_ss >> num_edges >> num_obj;  // Assuming: <number of edges> <number of objectives>
    vector<double> all_weights;
    while (getline(file, line)) {
        if (line.empty()) continue;
        all_weights.push_back(stod(line));
    }
    N = 0;
    int edge_count = num_edges;
    int total_weights = all_weights.size();
    int vertices_guess = static_cast<int>(sqrt(edge_count)) + 1;
    graph.assign(num_obj, vector<vector<pair<int, double>>>(vertices_guess));
    // Generate edges with weights
    for (int i = 0; i < edge_count; ++i) {
        int u = i % vertices_guess;
        int v = (i + 1) % vertices_guess; // just to simulate v
        vector<double> weights(num_obj);
        for (int j = 0; j < num_obj; ++j) {
            int idx = i * num_obj + j;
            if (idx < total_weights)
                weights[j] = all_weights[idx];
            else
                weights[j] = 1.0; // default weight if missing
        }
        orig.emplace_back(u, v, weights);
        for (int j = 0; j < num_obj; ++j)
            graph[j][u].emplace_back(v, weights[j]);
        N = max(N, max(u, v) + 1);
    }
    // Generate random inserted edges
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(1.0, 10.0);
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



int main() {
    // string dataset_path = "datasets/road_usa/road_usa.mtx";
    string dataset_path = "datasets/road_usa/road_usa.mtx";
    // string dataset_path = "datasets/test.data";
    int num_obj = 2, num_changes = 5, N, source = 1;
    vector<Edge> inserted, original;
    vector<vector<vector<pair<int, double>>>> graph;

    cout << "Initializing graph..." << endl;
    auto start_graph = high_resolution_clock::now();
    initialize_graph(dataset_path, graph, inserted, N, num_obj, original, num_changes);
    auto end_graph = high_resolution_clock::now();
    cout << "Graph initialized in " 
         << duration_cast<milliseconds>(end_graph - start_graph).count() << " ms." << endl;


    cout << "Number of vertices: " << N << endl;
    cout << "Number of edges: " << original.size() + inserted.size() << endl;
    cout << "Number of objects: " << num_obj << endl;
    cout << "Number of changes: " << inserted.size() << endl;
    cout << "Number of inserted edges: " << inserted.size() << endl;
    cout << "Number of original edges: " << original.size() << endl;

    vector<vector<int>> parents(num_obj, vector<int>(N));
    vector<vector<double>> distances(num_obj, vector<double>(N));

    cout << "Running Dijkstra's algorithm..." << endl;
    auto start_dijkstra = high_resolution_clock::now();
    // #pragma omp parallel for schedule(dynamic)
    // for(int j = 0; j < original.size(); j++)
    for (int i = 0; i < num_obj; i++) {
        dijkstra(graph[i], source, parents[i], distances[i], N);
    }
    auto end_dijkstra = high_resolution_clock::now();
    cout << "Dijkstra's algorithm completed in " 
         << duration_cast<milliseconds>(end_dijkstra - start_dijkstra).count() << " ms." << endl;

    cout << "Running MOSP update..." << endl;
    auto start_mosp = high_resolution_clock::now();
    mosp_update(graph, source, parents, distances, inserted, N, num_obj, original);
    auto end_mosp = high_resolution_clock::now();
    cout << "MOSP update completed in " 
         << duration_cast<milliseconds>(end_mosp - start_mosp).count() << " ms." << endl;

    return 0;
}