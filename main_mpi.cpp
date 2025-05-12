#include <mpi.h>
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>
#include <limits>
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
    Edge() : u(-1), v(-1), weights() {}
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
    const vector<vector<pair<int, double> > >& graph,
    int source,
    vector<int>& parent,
    vector<double>& dist,
    int N)
{
    dist.assign(N, numeric_limits<double>::infinity());
    parent.assign(N, -1);
    dist[source] = 0.0;

    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P> > pq;
    pq.push(make_pair(0.0, source));

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << "Rank " << rank << ": Starting Dijkstra for source=" << source << endl;

    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u]) continue;

        for (size_t i = 0; i < graph[u].size(); ++i) {
            int v = graph[u][i].first;
            double w = graph[u][i].second;
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.push(make_pair(dist[v], v));
            }
        }
    }
    cout << "Rank " << rank << ": Finished Dijkstra" << endl;
}

void sosp_update(
    vector<vector<pair<int, double> > >& graph,
    int source,
    vector<int>& parent,
    vector<double>& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices,
    int obj_idx)
{
    vector<vector<pair<int, double> > > rev(num_vertices);
    for (int u = 0; u < num_vertices; u++) {
        for (size_t i = 0; i < graph[u].size(); i++) {
            int v = graph[u][i].first;
            double w = graph[u][i].second;
            if (w == INF) continue;
            rev[v].push_back(make_pair(u, w));
        }
    }

    vector<vector<pair<int, double> > > byV(num_vertices);
    for (size_t i = 0; i < inserted_edges.size(); i++) {
        const Edge& e = inserted_edges[i];
        if (e.weights[obj_idx] != INF)
            byV[e.v].push_back(make_pair(e.u, e.weights[obj_idx]));
    }

    vector<int> affected;
    vector<char> marked(num_vertices, 0);

    if (inserted_edges.empty()) {
        distances.assign(num_vertices, INF);
        parent.assign(num_vertices, -1);
        distances[source] = 0;
        marked[source] = 1;
        affected.push_back(source);
        for (size_t i = 0; i < graph[source].size(); i++) {
            int v = graph[source][i].first;
            double w = graph[source][i].second;
            if (w != INF && !marked[v]) {
                marked[v] = 1;
                affected.push_back(v);
            }
        }
    } else {
        for (int v = 0; v < num_vertices; v++) {
            for (size_t i = 0; i < byV[v].size(); i++) {
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
        for (size_t i = 0; i < affected.size(); i++) {
            int v = affected[i];
            for (size_t j = 0; j < graph[v].size(); j++) {
                int nbr = graph[v][j].first;
                inN[nbr] = 1;
            }
        }

        vector<int> N;
        for (int v = 0; v < num_vertices; v++) {
            if (inN[v]) N.push_back(v);
        }

        vector<int> new_aff;
        for (size_t i = 0; i < N.size(); i++) {
            int v = N[i];
            for (size_t j = 0; j < rev[v].size(); j++) {
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
            for (size_t j = 0; j < graph[v].size(); j++) {
                int w = graph[v][j].first;
                double weight = graph[v][j].second;
                if (weight == INF) continue;
                double nd = distances[v] + weight;
                if (nd < distances[w]) {
                    distances[w] = nd;
                    parent[w] = v;
                    if (!marked[w]) {
                        marked[w] = 1;
                        new_aff.push_back(w);
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

void create_combined_graph(
    const vector<vector<int> >& parents,
    int num_vertices,
    int num_obj,
    const vector<Edge>& orig,
    const vector<Edge>& ins,
    vector<vector<pair<int, double> > >& comb,
    vector<vector<CombinedEdge> >& combE)
{
    comb.assign(num_vertices, vector<pair<int, double> >());
    combE.assign(num_vertices, vector<CombinedEdge>());
    set<pair<int, int> > S;

    for (int i = 0; i < num_obj; i++) {
        for (int v = 0; v < num_vertices; v++) {
            if (parents[i][v] >= 0)
                S.insert(make_pair(parents[i][v], v));
        }
    }

    unordered_map<pair<int, int>, vector<double>, pair_hash> edge_weights;

    for (size_t i = 0; i < orig.size(); i++) {
        const Edge& e = orig[i];
        edge_weights[make_pair(e.u, e.v)] = e.weights;
    }
    for (size_t i = 0; i < ins.size(); i++) {
        const Edge& e = ins[i];
        edge_weights[make_pair(e.u, e.v)] = e.weights;
    }

    for (set<pair<int, int> >::iterator it = S.begin(); it != S.end(); ++it) {
        int u = it->first;
        int v = it->second;

        int count = 0;
        for (int i = 0; i < num_obj; i++) {
            if (parents[i][v] == u) count++;
        }

        vector<double> ow(num_obj, 1e9);
        unordered_map<pair<int, int>, vector<double>, pair_hash>::iterator wit = edge_weights.find(make_pair(u, v));
        if (wit != edge_weights.end()) {
            ow = wit->second;
        }

        bool hasAny = false;
        for (size_t i = 0; i < ow.size(); i++) {
            if (ow[i] != 1e9) {
                hasAny = true;
                break;
            }
        }
        if (!hasAny) continue;

        double cw = (num_obj - count + 1);
        comb[u].push_back(make_pair(v, cw));
        combE[u].push_back(CombinedEdge(u, v, cw, ow));
    }
}

void mosp_update(
    vector<vector<vector<pair<int, double> > > >& graph,
    int source,
    vector<vector<int> >& parents,
    vector<vector<double> >& distances,
    const vector<Edge>& inserted_edges,
    int num_vertices,
    int num_obj,
    const vector<Edge>& orig_edges)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << "Rank " << rank << ": Starting mosp_update" << endl;

    for (size_t i = 0; i < inserted_edges.size(); i++) {
        const Edge& e = inserted_edges[i];
        for (int j = 0; j < num_obj; j++) {
            if (e.weights[j] != INF)
                graph[j][e.u].push_back(make_pair(e.v, e.weights[j]));
        }
    }

    for (int i = 0; i < num_obj; i++) {
        sosp_update(graph[i], source,
                    parents[i], distances[i],
                    inserted_edges, num_vertices, i);
    }

    vector<vector<pair<int, double> > > comb(num_vertices);
    vector<vector<CombinedEdge> > combE(num_vertices);
    create_combined_graph(parents, num_vertices, num_obj, orig_edges, inserted_edges, comb, combE);

    vector<int> cpar(num_vertices, -1);
    vector<double> cd(num_vertices, INF);
    sosp_update(comb, source, cpar, cd, vector<Edge>(), num_vertices, 0);

    cout << "Rank " << rank << ": Finished mosp_update" << endl;
}

void broadcast_edges(vector<Edge>& edges, int rank, MPI_Comm comm) {
    cout << "Rank " << rank << ": Starting broadcast_edges, edges.size()=" << edges.size() << endl;

    // Broadcast number of edges
    int size = edges.size();
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Rank " << rank << ": Broadcasted size=" << size << endl;

    // Validate size
    if (size < 0 || size > 1000000) {
        cout << "Rank " << rank << ": Error: Invalid edge size=" << size << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Broadcast weight sizes for each edge
    vector<int> weight_sizes(size);
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            weight_sizes[i] = edges[i].weights.size();
            if (weight_sizes[i] <= 0 || weight_sizes[i] > 1000) {
                cout << "Rank " << rank << ": Error: Invalid weight_size[" << i << "]=" << weight_sizes[i] << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        cout << "Rank " << rank << ": weight_sizes={";
        for (int i = 0; i < min(5, size); ++i) cout << weight_sizes[i] << ",";
        if (size > 5) cout << "...";
        cout << "}" << endl;
    }
    MPI_Bcast(weight_sizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Rank " << rank << ": Broadcasted weight_sizes, first=" << (size > 0 ? weight_sizes[0] : -1) << endl;

    // Validate weight sizes on non-root ranks
    if (rank != 0) {
        for (int i = 0; i < size; ++i) {
            if (weight_sizes[i] <= 0 || weight_sizes[i] > 1000) {
                cout << "Rank " << rank << ": Error: Received invalid weight_size[" << i << "]=" << weight_sizes[i] << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // Compute total weights
    int total_weights = 0;
    for (int i = 0; i < size; ++i) {
        total_weights += weight_sizes[i];
    }
    cout << "Rank " << rank << ": total_weights=" << total_weights << endl;

    // Validate total_weights
    if (total_weights < 0 || total_weights > 1000000) {
        cout << "Rank " << rank << ": Error: Invalid total_weights=" << total_weights << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Resize edges on non-root ranks
    if (rank != 0) {
        edges.resize(size);
        for (int i = 0; i < size; ++i) {
            edges[i].weights.resize(weight_sizes[i]);
        }
        cout << "Rank " << rank << ": Resized edges to " << size << endl;
    }

    // Broadcast u, v, and weights
    vector<int> u_vec(size), v_vec(size);
    vector<double> weights_flat(total_weights);
    if (rank == 0) {
        int weight_idx = 0;
        for (int i = 0; i < size; ++i) {
            u_vec[i] = edges[i].u;
            v_vec[i] = edges[i].v;
            for (size_t j = 0; j < edges[i].weights.size(); j++) {
                weights_flat[weight_idx++] = edges[i].weights[j];
            }
        }
        cout << "Rank " << rank << ": u_vec={";
        for (int i = 0; i < min(5, size); ++i) cout << u_vec[i] << ",";
        if (size > 5) cout << "...";
        cout << "}, v_vec={";
        for (int i = 0; i < min(5, size); ++i) cout << v_vec[i] << ",";
        if (size > 5) cout << "...";
        cout << "}" << endl;
    }

    MPI_Bcast(u_vec.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(v_vec.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    if (total_weights > 0) {
        MPI_Bcast(weights_flat.data(), total_weights, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    cout << "Rank " << rank << ": Broadcasted u_vec, v_vec, weights_flat" << endl;

    // Reconstruct edges on non-root ranks
    if (rank != 0) {
        int weight_idx = 0;
        for (int i = 0; i < size; ++i) {
            edges[i].u = u_vec[i];
            edges[i].v = v_vec[i];
            for (int j = 0; j < weight_sizes[i]; ++j) {
                if (weight_idx >= total_weights) {
                    cout << "Rank " << rank << ": Error: weight_idx=" << weight_idx << " exceeds total_weights=" << total_weights << endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                edges[i].weights[j] = weights_flat[weight_idx++];
            }
        }
        cout << "Rank " << rank << ": Reconstructed edges" << endl;
    }

    cout << "Rank " << rank << ": Finished broadcast_edges" << endl;
}

void broadcast_graph(vector<vector<vector<pair<int, double> > > >& graph, int num_obj, int N, int rank, MPI_Comm comm) {
    cout << "Rank " << rank << ": Starting broadcast_graph, num_obj=" << num_obj << ", N=" << N << endl;

    // Validate inputs
    if (num_obj <= 0 || N <= 0) {
        cout << "Rank " << rank << ": Error: Invalid num_obj=" << num_obj << " or N=" << N << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank != 0) {
        graph.assign(num_obj, vector<vector<pair<int, double> > >(N));
        cout << "Rank " << rank << ": Allocated graph" << endl;
    }

    for (int i = 0; i < num_obj; ++i) {
        vector<int> edge_counts(N);
        vector<int> v_flat;
        vector<double> w_flat;
        if (rank == 0) {
            for (int u = 0; u < N; ++u) {
                edge_counts[u] = graph[i][u].size();
                if (edge_counts[u] < 0 || edge_counts[u] > 1000000) {
                    cout << "Rank " << rank << ": Error: Invalid edge_counts[" << u << "]=" << edge_counts[u] << endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                for (size_t j = 0; j < graph[i][u].size(); j++) {
                    int v = graph[i][u][j].first;
                    double w = graph[i][u][j].second;
                    v_flat.push_back(v);
                    w_flat.push_back(w);
                }
            }
            cout << "Rank " << rank << ": Prepared edge_counts, v_flat.size()=" << v_flat.size() << endl;
        }

        MPI_Bcast(edge_counts.data(), N, MPI_INT, 0, comm);
        cout << "Rank " << rank << ": Broadcasted edge_counts" << endl;

        // Validate edge_counts on non-root ranks
        if (rank != 0) {
            for (int u = 0; u < N; ++u) {
                if (edge_counts[u] < 0 || edge_counts[u] > 1000000) {
                    cout << "Rank " << rank << ": Error: Received invalid edge_counts[" << u << "]=" << edge_counts[u] << endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        int total_edges = 0;
        if (rank == 0) {
            total_edges = v_flat.size();
        }
        MPI_Bcast(&total_edges, 1, MPI_INT, 0, comm);
        cout << "Rank " << rank << ": Broadcasted total_edges=" << total_edges << endl;

        // Validate total_edges
        if (total_edges < 0 || total_edges > 1000000) {
            cout << "Rank " << rank << ": Error: Invalid total_edges=" << total_edges << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (rank != 0) {
            v_flat.resize(total_edges);
            w_flat.resize(total_edges);
            cout << "Rank " << rank << ": Resized v_flat, w_flat to " << total_edges << endl;
        }

        MPI_Bcast(v_flat.data(), total_edges, MPI_INT, 0, comm);
        MPI_Bcast(w_flat.data(), total_edges, MPI_DOUBLE, 0, comm);
        cout << "Rank " << rank << ": Broadcasted v_flat, w_flat" << endl;

        if (rank != 0) {
            int idx = 0;
            for (int u = 0; u < N; ++u) {
                graph[i][u].resize(edge_counts[u]);
                for (int j = 0; j < edge_counts[u]; ++j) {
                    if (idx >= total_edges) {
                        cout << "Rank " << rank << ": Error: idx=" << idx << " exceeds total_edges=" << total_edges << endl;
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    graph[i][u][j] = make_pair(v_flat[idx], w_flat[idx]);
                    ++idx;
                }
            }
            cout << "Rank " << rank << ": Reconstructed graph[" << i << "]" << endl;
        }
    }

    cout << "Rank " << rank << ": Finished broadcast_graph" << endl;
}

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "Rank " << rank << ": Started main, argc=" << argc << endl;

    if (argc < 4) {
        if (rank == 0)
            cout << "Usage: " << argv[0] << " <num_obj> <num_changes> <source_node>" << endl;
        MPI_Finalize();
        return 1;
    }

    int num_obj = stoi(argv[1]);
    int num_changes = stoi(argv[2]);
    int source = stoi(argv[3]);
    int N = 0;

    cout << "Rank " << rank << ": Parsed args: num_obj=" << num_obj
         << ", num_changes=" << num_changes << ", source=" << source << endl;

    // Validate inputs
    if (num_obj <= 0 || num_changes < 0) {
        if (rank == 0)
            cout << "Error: num_obj must be positive, num_changes must be non-negative" << endl;
        MPI_Finalize();
        return 1;
    }

    vector<Edge> inserted, original;
    vector<vector<vector<pair<int, double> > > > graph;

    auto start_graph = high_resolution_clock::now();
    if (rank == 0) {
        cout << "Rank " << rank << ": Initializing hardcoded dataset" << endl;

        // Hardcode dataset equivalent to small_data.txt
        N = 4; // 4 nodes: 0, 1, 2, 3
        vector<pair<int, int> > edge_list = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}
        };

        // Initialize random number generator with fixed seed for reproducibility
        mt19937 gen(42); // Fixed seed
        uniform_real_distribution<> dist(1.0, 10.0);

        // Original edges
        graph.assign(num_obj, vector<vector<pair<int, double> > >(N));
        for (size_t i = 0; i < edge_list.size(); i++) {
            int u = edge_list[i].first;
            int v = edge_list[i].second;
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; ++j)
                weights[j] = dist(gen);
            original.push_back(Edge(u, v, weights));
            for (int j = 0; j < num_obj; ++j)
                graph[j][u].push_back(make_pair(v, weights[j]));
            cout << "Rank " << rank << ": Added edge (" << u << "," << v << ") with weights";
            for (double w : weights) cout << " " << w;
            cout << endl;
        }

        // Random inserted edges
        uniform_int_distribution<> node_dist(0, N - 1);
        for (int i = 0; i < num_changes; ++i) {
            int u = node_dist(gen);
            int v = node_dist(gen);
            while (u == v) v = node_dist(gen);
            vector<double> weights(num_obj);
            for (int j = 0; j < num_obj; ++j)
                weights[j] = dist(gen);
            inserted.push_back(Edge(u, v, weights));
            cout << "Rank " << rank << ": Added inserted edge (" << u << "," << v << ") with weights";
            for (double w : weights) cout << " " << w;
            cout << endl;
        }

        cout << "Rank " << rank << ": Dataset initialized: N=" << N
             << ", original.size()=" << original.size()
             << ", inserted.size()=" << inserted.size() << endl;
    }

    cout << "Rank " << rank << ": Broadcasting N and num_obj" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_obj, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Rank " << rank << ": Received N=" << N << ", num_obj=" << num_obj << endl;

    // Validate source node
    if (source >= N || source < 0) {
        if (rank == 0)
            cout << "Error: Source node " << source << " is out of range (N = " << N << ")" << endl;
        MPI_Finalize();
        return 1;
    }

    // Broadcast graph, inserted, and original
    broadcast_graph(graph, num_obj, N, rank, MPI_COMM_WORLD);
    broadcast_edges(inserted, rank, MPI_COMM_WORLD);
    broadcast_edges(original, rank, MPI_COMM_WORLD);

    // Log broadcasted data
    cout << "Rank " << rank << ": N=" << N << ", num_obj=" << num_obj
         << ", graph.size()=" << graph.size() << ", inserted.size()=" << inserted.size()
         << ", original.size()=" << original.size() << endl;

    auto end_graph = high_resolution_clock::now();
    if (rank == 0)
        cout << num_obj << "," << num_changes << "," << duration_cast<milliseconds>(end_graph - start_graph).count() << ",";

    // Distribute objectives among processes
    int chunk_size = (num_obj + size - 1) / size;
    int start_idx = rank * chunk_size;
    int end_idx = min(start_idx + chunk_size, num_obj);

    // Allocate local_parents and local_distances only if there are objectives to process
    vector<vector<int> > local_parents;
    vector<vector<double> > local_distances;
    if (start_idx < num_obj) {
        local_parents.assign(end_idx - start_idx, vector<int>(N));
        local_distances.assign(end_idx - start_idx, vector<double>(N));
    }
    cout << "Rank " << rank << ": Allocated local_parents.size()=" << local_parents.size()
         << ", start_idx=" << start_idx << ", end_idx=" << end_idx << endl;

    auto start_dijkstra = high_resolution_clock::now();
    for (int i = start_idx; i < end_idx; i++) {
        local_distances[i - start_idx].assign(N, INF);
        dijkstra(graph[i], source, local_parents[i - start_idx], local_distances[i - start_idx], N);
    }
    auto end_dijkstra = high_resolution_clock::now();

    // Gather results at rank 0
    vector<vector<int> > parents(num_obj, vector<int>(N));
    vector<vector<double> > distances(num_obj, vector<double>(N));

    // Broadcast results from each rank
    for (int i = 0; i < num_obj; ++i) {
        int src_rank = i / chunk_size;
        if (rank == src_rank && i >= start_idx && i < end_idx) {
            parents[i] = local_parents[i - start_idx];
            distances[i] = local_distances[i - start_idx];
        }
        MPI_Bcast(parents[i].data(), N, MPI_INT, src_rank, MPI_COMM_WORLD);
        MPI_Bcast(distances[i].data(), N, MPI_DOUBLE, src_rank, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        auto dijkstra_time = duration_cast<milliseconds>(end_dijkstra - start_dijkstra).count();
        cout << dijkstra_time << ",";

        auto start_mosp = high_resolution_clock::now();
        mosp_update(graph, source, parents, distances, inserted, N, num_obj, original);
        auto end_mosp = high_resolution_clock::now();
        cout << duration_cast<milliseconds>(end_mosp - start_mosp).count() << ", MPI" << endl;
    }

    MPI_Finalize();
    return 0;
}