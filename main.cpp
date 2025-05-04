#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>

using namespace std;

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
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq;
    distances.assign(num_vertices, INT_MAX);
    parent.assign(num_vertices, -1);
    distances[source] = 0;
    pq.push({0, source});
    
    while (!pq.empty()) {
        double dist = pq.top().first;
        int u = pq.top().second;
        pq.pop();
        
        if (dist > distances[u]) continue;
        
        for (const auto& [v, weight] : graph[u]) {
            if (distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
                parent[v] = u;
                pq.push({distances[v], v});
            }
        }
    }
}

void sosp_update(vector<vector<pair<int, double>>>& graph, int source, vector<int>& parent, 
                vector<double>& distances, const vector<Edge>& inserted_edges, int num_vertices, int objective_idx) {
    vector<vector<pair<int, double>>> grouped_edges(num_vertices);
    for (const auto& edge : inserted_edges) {
        grouped_edges[edge.v].push_back({edge.u, edge.weights[objective_idx]});
    }
    
    vector<int> affected;
    vector<int> marked(num_vertices, 0);
    
    for (int v = 0; v < num_vertices; ++v) {
        for (const auto& [u, weight] : grouped_edges[v]) {
            double new_dist = distances[u] + weight;
            if (new_dist < distances[v]) {
                distances[v] = new_dist;
                parent[v] = u;
                affected.push_back(v);
                marked[v] = 1;
            }
        }
    }
    
    while (!affected.empty()) {
        set<int> neighbors;
        for (int v : affected) {
            for (const auto& [u, w] : graph[v]) {
                neighbors.insert(u);
            }
        }
        
        vector<int> new_affected;
        for (int v : neighbors) {
            for (const auto& [u, weight] : graph[v]) {
                if (marked[u] == 1) {
                    double new_dist = distances[u] + weight;
                    if (new_dist < distances[v]) {
                        distances[v] = new_dist;
                        parent[v] = u;
                        marked[v] = 1;
                        new_affected.push_back(v);
                    }
                }
            }
        }
        affected = move(new_affected);
    }
}

void initialize_graph(vector<vector<vector<pair<int, double>>>>& graph, vector<Edge>& inserted_edges, 
                     int& num_vertices, int num_objectives, vector<Edge>& original_edges) {
    original_edges = {
        Edge(0, 1, {5.0, 10.0}),
        Edge(0, 2, {2.0, 8.0}),
        Edge(1, 3, {3.0, 4.0}),
        Edge(2, 3, {4.0, 3.0}),
        Edge(2, 4, {6.0, 5.0}),
        Edge(3, 5, {2.0, 6.0}),
        Edge(4, 5, {3.0, 2.0})
    };

    num_vertices = 6;
    graph.resize(num_objectives, vector<vector<pair<int, double>>>(num_vertices));
    for (const auto& edge : original_edges) {
        for (int i = 0; i < num_objectives; ++i) {
            graph[i][edge.u].push_back({edge.v, edge.weights[i]});
        }
    }

    inserted_edges = {
        Edge(1, 4, {2.0, 3.0}),
        Edge(0, 3, {3.5, 7.0})
    };
}

void create_combined_graph(const vector<vector<int>>& parents, int num_vertices, int num_objectives,
                          const vector<Edge>& original_edges, const vector<Edge>& inserted_edges,
                          vector<vector<pair<int, double>>>& combined_graph,
                          vector<vector<CombinedEdge>>& combined_edges) {
    combined_graph.resize(num_vertices);
    combined_edges.resize(num_vertices);
    
    set<tuple<int, int>> edge_set;
    for (int i = 0; i < num_objectives; ++i) {
        for (int v = 0; v < num_vertices; ++v) {
            if (parents[i][v] != -1) {
                edge_set.insert({parents[i][v], v});
            }
        }
    }
    
    for (const auto& [u, v] : edge_set) {
        int count = 0;
        vector<double> original_weights(num_objectives, -1.0);
        for (int i = 0; i < num_objectives; ++i) {
            if (parents[i][v] == u) {
                count++;
            }
        }
        // Check original edges
        for (const auto& edge : original_edges) {
            if (edge.u == u && edge.v == v) {
                original_weights = edge.weights;
                break;
            }
        }
        // Check inserted edges
        for (const auto& edge : inserted_edges) {
            if (edge.u == u && edge.v == v) {
                original_weights = edge.weights;
                break;
            }
        }
        if (original_weights[0] == -1.0) {
            cerr << "Warning: No weights found for edge (" << u << "," << v << ")\n";
            continue;
        }
        double combined_weight = num_objectives - count + 1;
        combined_graph[u].push_back({v, combined_weight});
        combined_edges[u].push_back(CombinedEdge(u, v, combined_weight, original_weights));
    }
}

void mosp_update(vector<vector<vector<pair<int, double>>>>& graph, int source, 
                 vector<vector<int>>& parents, vector<vector<double>>& distances,
                 const vector<Edge>& inserted_edges, int num_vertices, int num_objectives,
                 const vector<Edge>& original_edges) {
    // Update graph with inserted edges
    for (const auto& edge : inserted_edges) {
        for (int i = 0; i < num_objectives; ++i) {
            graph[i][edge.u].push_back({edge.v, edge.weights[i]});
        }
    }
    
    // Update SOSP trees
    for (int i = 0; i < num_objectives; ++i) {
        sosp_update(graph[i], source, parents[i], distances[i], inserted_edges, num_vertices, i);
    }
    
    vector<vector<pair<int, double>>> combined_graph;
    vector<vector<CombinedEdge>> combined_edges;
    create_combined_graph(parents, num_vertices, num_objectives, original_edges, inserted_edges, combined_graph, combined_edges);
    
    vector<int> combined_parent;
    vector<double> combined_distances;
    dijkstra(combined_graph, source, combined_parent, combined_distances, num_vertices);
    
    cout << "\nMOSP Path (vertex: [weights for objectives]):\n";
    for (int v = 0; v < num_vertices; ++v) {
        if (combined_parent[v] != -1 || v == source) {
            cout << "Vertex " << v << ": ";
            if (v != source) {
                int u = combined_parent[v];
                bool found = false;
                for (const auto& edge : combined_edges[u]) {
                    if (edge.v == v) {
                        cout << "[";
                        for (int i = 0; i < num_objectives; ++i) {
                            cout << edge.original_weights[i];
                            if (i < num_objectives - 1) cout << ", ";
                        }
                        cout << "]";
                        found = true;
                        break;
                    }
                }
                if (!found) cout << "[not found]";
            }
            cout << "\n";
        }
    }
}

int main() {
    vector<vector<vector<pair<int, double>>>> graph;
    vector<Edge> inserted_edges;
    vector<Edge> original_edges;
    int num_vertices;
    int source = 0;
    int num_objectives = 2;
    vector<vector<int>> parents(num_objectives);
    vector<vector<double>> distances(num_objectives);
    
    initialize_graph(graph, inserted_edges, num_vertices, num_objectives, original_edges);
    
    for (int i = 0; i < num_objectives; ++i) {
        parents[i].resize(num_vertices);
        distances[i].resize(num_vertices);
        dijkstra(graph[i], source, parents[i], distances[i], num_vertices);
    }
    
    cout << fixed << setprecision(6);
    for (int i = 0; i < num_objectives; ++i) {
        cout << "\nInitial Distances (Objective " << i + 1 << "):\n";
        for (int j = 0; j < num_vertices; ++j) {
            if (distances[i][j] != INT_MAX) {
                cout << "Vertex " << j << ": " << distances[i][j] << "\n";
            }
        }
        cout << "Initial Parents (Objective " << i + 1 << "):\n";
        for (int j = 0; j < num_vertices; ++j) {
            if (parents[i][j] != -1 || j == source) {
                cout << "Vertex " << j << ": " << parents[i][j] << "\n";
            }
        }
    }
    
    mosp_update(graph, source, parents, distances, inserted_edges, num_vertices, num_objectives, original_edges);
    
    for (int i = 0; i < num_objectives; ++i) {
        cout << "\nUpdated Distances (Objective " << i + 1 << "):\n";
        for (int j = 0; j < num_vertices; ++j) {
            if (distances[i][j] != INT_MAX) {
                cout << "Vertex " << j << ": " << distances[i][j] << "\n";
            }
        }
        cout << "Updated Parents (Objective " << i + 1 << "):\n";
        for (int j = 0; j < num_vertices; ++j) {
            if (parents[i][j] != -1 || j == source) {
                cout << "Vertex " << j << ": " << parents[i][j] << "\n";
            }
        }
    }
    
    return 0;
}