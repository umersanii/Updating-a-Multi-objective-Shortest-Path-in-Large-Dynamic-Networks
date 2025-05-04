#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>

using namespace std;

// Structure to represent an edge with multiple objectives
struct Edge {
    int u, v;
    vector<double> weights; // Weight vector for multiple objectives
    Edge(int _u, int _v, const vector<double>& _w) : u(_u), v(_v), weights(_w) {}
};

// Structure to represent an edge in the combined graph
struct CombinedEdge {
    int u, v;
    double weight; // Single weight for combined graph
    vector<double> original_weights; // Original multi-objective weights
    CombinedEdge(int _u, int _v, double _w, const vector<double>& _ow) 
        : u(_u), v(_v), weight(_w), original_weights(_ow) {}
};

// Function to compute initial SOSP tree using Dijkstra's algorithm for a single objective
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

// Function to update SOSP tree for a single objective
void sosp_update(vector<vector<pair<int, double>>>& graph, int source, vector<int>& parent, 
                vector<double>& distances, const vector<Edge>& inserted_edges, int num_vertices, int objective_idx) {
    // Step 0: Preprocessing - Group inserted edges by destination vertex
    vector<vector<pair<int, double>>> grouped_edges(num_vertices);
    for (const auto& edge : inserted_edges) {
        grouped_edges[edge.v].push_back({edge.u, edge.weights[objective_idx]});
    }
    
    // Step 1: Process Changed Edges
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
    
    // Step 2: Propagate the Update
    while (!affected.empty()) {
        // Gather unique neighbors of affected vertices
        set<int> neighbors;
        for (int v : affected) {
            for (const auto& [u, w] : graph[v]) {
                neighbors.insert(u);
            }
        }
        
        // Prepare for next iteration
        vector<int> new_affected;
        
        // Process each neighbor
        for (int v : neighbors) {
            for (const auto& [u, weight] : graph[v]) { // Predecessors of v
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
        
        // Update affected list for next iteration
        affected = move(new_affected);
    }
}

// Function to initialize graph with hardcoded dataset
void initialize_graph(vector<vector<vector<pair<int, double>>>>& graph, vector<Edge>& inserted_edges, 
                     int& num_vertices, int num_objectives) {
    // Hardcoded dataset: 6 vertices, edges with two objectives (time, fuel)
    vector<Edge> edges = {
        Edge(0, 1, {5.0, 10.0}),  // 0 -> 1: 5 time, 10 fuel
        Edge(0, 2, {2.0, 8.0}),   // 0 -> 2: 2 time, 8 fuel
        Edge(1, 3, {3.0, 4.0}),   // 1 -> 3: 3 time, 4 fuel
        Edge(2, 3, {4.0, 3.0}),   // 2 -> 3: 4 time, 3 fuel
        Edge(2, 4, {6.0, 5.0}),   // 2 -> 4: 6 time, 5 fuel
        Edge(3, 5, {2.0, 6.0}),   // 3 -> 5: 2 time, 6 fuel
        Edge(4, 5, {3.0, 2.0})    // 4 -> 5: 3 time, 2 fuel
    };

    // Determine number of vertices
    num_vertices = 6; // Vertices 0 to 5

    // Initialize graph for each objective
    graph.resize(num_objectives, vector<vector<pair<int, double>>>(num_vertices));
    for (const auto& edge : edges) {
        for (int i = 0; i < num_objectives; ++i) {
            graph[i][edge.u].push_back({edge.v, edge.weights[i]});
        }
    }

    // Hardcoded inserted edges for dynamic update
    inserted_edges = {
        Edge(1, 4, {2.0, 3.0}), // New edge 1 -> 4: 2 time, 3 fuel
        Edge(0, 3, {3.5, 7.0})  // New edge 0 -> 3: 3.5 time, 7 fuel
    };
}

// Function to create combined graph
void create_combined_graph(const vector<vector<int>>& parents, int num_vertices, int num_objectives,
                          const vector<vector<vector<pair<int, double>>>>& original_graph,
                          vector<vector<pair<int, double>>>& combined_graph,
                          vector<vector<CombinedEdge>>& combined_edges) {
    combined_graph.resize(num_vertices);
    combined_edges.resize(num_vertices);
    
    // Collect edges from all SOSP trees
    set<tuple<int, int>> edge_set;
    for (int i = 0; i < num_objectives; ++i) {
        for (int v = 0; v < num_vertices; ++v) {
            if (parents[i][v] != -1) {
                int u = parents[i][v];
                edge_set.insert({u, v});
            }
        }
    }
    
    // Assign weights to edges in combined graph
    for (const auto& [u, v] : edge_set) {
        int count = 0;
        vector<double> original_weights(num_objectives);
        for (int i = 0; i < num_objectives; ++i) {
            // Check if edge (u,v) is in SOSP tree i
            if (parents[i][v] == u) {
                count++;
                original_weights[i] = 0;
                for (const auto& [dest, w] : original_graph[i][u]) {
                    if (dest == v) {
                        original_weights[i] = w;
                        break;
                    }
                }
            } else {
                // Use original weight if not in SOSP tree
                for (const auto& [dest, w] : original_graph[i][u]) {
                    if (dest == v) {
                        original_weights[i] = w;
                        break;
                    }
                }
            }
        }
        double combined_weight = num_objectives - count + 1;
        combined_graph[u].push_back({v, combined_weight});
        combined_edges[u].push_back(CombinedEdge(u, v, combined_weight, original_weights));
    }
}

// Function to compute MOSP
void mosp_update(vector<vector<vector<pair<int, double>>>>& graph, int source, 
                 vector<vector<int>>& parents, vector<vector<double>>& distances,
                 const vector<Edge>& inserted_edges, int num_vertices, int num_objectives) {
    // Step 1: Update SOSP trees for each objective
    for (int i = 0; i < num_objectives; ++i) {
        sosp_update(graph[i], source, parents[i], distances[i], inserted_edges, num_vertices, i);
    }
    
    // Step 2: Create combined graph
    vector<vector<pair<int, double>>> combined_graph;
    vector<vector<CombinedEdge>> combined_edges;
    create_combined_graph(parents, num_vertices, num_objectives, graph, combined_graph, combined_edges);
    
    // Step 3: Find SOSP in combined graph
    vector<int> combined_parent;
    vector<double> combined_distances;
    dijkstra(combined_graph, source, combined_parent, combined_distances, num_vertices);
    
    // Reassign original weights to get MOSP
    cout << "\nMOSP Path (vertex: [weights for objectives]):\n";
    for (int v = 0; v < num_vertices; ++v) {
        if (combined_parent[v] != -1 || v == source) {
            cout << "Vertex " << v << ": ";
            if (v != source) {
                int u = combined_parent[v];
                for (const auto& edge : combined_edges[u]) {
                    if (edge.v == v) {
                        cout << "[";
                        for (int i = 0; i < num_objectives; ++i) {
                            cout << edge.original_weights[i];
                            if (i < num_objectives - 1) cout << ", ";
                        }
                        cout << "]";
                        break;
                    }
                }
            }
            cout << "\n";
        }
    }
}

int main() {
    // Initialize variables
    vector<vector<vector<pair<int, double>>>> graph; // Graph for each objective
    vector<Edge> inserted_edges;
    int num_vertices;
    int source = 0; // Choose vertex 0 as source
    int num_objectives = 2; // Two objectives (time, fuel)
    vector<vector<int>> parents(num_objectives);
    vector<vector<double>> distances(num_objectives);
    
    // Initialize graph with hardcoded dataset
    initialize_graph(graph, inserted_edges, num_vertices, num_objectives);
    
    // Compute initial SOSP trees for each objective
    for (int i = 0; i < num_objectives; ++i) {
        parents[i].resize(num_vertices);
        distances[i].resize(num_vertices);
        dijkstra(graph[i], source, parents[i], distances[i], num_vertices);
    }
    
    // Print initial state
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
    
    // Update MOSP
    mosp_update(graph, source, parents, distances, inserted_edges, num_vertices, num_objectives);
    
    // Print updated state
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