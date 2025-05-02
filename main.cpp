#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

// Structure to represent an edge
struct Edge {
    int u, v;
    double weight;
    Edge(int _u, int _v, double _w) : u(_u), v(_v), weight(_w) {}
};

// Function to compute initial SOSP tree using Dijkstra's algorithm
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

// Function to update SOSP tree
void sosp_update(vector<vector<pair<int, double>>>& graph, int source, vector<int>& parent, 
                vector<double>& distances, const vector<Edge>& inserted_edges, int num_vertices) {
    // Step 0: Preprocessing - Group inserted edges by destination vertex
    vector<vector<pair<int, double>>> grouped_edges(num_vertices);
    for (const auto& edge : inserted_edges) {
        grouped_edges[edge.v].push_back({edge.u, edge.weight});
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

// Function to initialize graph from an .edge file
void initialize_graph(vector<vector<pair<int, double>>>& graph, vector<Edge>& inserted_edges, int& num_vertices, const string& filename) {
    // Open the .edge file
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    // Read edges from the file
    vector<Edge> edges;
    int u, v;
    double weight;
    int max_vertex = -1;

    string line;
    while (file >> u >> v >> weight) {
        edges.emplace_back(u, v, weight);
        max_vertex = max(max_vertex, max(u, v));
    }

    file.close();

    if (edges.empty()) {
        cerr << "Error: No edges found in the file" << endl;
        exit(1);
    }

    // Determine number of vertices (max vertex ID + 1)
    num_vertices = max_vertex + 1;

    // Initialize graph
    graph.resize(num_vertices);
    for (const auto& edge : edges) {
        graph[edge.u].push_back({edge.v, edge.weight});
    }

    // Example inserted edges for dynamic update (same as in original code)
    inserted_edges = {
        Edge(0, 2, 2.0), // New edge from 0 to 2
        Edge(1, 4, 1.5)  // New edge from 1 to 4
    };
}

int main() {
    // Initialize variables
    vector<vector<pair<int, double>>> graph;
    vector<Edge> inserted_edges;
    int num_vertices;
    int source = 0; // Choose vertex 0 as source
    vector<int> parent;
    vector<double> distances;
    
    // Specify the .edge file
    string filename = "graph.edge"; // Replace with your .edge file path

    // Initialize graph from the .edge file
    initialize_graph(graph, inserted_edges, num_vertices, filename);
    
    // Compute initial SOSP tree
    dijkstra(graph, source, parent, distances, num_vertices);
    
    // Print initial state
    cout << fixed << setprecision(6);
    cout << "Initial Distances:\n";
    for (int i = 0; i < num_vertices; ++i) {
        if (distances[i] != INT_MAX) {
            cout << "Vertex " << i << ": " << distances[i] << "\n";
        }
    }
    cout << "Initial Parents:\n";
    for (int i = 0; i < num_vertices; ++i) {
        if (parent[i] != -1 || i == source) {
            cout << "Vertex " << i << ": " << parent[i] << "\n";
        }
    }
    
    // Update SOSP tree
    sosp_update(graph, source, parent, distances, inserted_edges, num_vertices);
    
    // Print updated state
    cout << "\nUpdated Distances:\n";
    for (int i = 0; i < num_vertices; ++i) {
        if (distances[i] != INT_MAX) {
            cout << "Vertex " << i << ": " << distances[i] << "\n";
        }
    }
    cout << "Updated Parents:\n";
    for (int i = 0; i < num_vertices; ++i) {
        if (parent[i] != -1 || i == source) {
            cout << "Vertex " << i << ": " << parent[i] << "\n";
        }
    }
    
    return 0;
}