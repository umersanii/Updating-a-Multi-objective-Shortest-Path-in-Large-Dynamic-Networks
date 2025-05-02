// MOSP_Parallel_Update.cpp
// Implements SOSP_Update and MOSP_Update based on the referenced paper

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <unordered_set>
#include <omp.h>

using namespace std;

const int INF = numeric_limits<int>::max();

struct Edge {
    int to;
    int weight; // For single-objective, extend to vector<int> for MOSP
};

struct NodeInfo {
    int distance = INF;
    int parent = -1;
    bool affected = false;
};

typedef vector<vector<Edge>> Graph;

void SOSP_Update(const Graph& G, vector<NodeInfo>& tree, const vector<pair<int, int>>& insertedEdges, const vector<int>& weights) {
    int V = G.size();
    vector<vector<pair<int, int>>> I(V); // I[v] = {(u, weight)}
    for (size_t i = 0; i < insertedEdges.size(); ++i) {
        auto [u, v] = insertedEdges[i];
        I[v].emplace_back(u, weights[i]);
    }

    vector<bool> marked(V, false);
    vector<int> affected;

    // Step 1: Process Changed Edges
    #pragma omp parallel for
    for (int v = 0; v < V; ++v) {
        for (auto [u, w] : I[v]) {
            if (tree[u].distance != INF && tree[v].distance > tree[u].distance + w) {
                #pragma omp critical
                {
                    tree[v].distance = tree[u].distance + w;
                    tree[v].parent = u;
                    tree[v].affected = true;
                    marked[v] = true;
                    affected.push_back(v);
                }
            }
        }
    }

    // Step 2: Propagate the Update
    while (!affected.empty()) {
        unordered_set<int> N;
        for (int u : affected) {
            for (const Edge& edge : G[u]) {
                N.insert(edge.to);
            }
        }

        vector<int> nextAffected;
        #pragma omp parallel for
        for (int i = 0; i < (int)N.size(); ++i) {
            int v = *next(N.begin(), i);
            for (const Edge& edge : G[v]) {
                int u = edge.to;
                if (!marked[u]) continue;
                if (tree[u].distance != INF && tree[v].distance > tree[u].distance + edge.weight) {
                    #pragma omp critical
                    {
                        tree[v].distance = tree[u].distance + edge.weight;
                        tree[v].parent = u;
                        tree[v].affected = true;
                        marked[v] = true;
                        nextAffected.push_back(v);
                    }
                }
            }
        }
        affected = nextAffected;
    }
}

// Example usage:
int main() {
    int V = 6;
    Graph G(V);

    G[0] = {{1, 2}, {2, 8}};   // u1
    G[1] = {{4, 5}};           // u2
    G[2] = {{4, 3}};           // u3
    G[3] = {{5, 2}};           // u4
    G[4] = {{5, 6}};           // u5
    // u6 has no outgoing edges

    vector<NodeInfo> tree(V);
    tree[0].distance = 0;

    vector<pair<int, int>> inserted = {{1, 2}, {2, 4}, {0, 4}};
    vector<int> weights = {7, 1, 4};

    SOSP_Update(G, tree, inserted, weights);

    cout << "Updated distances from source (u1):\n";
    for (int i = 0; i < V; ++i) {
        cout << "Node " << i << ": distance = " << tree[i].distance << ", parent = " << tree[i].parent << "\n";
    }

    return 0;
}
