#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <climits>
#include <iomanip>
#include <omp.h>
#include <limits>

using namespace std;

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
    priority_queue<pair<double,int>, vector<pair<double,int>>, greater<>> pq;
    distances.assign(num_vertices, INF);
    parent.assign(num_vertices, -1);
    distances[source] = 0;
    pq.push({0, source});
    while (!pq.empty()) {
        auto [dist,u] = pq.top(); pq.pop();
        if (dist > distances[u]) continue;
        for (auto &e : graph[u]) {
            int v = e.first; double w = e.second;
            if (w == INF) continue;
            if (distances[u] + w < distances[v]) {
                distances[v] = distances[u] + w;
                parent[v] = u;
                pq.push({distances[v], v});
            }
        }
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

    //////////////////////// debugging ////////////////////////
    cout << "\nAfter update:\n";
    for (int i = 0; i < num_obj; i++) {
        cout << "Obj " << i + 1 << ":\n";
        for (int v = 0; v < num_vertices; v++) {
            if (distances[i][v] != INF)
                cout << " v" << v + 1 << ":" << distances[i][v] << " from " << parents[i][v] + 1 << "\n";
        }
    }
    //////////////////////////////////////////////////////////

    vector<vector<pair<int,double>>> comb(num_vertices);
    vector<vector<CombinedEdge>> combE(num_vertices);
    create_combined_graph(parents, num_vertices, num_obj, orig_edges, inserted_edges, comb, combE);

    ///////////////////////// debugging //////////
    cout << "\nCombined graph edges:\n";
    for (int u = 0; u < num_vertices; u++) {
        for (auto &ce : combE[u]) {
            cout << "Edge: " << ce.u + 1 << "->" << ce.v + 1 << ": ";
            int max_w = 0;
            for (int i = 0; i < num_obj; i++) {
                if (ce.original_weights[i] != INF)
                    max_w = max(max_w, (int)ce.original_weights[i]);
            }
            cout << ce.weight;
            cout << "\n";
        }
    }
    ///////////////////////////////////////////////

    vector<int> cpar(num_vertices, -1);
    vector<double> cd(num_vertices, INF);
    sosp_update(comb, source, cpar, cd, {}, num_vertices, 0);

    /////////////////// debugging //////////
    cout << "\nCombined graph distances (SOSP):\n";
    for (int v = 0; v < num_vertices; v++) {
        if (cd[v] != INF)
            cout << "Edge: " << cpar[v] + 1 << "->" <<  v + 1 << " with " <<cd[v]<< "\n";
    }
    ///////////////////////////////////////
}

void initialize_graph(
    vector<vector<vector<pair<int,double>>>>& graph,
    vector<Edge>& ins, int& N, int num_obj,
    vector<Edge>& orig)
{
    //// example 1 ////
    orig = {
        {0, 1, {2, INF}},   
        {0, 2, {7, 1}},   
        {1, 3, {3, 4}},   
        {3, 4, {3, 2}},   
        {3, 5, {1, 7}},   
        {2, 1, {INF, 2}}, 
        {1, 4, {INF, 2}}, 
        {4, 5, {INF, 2}}
    };

    //// example 2 ////
    // orig = {
    //     {0, 2, {8, INF}},   
    //     {2, 1, {2, 1}},   
    //     {1, 3, {3, 4}},   
    //     {1, 4, {9, 2}},   
    //     {3, 5, {2, 7}},
    //     {4, 3, {2, 2}}, 
    //     {4, 5, {6, 2}}, 
    // };

    N = 6;
    graph.assign(num_obj, vector<vector<pair<int,double>>>(N));
    for(auto &e: orig)
      for(int i=0;i<num_obj;i++)
        if(e.weights[i] != INF)
          graph[i][e.u].push_back({e.v,e.weights[i]});
    

    //// example 1 ////      
    ins = {
    };
    //// example 2 ////
    // ins = {
    //     {0,1,{7,INF}},
    //     {0,4,{4,INF}},
    //     {2,4,{1,3}}
    // };
}

int main(){
    int source=0, num_obj=2, N;
    vector<Edge> inserted, original;
    vector<vector<vector<pair<int,double>>>> graph;
    initialize_graph(graph, inserted, N, num_obj, original);

    ////////////////////////////// debugging ////////////////////////////
    cout << "Original graph:\n";
    for (int i = 0; i < num_obj; i++) {
        cout << "Obj " << i + 1 << ":\n";
        for (int u = 0; u < N; u++) {
            cout << "Vertex " << u << ": ";
            for (auto &vp : graph[i][u]) {
                cout << "[" << vp.first << "," << vp.second << "] ";
            }
            cout << "\n";
        }
    }
    cout << "\nInserted edges:\n";
    for (auto &e : inserted) {
        cout << "Edge: " << e.u << "->" << e.v << ": ";
        for (int i = 0; i < num_obj; i++) {
            if (e.weights[i] != INF)
                cout << e.weights[i] << " ";
            else
                cout << "INF ";
        }
        cout << "\n";
    }
    cout << "\nOriginal edges:\n";
    for (auto &e : original) {
        cout << "Edge: " << e.u + 1<< "->" << e.v + 1<< ": ";
        for (int i = 0; i < num_obj; i++) {
            if (e.weights[i] != INF)
                cout << e.weights[i] << " ";
            else
                cout << "INF ";
        }
        cout << "\n";
    }
    //////////////////////////////////////////////////////////////////////

    vector<vector<int>> parents(num_obj, vector<int>(N));
    vector<vector<double>> distances(num_obj, vector<double>(N));

    for(int i=0;i<num_obj;i++){
        dijkstra(graph[i], source, parents[i], distances[i], N);
    }

    mosp_update(graph, source, parents, distances,
                inserted, N, num_obj, original);


    return 0;
}