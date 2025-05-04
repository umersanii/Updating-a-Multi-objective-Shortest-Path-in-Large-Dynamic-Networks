# Updating-a-Multi-objective-Shortest-Path-in-Large-Dynamic-Networks


- **Initialization**:
  - Reads graph from `graph.edge` file (format: `u v weight`).
  - Builds adjacency list (`graph`) and sets `num_vertices`.
  - Defines example dynamic updates (e.g., insert edge `(0→3, 3.0)` when vertex `1` processed, update `(1→2, 1.0)` when vertex `2` processed).
  - Defines post-Dijkstra inserted edges (e.g., `(0→2, 2.0)`, `(1→4, 1.5)`).

- **Run Dynamic Dijkstra (`dijkstra_dynamic`)**:
  - Initializes `distances` (all `INT_MAX` except source = 0), `parent` (all `-1`), and a `set` for vertices (starts with `{0, source}`).
  - Groups updates by trigger vertex (e.g., updates for vertex `1`).
  - While `set` is not empty:
    - Extracts vertex `u` with smallest distance.
    - Applies updates triggered by `u`:
      - Inserts new edge or updates edge weight in `graph`.
      - If the update reduces distance to destination `v`, updates `distances[v]`, `parent[v]`, and inserts `v` into `set`.
    - Processes neighbors of `u`:
      - For each neighbor `v`, if `distances[u] + weight < distances[v]`, updates `distances[v]`, `parent[v]`, and inserts `v` into `set`.
    - Ignores outdated `set` entries to ensure correctness.

- **Print Initial Results**:
  - Outputs `distances` and `parent` for reachable vertices after dynamic Dijkstra, reflecting changes made during execution.

- **Post-Dijkstra Update (`sosp_update`)**:
  - Applies predefined inserted edges (e.g., `(0→2, 2.0)`, `(1→4, 1.5)`).
  - Groups edges by destination vertex.
  - Updates `distances` and `parent` for directly affected vertices.
  - Propagates changes to neighbors iteratively until no further updates needed.

- **Print Final Results**:
  - Outputs updated `distances` and `parent` after post-Dijkstra updates.

  - **Compile and Run**:
    - Compile the program using: `g++ -o main.cpp`
    - Run the program with: `./main`
- **End**:
  - Program terminates, providing shortest paths before and after post-Dijkstra edge insertions.