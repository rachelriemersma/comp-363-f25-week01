class MST:
    def create_edgeless_copy(self, G: list[list[float]]) -> list[list[float]]:
        """
        Create an edgeless copy of the input graph G.
        The result has the same number of vertices but no edges.
        """
        # Get the no-edge sentinel value from the diagonal
        no_edge = G[0][0]
        
        # Get the number of vertices
        n = len(G)
        
        # Create a new matrix filled with no-edge values
        T = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append(no_edge)
            T.append(row)
        
        return T

    def reachability_of(self, s: int, G: list[list[float]]) -> list[int]:
        """
        Compute the vertices reachable from a given starting vertex
        in a directed graph represented by an adjacency matrix.
        """
        # Grab what the input graph uses to indicate absence of an edge.
        no_edge = G[0][0]
        
        # Initialize a list to return
        reach = []
        
        # List to remember which vertices to try next. Start from the
        # given vertex.
        visit_next = [s]

        # Explore vertices that we need to visit next. The algorithm stops
        # when the list of vertices to visit next is empty.
        while visit_next:
            # Grab some vertex that we wish to visit next
            v = visit_next.pop()
            
            # Check to see if the vertex we just grabbed is already
            # in the list of vertices reachable from s.
            if v not in reach:
                # The vertex we just grabbed is not in the list of
                # vertices accessible from s so let's add it.
                reach.append(v)
                
                # Now, find all of the neighbors of the vertex we just
                # added to the reach list.
                for u in range(len(G)):
                    if G[v][u] != no_edge:
                        visit_next.append(u)
        
        return reach

    def find_components(self, G: list[list[float]]) -> list[int]:
        """
        Find and label the connected components of a graph.
        """
        n = len(G)
        component = [-1] * n  # -1 means unassigned
        component_id = 0
        
        for vertex in range(n):
            # If this vertex hasn't been assigned to a component yet
            if component[vertex] == -1:
                # Find all vertices reachable from this vertex
                reachable = self.reachability_of(vertex, G)
                
                # Assign all reachable vertices to the same component
                for v in reachable:
                    component[v] = component_id
                
                # Move to next component number
                component_id += 1
        
        return component

    def count_components(self, G: list[list[float]]) -> int:
        """
        Count the number of connected components in a graph.
        """
        component_labels = self.find_components(G)
        return max(component_labels) + 1

    def find_safe_edge(self, G: list[list[float]], components: list[int]) -> tuple[int, int, float]:
        """
        Find the minimum weight edge between any two different components.
        """
        no_edge = G[0][0]
        min_weight = float('inf')
        best_edge = (-1, -1, float('inf'))
        
        n = len(G)
        
        # Look at every possible edge in the original graph
        for i in range(n):
            for j in range(n):
                # If there's an edge and vertices are in different components
                if G[i][j] != no_edge and components[i] != components[j]:
                    if G[i][j] < min_weight:
                        min_weight = G[i][j]
                        best_edge = (i, j, G[i][j])
        
        return best_edge

    def add_edge_to_graph(self, T: list[list[float]], i: int, j: int, weight: float) -> None:
        """
        Add an edge to the graph T (modify in place).
        Since our graph is undirected, we add the edge in both directions.
        """
        T[i][j] = weight
        T[j][i] = weight  # undirected graph

    def min_span_tree(self, G: list[list[float]]) -> list[list[float]]:
        """
        Find the minimum spanning tree (or forest) of a graph using Borůvka's algorithm.
        """
        # Step 1: Initialize T as an edgeless copy of G
        T = self.create_edgeless_copy(G)
        
        # Step 2: Continue until we can't add more edges
        continue_algorithm = True
        
        while continue_algorithm:
            # Find current components in our MST candidate
            components = self.find_components(T)
            num_components = self.count_components(T)
            
            print(f"Current number of components: {num_components}")
            
            # Check if we should continue
            if num_components <= 1:
                print("MST complete - only 1 component remaining")
                continue_algorithm = False
            else:
                # Find the safe edge between components
                safe_edge = self.find_safe_edge(G, components)
                
                # If no safe edge exists, we can't connect more components
                if safe_edge[0] == -1:
                    print("No more safe edges found - graph is disconnected")
                    continue_algorithm = False
                else:
                    # Add the safe edge to our MST
                    print(f"Adding edge: {safe_edge[0]} -> {safe_edge[1]} (weight {safe_edge[2]})")
                    self.add_edge_to_graph(T, safe_edge[0], safe_edge[1], safe_edge[2])
        
        return T


if __name__ == "__main__":
    _ = float('inf')
    G = [
        [_, 1, _, 1, _, _],  # vertex 0's neighbors
        [1, _, _, 1, _, 1],  # vertex 1's neighbors
        [_, _, _, _, _, _],  # vertex 2's neighbors
        [1, 1, _, _, _, 1],  # vertex 3's neighbors
        [_, _, _, _, _, 1],  # vertex 4's neighbors
        [_, 1, _, 1, 1, _]  # vertex 5's neighbors
    ]

    # Create an instance of the MST class
    mst = MST()
    
    # Test the complete Borůvka algorithm
    print("TESTING COMPLETE BORŮVKA ALGORITHM")
    print("="*50)
    
    final_mst = mst.min_span_tree(G)
    print("\nFinal MST result:")
    for row in final_mst:
        print(row)
    
    # Verify the result
    final_components = mst.find_components(final_mst)
    final_num_components = mst.count_components(final_mst)
    print(f"\nFinal number of components: {final_num_components}")