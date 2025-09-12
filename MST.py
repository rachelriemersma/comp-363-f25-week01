class MST:
    """
    Minimum spanning tree implementation using boruvkas algorithm.
    """
    def min_span_tree(self, G: list[list[int]]) -> list[list[int]]:
        """
        Finds the minimum spanning tree of a given graph.
        
        Args:
            G: Adjacency matrix represenatation of a weighted graph(undirected).

        Returns:
            Adjacency matrix of the minimum spanning tree.  
        """
        return self.boruvka(G)
    
    def reachability_of(self, s: int, G: list[list[int]]) -> list[int]:
        """
        Finds all verticies reachable from a given vertex using DFS.

        Args:
            s: Starting vertex index.
            G: Adjacency matrix of the graph.

        Returns:
            List of all verticies reachable from vertex s(inclusive).    
        """
        # no edge placeholder from diagonal
        no_edge: int = G[0][0]
        # store verticies and 'stack' for which verticies to go to next 
        reach: list[int] = []
        visit_next: list[int] = [s]
        # explore until there are no more verticies to visit 
        while visit_next:
            v: int = visit_next.pop()
            if v not in reach:
                reach.append(v)
                # add neighbors of v to 'stack'
                for u in range(len(G)):
                    if G[v][u] != no_edge:
                        visit_next.append(u)
        return reach

    def count_and_label(self, graph: list[list[int]]) -> tuple[int, list[int]]:
        """
        Counts the connected components in a graph and labels each vertex with its component number.

        Args:
            graph: Adjacency matrix of the graph.

        Returns:
            Tuple with:
                count: Number of connected components 
                comp_labels: List where comp_labels is the component number for vertex i.    
        """
        n: int = len(graph)
        visited: list[int] = []
        count: int = 0
        # component label for each vertex 
        comp_labels: list[int] = [0] * n
        # examine every vertex in graph 
        for u in range(n):
            if u not in visited:
                # new component found 
                count += 1
                # find all verticiees reachable from u 
                reachable_from_u: list[int] = self.reachability_of(u, graph)
                # label all verticies in this component with the component number 
                for vertex in reachable_from_u:
                    comp_labels[vertex] = count
                visited.extend(reachable_from_u)    
        return count, comp_labels
    
    def count_components(self, component: list[int]) -> int:
        """
        Counts the number of components in a component array.

        Args:
            component: Array where component[i] is the component ID of vertex i

        Returns:
            Number of distinct component IDs    
        """
        # convert to set to get unique values then count them 
        return len(set(component))

    def create_edgeless_copy(self, G: list[list[int]]) -> list[list[int]]:
        """
        Create and adjacency matrix with the same dimensions as G but no edges.

        Args:
            G: Original adj matrix

        Returns:
            New adj matrix with the same dimensions but no edges     
        """
        n: int = len(G)
        no_edge: int = G[0][0]
        # create empty matrix 
        T: list[list[int]] = []
        for i in range(n):
            row: list[int] = []
            # fill with the no edge placeholder 
            for j in range(n):
                row.append(no_edge)
            T.append(row)
        return T  

    def find_safe_edges(self, G: list[list[int]], component: list[int]) -> list[tuple[int, int, int]]:
        """
        Finds the safe edges: the minimum weight edge that connects a component to another component.
        Finds one safe edge per component and removes duplicates.

        Args:
            G: Original graph adjacency matrix
            component: Current component assignment for each vertex 

        Returns:
            List of tuples (u, v, weight) representing uniqie safe edges     
        """
        n: int = len(G)
        no_edge: int = G[0][0]
        # dictionary to track the safest edge for each component (vertex_u, vertex_v, weight)
        safe_for_component: dict[int, tuple[int, int, int]] = {}
        # examine all possible edges in graph 
        for i in range(n):
            for j in range(i + 1, n):
                if G[i][j] != no_edge: # edge exists 
                    # find which component each vertex belongs to 
                    comp_i: int = component[i]
                    comp_j: int = component[j]
                    # only consider edges between different components 
                    if comp_i != comp_j:
                        # get the weight of that edge 
                        weight: int = G[i][j]
                        # check if this is safest edge for component i
                        if comp_i not in safe_for_component or weight < safe_for_component[comp_i][2]:
                            # update the safe edge for component i 
                            safe_for_component[comp_i] = (i, j, weight)
                        # check if this is the safest edge for component j  
                        if comp_j not in safe_for_component or weight < safe_for_component[comp_j][2]:
                            # update the safe edge for component j
                            safe_for_component[comp_j] = (i, j, weight) 

        # remove duplicate edges      
        # list of edges to return - tracks which edges have already been added             
        unique_edges: list[tuple[int, int, int]] = []
        # process each edge 
        for edge in safe_for_component.values():
            u, v, weight = edge    
            # check if edge is new
            if edge not in unique_edges:
                # add edge to final result 
                unique_edges.append(edge)
        return unique_edges
    
    def merge_components(self, component: list[int], u: int, v: int) -> None:
        """
        Merge two components by updating the component array. 
        Changes all verticies in one component to have the same component ID as the other.

        Args:
            component: Component array to modify 
            u: Vertex in first component
            v: Vertex in second component
        """
        comp_u: int = component[u]
        comp_v: int = component[v]
        # only merge if in different components 
        if comp_u != comp_v:
            # change all verticies with comp_v to have comp_u - merge component v into component u 
            for i in range(len(component)):
                if component[i] == comp_v:
                    component[i] = comp_u

    def boruvka(self, G: list[list[int]]) -> list[list[int]]:
        """
        Boruvkas mst algorithm:
            1. Start with each vertex as its own component
            2. For each component find its safe edge 
            3. Add all safeedges simultaneously 
            4. Repeat until only one component is left or no safe edges exist (either connected graph or disconnected graph)

            Args:
                G: Input graph as adjacency matric where G[i][j] is edge weight ot no edge 

            Returns:
                Minimum spanninh tree as adjacency matrix     
        """
        n: int = len(G)
        # each vertex starts as its own component
        component: list[int] = list(range(n))
        # create mst structure starting with no edges (same size as input graph)
        T: list[list[int]] = self.create_edgeless_copy(G)
        # finishing condition: there are no safe edges (disconnected graph)
        found_safe_edges: bool = True 
        # executes while there are multiple components and a safe edge was found 
        while self.count_components(component) > 1 and found_safe_edges:
            # find safe edges for current components 
            safe_edges: list[tuple[int, int, int]] = self.find_safe_edges(G, component)
            # check if safe edges were found
            found_safe_edges = len(safe_edges) > 0
            if found_safe_edges:
                # add all safe edges to the mst and merge components 
                for u, v, weight in safe_edges:
                    # make sure verticies are still in different components (multiple safe edges could connect same pair of components)
                    if component[u] != component[v]:
                        # add edge to mst 
                        T[u][v] = weight
                        T[v][u] = weight 
                        # merge the components connected by the edge 
                        self.merge_components(component, u, v)
        return T
    
def test():
    mst = MST()
    G = [
        [float('inf'), 2, 3],
        [2, float('inf'), 1], 
        [3, 1, float('inf')]
    ]
    print("Before:")
    for row in G:
        print(row) 
    T = mst.min_span_tree(G)
    print("\nAfter:")
    for row in T:
        print(row)

test()
 
    


