#ifndef TOOLS
#define TOOLS

#include "graphs.h"


#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

int naive_closest_vertex(default_Graph g, V p);

int naive_closest_test(const default_Graph & g, V p ){
	// this is a really shoddy test that can be expanded. 
	// currently, it was written solely to check that naive closest was returning the correct index. 
	// surprisingly, v(19) is closer to (5,5,5) than v(24).
	std::cout << diff_V(g[19], p) << std::endl;
	std::cout << diff_V(g[24], p) << std::endl;
	return 0;
}
int naive_closest_vertex(default_Graph g, V p){
	// cycle through each vertex and save the one closest to the point
	// return the index for that vertex

	double min_norm = 100;
	double tmp_norm;
	int min_ind = -1;
	for (int i = 0; i < num_vertices(g); i++){
		tmp_norm = diff_V(p, g[i]);
		if (tmp_norm < min_norm){
			min_norm = tmp_norm;
			min_ind = i;
		}
	}


	std::cout << std::fixed;
	std::cout << std::setprecision(4);

	//std::cout << "\nminimum difference from " << p.x << ", " << p.y << ", " << p.z << std::endl
	//	<< "is: " << min_norm
	//	<< "\nFound from node: " << min_ind
	//	<< "\nAt location:" << g[min_ind].x << ", " << g[min_ind].y << ", " << g[min_ind].z << std::endl;

	return min_ind;
}



int find_shortest_path(default_Graph & g, V start, V end, 
						std::vector<vertex_descriptor> &predecessors, 
						std::vector<Weight> &distances, 
						std::vector<vertex_descriptor> &vertex_separator){

	//split this into two smaller functions
	//find vertex index for start andend 
	vertex_descriptor vert_1 = naive_closest_vertex(g, start);
	vertex_descriptor vert_2 = naive_closest_vertex(g, end);



	// Create maps for output for Dijkstra
	IndexMap indexMap = get(vertex_index, g);
	PredecessorMap predecessorMap(&predecessors[0], indexMap);
	DistanceMap distanceMap(&distances[0], indexMap);

	// Compute shortest paths from vert_1 to all vertices
	dijkstra_shortest_paths(g, vert_1, distance_map(distanceMap).predecessor_map(predecessorMap));
	//std::cout << "distance of (" << indexMap[vert_1] << ", " << indexMap[vert_2] << ") = " << distanceMap[vert_2] << std::endl;
	
	/* Create vertex separator - the shortest path from vert_1 to vert_2 */
	std::vector<vertex_descriptor>::iterator it;
	int v = vert_2;
	vertex_separator.push_back(v);
	while (v != predecessorMap[v]){
		v = predecessorMap[v];
		vertex_separator.push_back(v);
	}

	/* Print vertex separator */
	//for (it = vertex_separator.begin(); it != vertex_separator.end(); it++){
	//	std::cout << *it << " -> ";
	//}
	//std::cout << std::endl;


	return 0;

}


FilteredGraphType filter_separator(default_Graph & g, std::vector<vertex_descriptor> &vertex_separator){
	/* Create a filtered graph of original graph without the vertex separator */

	std::cout << "Filter out vertex separator" << std::endl;

	//initialize filter -- turning this into a class would make this a little more automatic.
	vertex_id_filter<default_Graph> filter;
	filter.vertex_separator = vertex_separator;
	filter.size_of_vs = vertex_separator.size();


	//creating filtered  graph

	FilteredGraphType filtered_graph(g, boost::keep_all(), filter); // (graph, EdgePredicate, VertexPredicate)


	print_graph(filtered_graph);
	return filtered_graph;
}

int find_connected_components(FilteredGraphType filtered_graph, std::vector<int> & component)
{
	/* Find all connected components of subgraph w/0 vertex separator */

	int num = connected_components(filtered_graph, &component[0]);

	return num;
}
int print_connected_components(int num, std::vector<int> & component)
{
	std::cout << "There are " << num << " connected components" << std::endl;

	std::vector<int>::size_type i;
	std::cout << "Total number of components: " << num << std::endl;
	for (i = 0; i != component.size(); ++i)
		std::cout << "Vertex " << i << " is in component " << component[i] << std::endl;
	std::cout << std::endl;


	return 0;
}



#endif