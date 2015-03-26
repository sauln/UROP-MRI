// adj_list_examples.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include "adjacency_list_io.h"


#include <boost/config.hpp>

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
#error adjacency_list_io.hpp has not been ported to work with VC++
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>

#include <boost/graph/copy.hpp>
#include <boost/array.hpp>
//#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/graph/graph_utility.hpp>


#include "graphs.h"
#include "tools.h"


using namespace boost;




int find_vertex_separator(default_Graph & g, V start, V end, std::vector<vertex_descriptor> &vertex_separator){
	
	//find vertex index for start andend 
	vertex_descriptor vert_1 = naive_closest_vertex(g, start);
	vertex_descriptor vert_2 = naive_closest_vertex(g, end);


	// Create output storage for Dijkstra
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances

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


int main(int argc, char* argv[])
{
	default_Graph g;
	char* filename = "../surface/surface.off";

	/* Create a graph from the file above */
	g = create_graph(filename);
	
	/* Two points on opposite sides of our test surface */
	V start, end;
	write_point(start, -5, -5, -5);
	write_point(end, 5, 5, 5);

	/*  */
	std::vector<vertex_descriptor> vertex_separator;
	find_vertex_separator(g, start, end, vertex_separator);


	FilteredGraphType filtered_graph = filter_separator(g, vertex_separator);

	std::vector<int> component(num_vertices(filtered_graph));
	int num = find_connected_components(filtered_graph, component);

	print_connected_components(num, component);



	/* Now we want to create separate graphs for each of these */
	/* not actually sure what to do next so we will try to fit the vertex separator onto the straight line of our shape
	then once we know how that works, we can figure out how to through the connected components onto it
	immediately, we will clean the shit out of the code */



	/* Add vertex separator back to connected components as boundary */





	int h;
	std::cin >> h;
	return 0;
}

