// graph_manipulations.cpp : Defines the entry point for the console application.
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




int main(int argc, char* argv[])
{
	default_Graph g;
	char* filename = "../../surface/surface.off";

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

