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


#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>



#include "graphs.h"
#include "tools.h"


using namespace boost;




int dijkstra_test(default_Graph & g, V start, V end)
{
	int vert_1 = naive_closest_vertex(g, start);
	int vert_2 = naive_closest_vertex(g, end);


	// Create things for Dijkstra
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances

	IndexMap indexMap = get(vertex_index, g);
	PredecessorMap predecessorMap(&predecessors[0], indexMap);
	DistanceMap distanceMap(&distances[0], indexMap);

	// Compute shortest paths from v0 to all vertices, and store the output in predecessors and distances
	// boost::dijkstra_shortest_paths(g, v0, boost::predecessor_map(predecessorMap).distance_map(distanceMap));
	// This is exactly the same as the above line - it is the idea of "named parameters" - you can pass the
	// prdecessor map and the distance map in any order.
	dijkstra_shortest_paths(g, vert_1, distance_map(distanceMap).predecessor_map(predecessorMap));



	// Output results
	std::cout << "distances and parents:" << std::endl;
	//NameMap nameMap = boost::get(boost::vertex_name, g);

	BGL_FORALL_VERTICES(v, g, default_Graph)
	{
		std::cout << "distance(" << indexMap[vert_1] << ", " << indexMap[v] << ") = " << distanceMap[v] << ", ";
		std::cout << "predecessor(" << indexMap[v] << ") = " << indexMap[predecessorMap[v]] << std::endl;
	}

	std::cout << "distance of (" << indexMap[vert_1] << ", " << indexMap[vert_2] << ") = " << distanceMap[vert_2] << std::endl;
	return 0;





}






/*int demo_Graph()
{

	
	E edge_array[] = { E(0, 1), E(1, 2), E(0, 2), E(1,3), E(0,3) };
	const std::size_t m = sizeof(edge_array) / sizeof(E);

	default_Graph g(edge_array, edge_array+m, 10);


	g[0].x = 0;
	g[1].x = 2.5;
	std::cout << "test:";
	std::cout << g[1].x;
	write_vertex(g, 2, 0.1, 0.2, 0.3);
	read_vertex(g, 2);
	std::cout << "number of vertices: " << num_vertices(g) << std::endl;
	std::cout << "number of edges: " << num_edges(g) << std::endl;

	//jeeze BGL is a pain in the ass.  
	//this prints out all of the edges, twice - not that big of deal though, right?

	graph_traits<default_Graph>::out_edge_iterator e, e_end;
	graph_traits<default_Graph>::vertex_descriptor s;
	
	for (unsigned int i = 0; i < num_vertices(g); i++){
		s = vertex(i, g);
		for (boost::tie(e, e_end) = out_edges(s, g); e != e_end; ++e)
			std::cout << "(" << source(*e, g)
			<< "," << target(*e, g) << ")" << std::endl;
	} 

	//edges(g)








	/*
	const std::size_t n = 3;
	E edge_array[] = { E(0, 1), E(0, 2), E(0, 1) };
	
	V vert_array[] = { V(0.0, 0.0, 0.0), V(1.0, 0.0, 1.0), V(1.0, 0.0, 0.0) };
	const std::size_t m = sizeof(edge_array) / sizeof(E);

	//my_Graph g1;
	Graph g(edge_array, edge_array + m, n);
	//std::cout << edges(g).first[1] << std::endl;
	//for (std::size_t i = 0; i < m; ++i)
	//	std::cout << edges(g).first[i] << " ";
	std::cout << "graph g1 from file " <<  "\n"
		<< write(g, property<int, edge_weight_t>())
		<< std::endl;
	
	return 0;
}*/





int main(int argc, char* argv[])
{
	default_Graph g;
	char* filename = "../surface/surface.off";

	g = create_graph(filename);
	

	V start, end;
	write_point(start, -5, -5, -5);
	write_point(end, 5, 5, 5);
	dijkstra_test(g, start, end);




	int h;
	std::cin >> h;
	return 0;
}

