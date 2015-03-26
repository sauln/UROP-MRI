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

int test_line_parameterization(default_Graph & g){



	//create a deep copy of the graph
	//(eventually make a map to associate the two)
	//this graph we will slowly deform 
	default_Graph g_new(g);






	/* we want to try a simple paramterization*/
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances


	//This code is going to look really ugly - 
	std::vector<V> corners_orig(4);
	std::vector<V> corners_new(4);
	//upper_left, upper_right, lower_left, lower_right;

	//lets define our unit square.

	//std::map<V, V> param_map;
	

	write_point(corners_orig[0], 2.1, 2.1, 1.7);//upper_right
	write_point(corners_orig[1], 2.0, 1.1, 2);//upper_left
	write_point(corners_orig[2], 0.5, 2.4, 2.5); //lower_right
	write_point(corners_orig[3], -0.3, 1.4, 1.5); //lower_left

	write_point(corners_new[0], 1.0, 1.0, 0.0);//upper_right
	write_point(corners_new[1], 1.0, 0.0, 0.0);//upper_left
	write_point(corners_new[2], 0.0, 1.0, 0.0); //lower_right
	write_point(corners_new[3], 0.0, 0.0, 0.0); //lower_left

	//param_map.insert(std::pair<V, V>(corners_orig[1], corners_new[1]));
	//param_map[corners_orig[1]] = corners_new[1];
	//for (int i = 0; i < 4; i++)
	//	param_map[corners_new[i]] =  corners_orig[i];

	//to find the new location of a node:   
	//*******************
	//  Let A,B be two corner points in our original space.
	//  Let a,b be two corner points in our new space
	//  a point C in [A,B] can be represented as a point c in [a,b] 
	//  and c =  ( d(A,C)/d(A,B) * (b-a)/d(a,b) ) + a
	//
	//*******************

	std::vector<vertex_descriptor> vertex_separator;
	find_shortest_path(g, corners_orig[0], corners_orig[1], predecessors, distances, vertex_separator);
	int n = vertex_separator.size();
	for (int i = 0; i < n; i++)
		std::cout << vertex_separator[i] << std::endl;

	//sets our new corners in place
	g_new[vertex_separator[0]] = corners_new[0];
	g_new[vertex_separator[n - 1]] = corners_new[1];

	//g_new[vertex_separator[1]] = 
	std::cout << "vertex_separator[0] is index:" << vertex_separator[0] << std::endl;
	std::cout << "vertex_separator[n-1] is index:" << vertex_separator[n - 1] << std::endl;

	std::cout << "CHECK " << std::endl << std::endl << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << distances[vertex_separator[i]] <<std::endl;
	std::cout << std::endl;




	double scale = distances[vertex_separator[0]]; // = d(A,B)
	V direction = diff_direction(g_new[vertex_separator[0]], g_new[vertex_separator[n - 1]]); // = direction(A,B)
	double dir_scale = diff_V(g_new[vertex_separator[n - 1]], g_new[vertex_separator[0]]);
	V tmp;


	for (int i = 0; i < n-1; i++)
	{
		double portion = distances[vertex_separator[i]]; // = d(A,C)
		tmp = project_ish(scale, portion, direction, dir_scale, g_new[vertex_separator[n - 1]]);
		g_new[vertex_separator[i]] = tmp;
	}
	//std::cout << std::endl;
	//std::cout << "point number " << vertex_separator[n - 1] << " is located at" << std::endl;
	//print_point(g_new[vertex_separator[n - 1]]);
	//std::cout << "point number " << vertex_separator[0] << " is located at" << std::endl;
	//print_point(g_new[vertex_separator[0]]);

	//std::cout << "point number " << vertex_separator[1] << " is located at" << std::endl;
	//print_point(g_new[vertex_separator[1]]);
	
	for (int i = 0; i < n; i++){
		print_point(g_new[vertex_separator[i]]);
	}






	//vertex_descriptor up_r = naive_closest_vertex(g, upper_right);
	//vertex_descriptor up_l = naive_closest_vertex(g, upper_left);
	//vertex_descriptor low_r = naive_closest_vertex(g, lower_right);
	//vertex_descriptor low_l = naive_closest_vertex(g, lower_left);
	//std::cout << "Closest vertices: " << up_r << ", " << up_l << ", " << low_r << ", " << low_l << std::endl;
	//print_point(g[up_r]);print_point(g[up_l]);print_point(g[low_r]);print_point(g[low_l]);

	

	//now we take 2 points, find the shortest path between the two.
	//take this shortest path and map it to the interval [0,1]

	/*
	for (int i = 0; i < vertex_separator.size(); i++){
		std::cout << distances[vertex_separator[i]] << std::endl;
	}


	double max = distances[vertex_separator[0]];
	// we neeed to relocate all of the positions.  

	//typedef adjacency_list < setS, vecS, undirectedS,
	//	V, WeightProperty, int, vecS > test_Graph;

	default_Graph relocated(vertex_separator.size());


	//typedef property_map < default_Graph, vertex_index_t >::type IndexMap;

	for (int i = 0; i < vertex_separator.size(); i++){
		write_vertex(relocated, i, 0, 0, distances[vertex_separator[i]] / max);
		std::cout << "The " << vertex_separator[i] << " vertex is now at " << i << "." << std::endl;
	}
	/* Try to now map to an entire boundary edge.   */





	//for (int i = 0; i < vertex_separator.size() - 1; i++){
	//	add_edge_N(i, i + 1, relocated);
	//}

	//print_all_edges(relocated);
	//print_all_vertices(relocated);

	//print_point(g[vertex_separator[0]]);
	//std::cout << vertex_separator[0] << std::endl;
	//std::cout << distances[vertex_separator[0]] << std::endl;
	
	return 0;
}


int prototype_components(default_Graph& g){

	/* Two points on opposite sides of our test surface */
	V start, end;
	write_point(start, -5, -5, -5);
	write_point(end, 5, 5, 5);

	/*  */

	// Create output storage for Dijkstra + separator
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances
	std::vector<vertex_descriptor> vertex_separator;

	find_shortest_path(g, start, end, predecessors, distances, vertex_separator);

	//create a filtered graph that does not have the vertex_serparator in it.  
	FilteredGraphType filtered_graph = filter_separator(g, vertex_separator);

	std::vector<int> component(num_vertices(filtered_graph));
	int num = find_connected_components(filtered_graph, component);

	print_connected_components(num, component);



	/* Now we want to create separate graphs for each of these */
	/* not actually sure what to do next so we will try to fit the vertex separator onto the straight line of our shape
	then once we know how that works, we can figure out how to through the connected components onto it
	immediately, we will clean the shit out of the code */

	/* Add vertex separator back to connected components as boundary */
	return 0;

}

int main(int argc, char* argv[])
{
	default_Graph g;
	char* filename = "../../surface/surface.off";

	/* Create a graph from the file above */
	g = create_graph(filename);


	test_line_parameterization(g);
	//prototype_components(g);

	int h;
	std::cin >> h;
	return 0;
}

