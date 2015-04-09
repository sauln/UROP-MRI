#ifndef TOOLS
#define TOOLS

#include "graphs.h"


#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <CGAL/boost/graph/dijkstra_shortest_paths.hpp>

//These need to be templated into one function
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


	//std::cout << std::fixed;
	//std::cout << std::setprecision(4);

	return min_ind;
}
int naive_closest_vertex(Polyhedron & g, const Point & p, std::vector<vertex_descriptor_mesh> &vert_descript){
	// cycle through each vertex and save the one closest to the point
	// return the index for that vertex
	//std::cout << "entry point of finding closest vertex, mesh " << std::endl;


	double  min_dist = 1000000;
	int closest_id = -1;
	Vector r;

	//std::cout << "Starting point: " << p << std::endl;
	vertex_descriptor_mesh tmp;
	vertex_iterator_mesh v, ve;
	//vertex_iterator_fin v, ve;
	//boost::tie(v, ve) = vertices(g);

	for (boost::tie(v, ve) = vertices(g); v != ve; ++v){
		//for (vertex_iterator_mesh_2 v = g.vertices_begin(); v != g.vertices_end(); ++v){
		r = (*v)->point() - p;
		//std::cout << "Point: " << (*v)->point() << " id: " << (*v)->id() << " and distance: " << r.squared_length() << std::endl;
		if (r.squared_length() < min_dist){
			tmp = (*v);
			min_dist = r.squared_length();
			closest_id = (*v)->id();
		}

	}
	vert_descript.push_back(tmp);
	//std::cout << "closest index is: " << closest_id << " and that distance is: " << min_dist << std::endl;

	return closest_id;
}
int naive_closest_vertex(Finite_Polyhedron & g, const Point & p){
	// cycle through each vertex and save the one closest to the point
	// return the index for that vertex
	std::cout << "entry point of finding closest vertex, mesh " << std::endl;


	double  min_dist = 1000000;
	int closest_id = -1;
	Vector r;
	
	//std::cout << "Starting point: " << p << std::endl;
	
	vertex_iterator_mesh v, ve;
	//vertex_iterator_fin v, ve;
	//boost::tie(v, ve) = vertices(g);

	for (boost::tie(v, ve) = vertices(g.m_g); v != ve; ++v){
	//for (vertex_iterator_mesh_2 v = g.vertices_begin(); v != g.vertices_end(); ++v){
		r = (*v)->point() - p;
		//std::cout << "Point: " << (*v)->point() << " id: " << (*v)->id() << " and distance: " << r.squared_length() << std::endl;
		if (r.squared_length() < min_dist){
			min_dist = r.squared_length();
			closest_id = (*v)->id();
		}

	}
	
	//std::cout << "closest index is: " << closest_id << " and that distance is: " << min_dist << std::endl;

	return closest_id;
}


Point V2Point(const V & v){
	Point p(v.x, v.y, v.z);
	return p;
}

int find_shortest_path(Finite_Polyhedron & g, Point start, Point end,
	std::vector<vertex_descriptor_mesh> &predecessors,
	std::vector<Weight> &distances,
	std::vector<vertex_descriptor_mesh> &vertex_separator){


	//We need to define a weight map that describes the d(a,b) = (a-b).squared_length()


	//std::cout <<" entry point of find shortest path, mesh. " <<std::endl;

	int vert_1 = naive_closest_vertex(g, start);
	int vert_2 = naive_closest_vertex(g, end);
	//std::cout << "closest vertices are " << vert_1 << " and " << vert_2 << std::endl;


	// Create maps for output for Dijkstra
	IndexMap_mesh indexMap = get(vertex_index, g.m_g);
	PredecessorMap_mesh predecessorMap(&predecessors[0], indexMap);
	DistanceMap_mesh distanceMap(&distances[0], indexMap);
	// Compute shortest paths from vert_1 to all vertices
	//dijkstra_shortest_paths(g.m_g, vert_1, distance_map(distanceMap).predecessor_map(predecessorMap));
	//std::cout << "distance of (" << indexMap[vert_1] << ", " << indexMap[vert_2] << ") = " << distanceMap[vert_2] << std::endl;
	// Create vertex separator - the shortest path from vert_1 to vert_2 */
	//std::vector<vertex_descriptor>::iterator it;
	//int v = vert_2;
	//vertex_separator.push_back(v);
	//while (v != predecessorMap[v]){
	//	v = predecessorMap[v];
	//	vertex_separator.push_back(v);
	//}
	//Print vertex separator 
	//for (it = vertex_separator.begin(); it != vertex_separator.end(); it++){
	//	std::cout << *it << " -> ";
	//}
	//std::cout << std::endl;


	return 0;

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


FilteredMeshType filter_out_separator(Polyhedron & mesh, std::vector<vertex_descriptor_mesh> &vertex_separator){
	/* Create a filtered graph of original graph without the vertex separator */

	std::cout << "Filter out vertex separator" << std::endl;

	//initialize filter -- turning this into a class would make this a little more automatic.
	vertex_id_filter_mesh<Polyhedron> filter;
	filter.vertex_separator = vertex_separator;
	filter.size_of_vs = vertex_separator.size();


	//creating filtered  graph

	FilteredMeshType filtered_graph(mesh, boost::keep_all(), filter); // (graph, EdgePredicate, VertexPredicate)


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





int create_vertex_separator(default_Graph &g, V a, V b, std::vector<vertex_descriptor> & vertex_separator){
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances

	find_shortest_path(g, a, b, predecessors, distances, vertex_separator);
	return 0;
}

int separate_graph(default_Graph & g, V a, V b){


	std::vector<vertex_descriptor> vertex_separator;
	create_vertex_separator(g, a, b, vertex_separator);


	//create a filtered graph that does not have the vertex_serparator in it.
	FilteredGraphType filtered_graph = filter_separator(g, vertex_separator);

	std::vector<int> component(num_vertices(filtered_graph));
	int num = find_connected_components(filtered_graph, component);

	print_connected_components(num, component);


	return 0;
}


int test_line_parameterization(default_Graph & g){



	//*******************************
	//  Takes a linear vertex separator and straightens it out to the unit line
	//  Tried to setup some basic border parameterization.
	//  ended up being very complicated and it's better
	//  to just use the CGAL methods
	//*******************************

	//create a deep copy of the graph
	//(eventually make a map to associate the two)
	//this graph we will slowly deform
	default_Graph g_new(g);

	// we want to try a simple paramterization
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

	std::cout << "vertex_separator[0] is index:" << vertex_separator[0] << std::endl;
	std::cout << "vertex_separator[n-1] is index:" << vertex_separator[n - 1] << std::endl;

	std::cout << "CHECK " << std::endl << std::endl << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << distances[vertex_separator[i]] << std::endl;
	std::cout << std::endl;


	double scale = distances[vertex_separator[0]]; // = d(A,B)
	V direction = diff_direction(g_new[vertex_separator[0]], g_new[vertex_separator[n - 1]]); // = direction(A,B)
	double dir_scale = diff_V(g_new[vertex_separator[n - 1]], g_new[vertex_separator[0]]);
	V tmp;

	//this takes each of the points and places them on a straight line between points a and b.
	for (int i = 0; i < n - 1; i++){
		double portion = distances[vertex_separator[i]]; // = d(A,C)
		tmp = project_ish(scale, portion, direction, dir_scale, g_new[vertex_separator[n - 1]]);
		g_new[vertex_separator[i]] = tmp;
	}

	for (int i = 0; i < n; i++)
		print_point(g_new[vertex_separator[i]]);


	return 0;
}
int dijkstra_with_graph(){
	char* filename = "../../surface/surface.off";
	default_Graph g;
	g = read_graph(filename);


	V start, end;
	write_point(start, -5, -5, -5);
	write_point(end, 5, 5, 5);


	std::vector<vertex_descriptor> vertex_separator;
	std::vector<vertex_descriptor> predecessors(num_vertices(g)); // To store parents
	std::vector<Weight> distances(num_vertices(g)); // To store distances

	int vert_1 = naive_closest_vertex(g, start);
	int vert_2 = naive_closest_vertex(g, end);

	IndexMap indexMap = get(vertex_index, g);
	PredecessorMap predecessorMap(&predecessors[0], indexMap);
	DistanceMap distanceMap(&distances[0], indexMap);


	dijkstra_shortest_paths(g, vert_1, distance_map(distanceMap).predecessor_map(predecessorMap));

	std::cout << "distance of (" << indexMap[vert_1] << ", " << indexMap[vert_2] << ") = " << distanceMap[vert_2] << std::endl;

	int v = vert_2;

	while (v != predecessorMap[v]){
		std::cout << v << " -> ";
		v = predecessorMap[v];
	}
	std::cout << v << std::endl;

	v = predecessorMap[vert_2];
	double prev = distanceMap[vert_2];

	while (v != predecessorMap[v]){
		std::cout << prev - distanceMap[v] << " -> ";
		prev = distanceMap[v];
		v = predecessorMap[v];
	}

	std::cout << distanceMap[v] << std::endl << std::endl;

	return 0;
}


#endif