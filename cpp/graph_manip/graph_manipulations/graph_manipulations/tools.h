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
int naive_closest_vertex(Polyhedron & g, const Point & p, vertex_descriptor_mesh & st){
	// cycle through each vertex and save the one closest to the point
	// return the index for that vertex
	//std::cout << "entry point of finding closest vertex, mesh " << std::endl;


	double  min_dist = 1000000;
	int closest_id = -1;
	Vector r;

	//std::cout << "Starting point: " << p << std::endl;

	vertex_iterator_mesh v, ve;
	//vertex_iterator_fin v, ve;
	//boost::tie(v, ve) = vertices(g);

	for (boost::tie(v, ve) = vertices(g); v != ve; ++v){
		//for (vertex_iterator_mesh_2 v = g.vertices_begin(); v != g.vertices_end(); ++v){
		r = (*v)->point() - p;
		//std::cout << "Point: " << (*v)->point() << " id: " << (*v)->id() << " and distance: " << r.squared_length() << std::endl;
		if (r.squared_length() < min_dist){
			st = (*v);
			min_dist = r.squared_length();
			closest_id = (*v)->id();
		}

	}

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



#endif