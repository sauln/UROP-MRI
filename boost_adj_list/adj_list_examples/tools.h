#ifndef TOOLS
#define TOOLS

#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>




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

	std::cout << "\nminimum difference from " << p.x << ", " << p.y << ", " << p.z << std::endl
		<< "is: " << min_norm
		<< "\nFound from node: " << min_ind
		<< "\nAt location:" << g[min_ind].x << ", " << g[min_ind].y << ", " << g[min_ind].z << std::endl;

	return min_ind;
}



#endif