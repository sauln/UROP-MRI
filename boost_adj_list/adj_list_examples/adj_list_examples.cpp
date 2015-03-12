// adj_list_examples.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include "adjacency_list_io.h"


#include <boost/config.hpp>

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
#error adjacency_list_io.hpp has not been ported to work with VC++
#endif

#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

//#include "adjacency_list.h"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
using namespace boost;


struct V {double x;double y;double z;};
typedef std::pair<std::size_t, std::size_t> E;
struct Weight { double weight; };


typedef adjacency_list < setS, vecS, undirectedS,
	V, property < edge_weight_t, double >, no_property, vecS > default_Graph;
typedef graph_traits<default_Graph>::edge_iterator EdgeIterator;
typedef std::pair<EdgeIterator, EdgeIterator> EdgePair;


int write_point(V &p, double x, double y, double z){
	p.x = x; p.y = y; p.z = z;
	return 0;
}
int write_vertex(default_Graph &g, int position, double x, double y, double z){
	write_point(g[position], x, y, z);
	return 0;
}

double diff_V(V p, V q){
	double sum = sqrt(pow((p.x - q.x), 2) + pow((p.y - q.y), 2) + pow((p.z - q.z), 2));
	return sum;

}

int add_edge_N(int in, int out, default_Graph &g){
	double w = diff_V(g[in], g[out]);

	if (in <= out)
		add_edge(in, out, w, g);
	else
		add_edge(out, in, w, g);
	return 0;

}

int read_vertex(default_Graph g, int position){
	std::cout << "Vertex is located at: (" << g[position].x << "," << g[position].y << "," <<g[position].z << ")\n";
	return 0;
}


int print_all_edges(default_Graph g){

	graph_traits<default_Graph>::out_edge_iterator e, e_end;
	graph_traits<default_Graph>::vertex_descriptor s;


	for (unsigned int i = 0; i < num_vertices(g); i++){
		s = vertex(i, g);
		for (boost::tie(e, e_end) = out_edges(s, g); e != e_end; ++e)
			std::cout << "(" << source(*e, g)
			<< "," << target(*e, g) << ")" << std::endl;
	}

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

	std::cout << "\nminimum difference from " << p.x << ", " << p.y << ", " << p.z << std::endl
		<< "is: " << min_norm
		<< "\nFound from node: " << min_ind
		<< "\nAt location:" << g[min_ind].x << ", " << g[min_ind].y << ", " << g[min_ind].z << std::endl;

	return min_ind;
}

int demo_graph_io()
{
	/*  This function will be a test for reading in the edges from a file and creating a graph out of it   */
	int num_verts;

	char* filename = "../surface/surface.off";
	for (int i = 0; i < 9; i++)
	{
		std::cout << filename[i];
	}
	std::cout << "filename is " << filename << "\n";

	std::ifstream infile(filename);
	std::string line; getline(infile, line);


	double a, b, c, d;
	infile >> num_verts >>b >> c;
	std::cout << "num vertices: " << num_verts << " and triangles: " << b << std::endl;
	
	default_Graph g(num_verts);

	//read and save vertex properties
	for (int i = 0; i < num_verts; i++){
		infile >> a >> b >> c;
		write_vertex(g, i, a, b, c);
	}

	int face_check = 0;
	//read and save each triangle as 3 edges
	int q, w, e, r;

	while (infile >> q >> w >> e >> r){
		face_check++;
		add_edge_N(w, e, g);
		add_edge_N(e, r, g);
		add_edge_N(r, w, g);
	}

	std::cout << "number of edges: " << num_edges(g) << "\nNumber of vertices: " << num_vertices(g) << std::endl;
	
	// now we want to choose a point and find the closest point to a landmark
	V landmark_1, landmark_2;
	write_point(landmark_1, -5, -5, -5);
	write_point(landmark_2, 5, 5, 5);
	int vert_1 = naive_closest_vertex(g, landmark_1);
	int vert_2 = naive_closest_vertex(g, landmark_2);

 
	//now find the shortest path between these points: use Djikstra
	//dijkstra_shortest_path
	EdgePair ep;

	Weight* weight = nullptr;
	//for (ep = edges(g); ep.first != ep.second; ++ep.first)
	//{
		//weight = get(ep.first, g);
		//std::cout << " \n All of the edge weights \n";
		//std::cout << weight->weight << std::endl;
		// Get the two vertices that are joined by this edge...
	//}

	return 0;
}


int test_add_edge(){
	//it seems that the add edge for the undirected graph is still adding multiple edges..
	default_Graph g(10);

	std::cout << num_vertices(g);

	add_edge(0, 1, g);

	std::cout << num_edges(g);

	add_edge(0, 1, g);
	std::cout << num_edges(g);
	return 0;
}




int demo_Graph()
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
	*/
	return 0;
}





int main(int argc, char* argv[])
{
	std::cout << "There were " << argc << " parameters!!\n";
	// argv is an array containing same;
	for (int i = 0; i<argc; i++)
		std::cout << "Parameter " << i << " was " << argv[i] << "\n";

	//io_Graph(argv[1]);
	//demo_Graph();
	demo_graph_io();
	//test_add_edge();

	int h;
	std::cin >> h;
	return 0;
}

