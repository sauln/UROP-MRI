// graph_manipulations.cpp : Defines the entry point for the console application.
//


#define _CRT_SECURE_NO_WARNINGS

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



//*******************
// BOOST libraries --- some of these are not used anymore
//*******************
#include <boost/graph/copy.hpp>
#include <boost/array.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>


//************************
// CGAL Libraries
//************************
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>



#include "graphs.h"
#include "tools.h"


using namespace boost;


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;


typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
Parameterization_polyhedron_adaptor;


static bool write_file_eps(const Parameterization_polyhedron_adaptor& mesh_adaptor,
	const char *pFilename,
	double scale = 500.0)
{
	const Polyhedron& mesh = mesh_adaptor.get_adapted_mesh();
	std::ofstream out(pFilename);
	if (!out)
		return false;
	CGAL::set_ascii_mode(out);
	// compute bounding box
	double xmin, xmax, ymin, ymax;
	xmin = ymin = xmax = ymax = 0;
	Polyhedron::Halfedge_const_iterator pHalfedge;
	for (pHalfedge = mesh.halfedges_begin();
		pHalfedge != mesh.halfedges_end();
		pHalfedge++)
	{
		double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
		double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
		double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
		double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
		xmin = (std::min)(xmin, x1);
		xmin = (std::min)(xmin, x2);
		xmax = (std::max)(xmax, x1);
		xmax = (std::max)(xmax, x2);
		ymax = (std::max)(ymax, y1);
		ymax = (std::max)(ymax, y2);
		ymin = (std::min)(ymin, y1);
		ymin = (std::min)(ymin, y2);
	}
	out << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
	out << "%%BoundingBox: " << int(xmin + 0.5) << " "
		<< int(ymin + 0.5) << " "
		<< int(xmax + 0.5) << " "
		<< int(ymax + 0.5) << std::endl;
	out << "%%HiResBoundingBox: " << xmin << " "
		<< ymin << " "
		<< xmax << " "
		<< ymax << std::endl;
	out << "%%EndComments" << std::endl;
	out << "gsave" << std::endl;
	out << "0.1 setlinewidth" << std::endl;
	// color macros
	out << std::endl;
	out << "% RGB color command - r g b C" << std::endl;
	out << "/C { setrgbcolor } bind def" << std::endl;
	out << "/white { 1 1 1 C } bind def" << std::endl;
	out << "/black { 0 0 0 C } bind def" << std::endl;
	// edge macro -> E
	out << std::endl;
	out << "% Black stroke - x1 y1 x2 y2 E" << std::endl;
	out << "/E {moveto lineto stroke} bind def" << std::endl;
	out << "black" << std::endl << std::endl;
	// for each halfedge
	for (pHalfedge = mesh.halfedges_begin();
		pHalfedge != mesh.halfedges_end();
		pHalfedge++)
	{
		double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
		double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
		double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
		double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
		out << x1 << " " << y1 << " " << x2 << " " << y2 << " E" << std::endl;
	}
	/* Emit EPS trailer. */
	out << "grestore" << std::endl;
	out << std::endl;
	out << "showpage" << std::endl;
	return true;
}











int figure_out_how_to_use_CGAL(char* filename){
	//try to do almost exactly what I've done with boost//
	//read in the surface, but as a CGAL structure.  Then let's try to use
	//the shortest path just like we did with the boost graph.
	//then we will try to use CGAL functions to parameterize.  


	std::ifstream stream(filename);
	//***************************************
	// Read the mesh
	//***************************************

	std::cout << "Read in file" << std::endl;
	Polyhedron mesh;
	stream >> mesh;
	if (!stream || !mesh.is_valid() || mesh.empty())
	{
		std::cerr << "Error: cannot read OFF file " << filename << std::endl;
		return EXIT_FAILURE;
	}


	//***************************************
	// Create Polyhedron adaptor
	// Note: no cutting => we support only
	// meshes that are topological disks
	//***************************************


	Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

	//***************************************
	// Floater Mean Value Coordinates parameterization
	// (defaults are circular border and OpenNL solver)
	//***************************************







	typedef CGAL::Parameterizer_traits_3<Parameterization_polyhedron_adaptor>
		Parameterizer;  // Type that defines the error codes

	Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor);
	switch (err) {
	case Parameterizer::OK: // Success
		break;
	case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
	case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
	case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
	case Parameterizer::ERROR_BORDER_TOO_SHORT:
		std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
		return EXIT_FAILURE;
		break;
	default: // Error
		std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
		return EXIT_FAILURE;
		break;
	};


	//***************************************
	// Output
	//***************************************

	// Raw output: dump (u,v) pairs




	Polyhedron::Vertex_const_iterator pVertex;
	for (pVertex = mesh.vertices_begin();
		pVertex != mesh.vertices_end();
		pVertex++)
	{
		// (u,v) pair is stored in any halfedge
		double u = mesh_adaptor.info(pVertex->halfedge())->uv().x();
		double v = mesh_adaptor.info(pVertex->halfedge())->uv().y();
		//double z = mesh_adaptor.info(pVertex->halfedge())->uv().z();
		std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
	}


	//***************************************
	// Save file for visualization
	//***************************************

	Polyhedron new_mesh = mesh_adaptor.get_adapted_mesh();

	// write the polyhedron out as a .OFF file
	std::ofstream os("dump.off");
	os << new_mesh;
	os.close();


	//***************************************
	// Output
	//***************************************
	// Write Postscript file
	if (!write_file_eps(mesh_adaptor, "dump2.eps"))
	{
		std::cerr << "Error: cannot write file " <<  "dump2.eps" << std::endl;
		return EXIT_FAILURE;
	}




	return 0;
}


/*
// barebones .OFF file reader, throws away texture coordinates, normals, etc.
// stores results in input coords array, packed [x0,y0,z0,x1,y1,z1,...] and
// tris array packed [T0a,T0b,T0c,T1a,T1b,T1c,...]
//http://jamesgregson.blogspot.com/2012/05/example-code-for-building.html
void load_obj(const char *filename, std::vector<double> &coords, std::vector<int> &tris){
	double x, y, z;
	char line[1024], v0[1024], v1[1024], v2[1024];

	// open the file, return if open fails
	FILE *fp = fopen(filename, "r");
	if (!fp) return;

	// read lines from the file, if the first character of the
	// line is 'v', we are reading a vertex, otherwise, if the
	// first character is a 'f' we are reading a facet
	while (fgets(line, 1024, fp)){
		if (line[0] == 'v'){
			sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
			coords.push_back(x);
			coords.push_back(y);
			coords.push_back(z);
		}
		else if (line[0] == 'f'){
			sscanf(line, "%*s%s%s%s", v0, v1, v2);
			tris.push_back(get_first_integer(v0) - 1);
			tris.push_back(get_first_integer(v1) - 1);
			tris.push_back(get_first_integer(v2) - 1);
		}
	}
	fclose(fp);
}
*/
/*int loadObject() {
	// two vectors to hold point coordinates and
	// triangle vertex indices
	//std::vector<double> coords;
	std::vector<int>    tris;

	// load the input file
	load_obj("input.obj", coords, tris);
	if (coords.size() == 0)
		return 1;

	// build a polyhedron from the loaded arrays
	Polyhedron P;
	polyhedron_builder<HalfedgeDS> builder(coords, tris);
	P.delegate(builder);

	// write the polyhedron out as a .OFF file
	std::ofstream os("dump.off");
	os << P;
	os.close();
}*/

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

	//this takes each of the points and places them on a straight line between points a and b.
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

	figure_out_how_to_use_CGAL(filename);
	//test_line_parameterization(g);
	//prototype_components(g);

	int h;
	std::cin >> h;
	return 0;
}

