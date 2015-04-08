// graph_manipulations.cpp : 


#define _CRT_SECURE_NO_WARNINGS

#include "stdafx.h"


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

//#include <boost/graph/dijkstra_shortest_paths.hpp>

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
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_patch_3.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>


#include <CGAL/boost/graph/dijkstra_shortest_paths.h>


#include <CGAL/Polyhedron_items_with_id_3.h>

//#include <CGAL/Two_vertices_parameterizer_3.h>
//#include <CGAL/LSCM_parameterizer_3.h>
//#include <CGAL/OpenNL/linear_solver.h>
//#include "Polyhedron_ex.h"
//#include "Mesh_cutter.h"
//#include "Parameterization_polyhedron_adaptor_ex.h"


#include "graphs.h"
#include "tools.h"


using namespace boost;




// read the mesh from our file 
//int read_mesh(char* filename, Polyhedron & mesh);
//Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron & mesh);



//  This function will take a graph and a set of points
//  It will generate the vertex separator
//int create_vertex_separator(default_Graph &g, V a, V b, std::vector<vertex_descriptor> & vertex_separator);
//int print_mesh(const Polyhedron & mesh, const Parameterization_polyhedron_adaptor & mesh_adaptor);



/*
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
	// Emit EPS trailer. 
	out << "grestore" << std::endl;
	out << std::endl;
	out << "showpage" << std::endl;
	return true;
}
*/

//Read in a mesh from a file.
int read_mesh(char* filename, Polyhedron & mesh){

	std::ifstream stream(filename);

	//***************************************
	// Read the mesh as a polyhedron
	//***************************************

	std::cout << "Read in file" << std::endl;

	stream >> mesh;
	if (!stream || !mesh.is_valid() || mesh.empty())
	{
		std::cerr << "Error: cannot read OFF file " << filename << std::endl;
		return EXIT_FAILURE;
	}
	return 0;
}


/*
//Parameterize the map and return the mesh adaptor
Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron & mesh){
	std::cout << "parameterize mesh and return parameterization" << std::endl;
	// Create Polyhedron adaptor mesh
	// Assume mesh is a topological disk
	Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

	// Parameterize our mesh according to Parameterizer specs
	Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor, Parameterizer());


	// Check to make sure it worked just find -- this is the standard checker
	switch (err) {
	case Parameterizer::OK: // Success
		break;
	case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
	case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
	case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
	case Parameterizer::ERROR_BORDER_TOO_SHORT:
		std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
		//return EXIT_FAILURE;
		break;
	default: // Error
		std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
		//return EXIT_FAILURE;
		break;
	};

	return mesh_adaptor;
}
*/
/*
//Save the mesh as .off and .eps (i think the .off is not printing quite right.
int save_mesh(Parameterization_polyhedron_adaptor & mesh_adaptor)
{
	std::cout << "save mesh " << std::endl;
	//***************************************
	// Save file for visualization
	//***************************************

	Polyhedron new_mesh = mesh_adaptor.get_adapted_mesh();

	// write the polyhedron out as a .OFF file
	// This seems to output the original mesh
	std::ofstream os("dump.off");
	os << new_mesh;
	os.close();

	// Write Postscript file
	// this outputs the new 2d mesh
	if (!write_file_eps(mesh_adaptor, "dump2.eps"))
	{
		std::cerr << "Error: cannot write file " <<  "dump2.eps" << std::endl;
		return EXIT_FAILURE;
	}
	
	return 0; 

}*/
/*
//Print the mesh to the console line
int print_mesh(const Polyhedron & mesh, const Parameterization_polyhedron_adaptor & mesh_adaptor){
	// Raw output: dump (u,v) pairs
	std::cout << "Printing mesh and mesh_adaptor" << std::endl;
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

	return 0;
}
*/

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






//***************************************
//
//  This function does everything we need from start to finish
//
//***************************************

bool isin(vertex_descriptor_mesh vit, std::vector<vertex_descriptor_mesh> vertex_separator){
	return std::find(vertex_separator.begin(), vertex_separator.end(), vit) != vertex_separator.end();
}

int go_deeper(		const Polyhedron & mesh, 
					const std::vector<vertex_descriptor_mesh> & vertex_separator, 
					std::map<vertex_descriptor_mesh, int>& component, 
					int & component_number, vertex_descriptor_mesh vit){


	if (!isin(vit, vertex_separator) && component.find(vit) == component.end())
	{
		component[vit] = component_number;
		//std::cout << "enter node " << vit->id() << ". Set to " << component_number << std::endl;

		edge_it_mesh beg, eeg;
		vertex_descriptor_mesh dog = vit;

		//for each of vertex on the end of an out going edge,
		for (boost::tie(beg, eeg) = out_edges(dog, mesh); beg != eeg; ++beg){
			
			vertex_descriptor_mesh fish = target(*beg, mesh);
			//std::cout << "From " << vit->id() << " into " << fish->id() << std::endl;

			if (!isin(fish, vertex_separator) && component.find(fish) == component.end()){
				//std::cout << fish->id() << "+";
				go_deeper(mesh, vertex_separator, component, component_number, fish);
			}

		}
		//std::cout << std::endl;
		//std::cout << (*vit)->vertex_begin()->next.id() <<std::endl;// << (*vit)->vertex_begin()->next.id();
		return 1;
	}
	else{
		return 0;
	}
}

int label_components(const Polyhedron & mesh, const std::vector<vertex_descriptor_mesh> & vertex_separator, std::map<vertex_descriptor_mesh, int> &component){
	vertex_iterator_mesh vit, ve;

	//for every vertex in the mesh
	int component_number = 0;
	int check = 0;
	for (boost::tie(vit, ve) = vertices(mesh); vit != ve; ++vit)
	{
		//std::cout << "Begin again at " << (*vit)->id() << std::endl;
		check = go_deeper(mesh, vertex_separator, component, component_number, *vit);
		//if the vertex is not in the separator and  has not already been dealt with
		if (check == 1)
			component_number += 1;
	}
	return 0;
}



int dijkstra_with_mesh(){


	//this will be a completely self contained version of dijkstra with mesh

	Point start(-5, -5, -5);
	Point end(5, 5, 5);

	char* filename = "../../surface/surface.off";
	Polyhedron mesh;
	read_mesh(filename, mesh);

	// Associate indices to the vertices
	vertex_iterator_mesh vit, ve;
	int index = 0;
	for (boost::tie(vit, ve) = boost::vertices(mesh); vit != ve; ++vit){
		(*vit)->id() = index++;
	}

	//std::vector<vertex_descriptor_mesh> vertex_separator;
	std::vector<vertex_descriptor_mesh> predecessors(num_vertices(mesh)); // To store parents
	std::vector<Weight> distances(num_vertices(mesh)); // To store distances

	vertex_descriptor_mesh st, en;
	int vert_1 = naive_closest_vertex(mesh, start, st);
	int vert_2 = naive_closest_vertex(mesh, end, en);

	IndexMap_mesh indexMap = get(vertex_index, mesh);
	PredecessorMap_mesh predecessorMap(&predecessors[0], indexMap);
	DistanceMap_mesh distanceMap(&distances[0], indexMap);

	dijkstra_shortest_paths(mesh, st, distance_map(distanceMap).predecessor_map(predecessorMap).distance_combine(dist_combine()).distance_compare(dist_compare()));

	std::cout << "distance of (" << indexMap[st] << ", " << indexMap[en] << ") = " << distanceMap[en] << std::endl;


	//print the shortest path
	vertex_descriptor_mesh it = en;
	do {
		std::cout << (*it).id() << " -> ";
		it = predecessorMap[it];
	} while (it != st);
	std::cout << (*it).id() << std::endl;
	it = en;
	do {
		std::cout << distanceMap[it] << "(" << it->id()  << ")" << " -> ";
		//prev = distanceMap[it];
		it = predecessorMap[it];
	} while (it != st);
	std::cout << distanceMap[it] << std::endl << std::endl<< std::endl;

	//isolate the separator
	std::vector<vertex_descriptor_mesh> vertex_separator;


	it = en;
	while (it != st){ 
		vertex_separator.push_back(it);
		it = predecessorMap[it];
	}
	vertex_separator.push_back(st);

	//print the separator
	std::cout << "Check vertex_separator: " << (*vertex_separator.begin())->id();
	for (std::vector<vertex_descriptor_mesh>::iterator jt = ++vertex_separator.begin(); jt != vertex_separator.end(); ++jt){
		std::cout << " -> " << (*jt)->id();
	}
	std::cout << std::endl;


	//do without the boost filtering
	std::map<vertex_descriptor_mesh, int> component_map;

	

	label_components(mesh, vertex_separator, component_map);

	//set all vertices in v_sep to be their own class
	for (boost::tie(vit, ve) = vertices(mesh); vit != ve; ++vit)		
		if (isin(*vit, vertex_separator))
			component_map[*vit] = -1;

	//print each vertex with its component
	for (boost::tie(vit, ve) = vertices(mesh); vit != ve; ++vit)
		std::cout << (*vit)->id() << "->" << component_map[*vit] << std::endl;





	//now we have a component map, we want to create new polyhedron for each of the components.
	Polyhedron mesh_p1;
	Polyhedron mesh_p2; //one for each connected component

	//somehow copy the mesh into these meshes.



	return 0;
}



int main(int argc, char* argv[])
{

	char* filename = "../../surface/surface.off";

	//compare_differences_between_graph_and_mesh(filename);
	dijkstra_with_mesh();
	//dijkstra_with_graph();


	//take the shortest path found by dijkstra and find the connected components of mesh\separator


	//Parameterization_polyhedron_adaptor new_mesh = parameterize_mesh(mesh);
	//save_mesh(new_mesh);

	int h;
	std::cin >> h;
	return 0;
}

