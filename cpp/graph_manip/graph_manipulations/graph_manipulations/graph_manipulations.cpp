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
#include <stdio.h>


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
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>

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

#include <CGAL/Polyhedron_incremental_builder_3.h>

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

//typedef std::pair<int, int> Ed;


// read the mesh from our file 
int read_mesh(char* filename, Polyhedron & mesh);
Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron & mesh);
int print_mesh(const Polyhedron & mesh, const Parameterization_polyhedron_adaptor & mesh_adaptor);


//  This function will take a graph and a set of points
//  It will generate the vertex separator
//int create_vertex_separator(default_Graph &g, V a, V b, std::vector<vertex_descriptor> & vertex_separator);




Polyhedron partial_mesh_builder(Polyhedron &mesh, std::map<vertex_descriptor_mesh, int> &component_map, int k);

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

	//associate indices to each vertex

	vertex_iterator_mesh vit, ve;
	int index = 0;
	for (boost::tie(vit, ve) = boost::vertices(mesh); vit != ve; ++vit){
		(*vit)->id() = index++;
	}

	return 0;
}

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


//Save the mesh as .off and .eps (i think the .off is not printing quite right.
int save_mesh(Parameterization_polyhedron_adaptor & mesh_adaptor, char * filename)
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
	if (!write_file_eps(mesh_adaptor, filename))
	{
		std::cerr << "Error: cannot write file " <<  "dump2.eps" << std::endl;
		return EXIT_FAILURE;
	}
	
	return 0; 

}

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







//***************************************
//
//  This function does everything we need from start to finish
//
//***************************************
bool isin(std::set<int> keys, int item){
	return std::find(keys.begin(), keys.end(), item) != keys.end();
}

bool isin(vertex_descriptor_mesh vit, std::set<vertex_descriptor_mesh> sep){
	return std::find(sep.begin(), sep.end(), vit) != sep.end();
}

int print_vertices_of_mesh(Polyhedron& mesh){
	vertex_iterator_mesh vit, ve;

	std::cout << "Vertices in first new mesh:" << num_vertices(mesh) << ". Edges in mesh: " << num_edges(mesh) << std::endl;
	for (boost::tie(vit, ve) = boost::vertices(mesh); vit != ve; ++vit)
		std::cout << (*vit)->id() << ", ";
	std::cout << std::endl;
	return 0;
}



int go_deeper(		const Polyhedron & mesh, 
					const std::set<vertex_descriptor_mesh> & vertex_separator, 
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

int label_components(	const Polyhedron & mesh, 
						const std::set<vertex_descriptor_mesh> & vertex_separator, 
						std::map<vertex_descriptor_mesh, int> &component,
						std::set<int> &values){


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


	//set all vertices in v_sep to be their own class
	for (boost::tie(vit, ve) = boost::vertices(mesh); vit != ve; ++vit)
		if (isin(*vit, vertex_separator))
			component[*vit] = -1;

	//isolate values for easy existance checking

	for (std::map<vertex_descriptor_mesh, int>::iterator k = component.begin(); k != component.end(); ++k){
		if (k->second != -1){
			values.insert(k->second);
		}
	}


	std::cout << "There are " << values.size() << " connected components" << std::endl;

	return 0;
}



int create_frame(std::vector<Point> &points, std::vector<E> &edges){
	//this will simulate the act of reading the landmarks from a graph
	//take a filename, return a list of landmarks
	
	points.push_back(Point(5, 5, 5));
	points.push_back(Point(-5, 5, 5));
	points.push_back(Point(-5, -5, 5));
	points.push_back(Point(5, -5, 5));
	points.push_back(Point(5, 5, -5));
	points.push_back(Point(-5, 5, -5));
	points.push_back(Point(-5, -5, -5));
	points.push_back(Point(5, -5, -5));

	edges.push_back(E(0, 1));
	edges.push_back(E(1, 2));
	edges.push_back(E(2, 3));
	edges.push_back(E(3, 0));
	edges.push_back(E(4, 5));
	edges.push_back(E(5, 6));
	edges.push_back(E(6, 7));
	edges.push_back(E(7, 4));
	edges.push_back(E(0, 4));
	edges.push_back(E(1, 5));
	edges.push_back(E(2, 6));
	edges.push_back(E(3, 7));




	return 0;
}



int dijkstra_create_separator(Polyhedron &mesh, 
	std::set<vertex_descriptor_mesh> &vertex_separator, 
	std::vector<vertex_descriptor_mesh>& vert_descript, 
	std::vector<E> &edges){

	
	for (std::vector<E>::iterator it = edges.begin(); it != edges.end(); ++it){


		std::cout << "finding separator for edge (" << (*it).first << ", " << (*it).second << ")" << std::endl;

		std::vector<vertex_descriptor_mesh> predecessors(num_vertices(mesh)); // To store parents
		std::vector<Weight> distances(num_vertices(mesh)); // To store distances

		IndexMap_mesh indexMap = get(vertex_index, mesh);
		PredecessorMap_mesh predecessorMap(&predecessors[0], indexMap);
		DistanceMap_mesh distanceMap(&distances[0], indexMap);

		dijkstra_shortest_paths(mesh, vert_descript[(*it).first],
			distance_map(distanceMap).predecessor_map(predecessorMap).distance_combine(dist_combine()).distance_compare(dist_compare()));


		//add the path to the vertex_separator
		vertex_descriptor_mesh x = vert_descript[(*it).second];

		std::vector<vertex_descriptor_mesh> tmp_sep;

		while (x != vert_descript[(*it).first])
		{
			tmp_sep.push_back(x);
			x = predecessorMap[x];
		}
		
		tmp_sep.push_back(vert_descript[(*it).first]);

		//std::cout << "Temporary vertex separator:" << std::endl;
		for (std::vector<vertex_descriptor_mesh>::iterator tm_it = tmp_sep.begin(); tm_it != tmp_sep.end(); ++tm_it)
		{
			vertex_separator.insert((*tm_it));
			//std::cout << (*tm_it)->id() << "->";
		}
		//std::cout << std::endl;
		

	}//end for

	return 0;
} // dijkstra_create_separator


//now I will compartmentalize everything so it can be generalized.
int split_and_parameterize_mesh(Polyhedron &mesh, std::vector<Point> &points, std::vector<E> &edges){
	
	
	vertex_iterator_mesh vit, ve;
	vertex_descriptor_mesh st, en, it;


	// Associate indices to the vertices


	//we want both the index and a reference to the actual vertex
	//the verts vector is unused - remove and put vert_descript in its place
	std::vector<int> verts;
	std::vector<vertex_descriptor_mesh> vert_descript;
	for (std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit){
		std::cout << "finding closest point to (" << (*pit).x() << ", " << (*pit).y() << ", " << (*pit).x() << ")" << std::endl;
		verts.push_back(naive_closest_vertex(mesh, (*pit), vert_descript));
	}






	//for (int i = 0; i < vert_descript.size(); i++)
	//	std::cout << "Item " << i << " is mesh vertex number " << vert_descript[i]->id() << std::endl;



	//create vertex separator
	std::set<vertex_descriptor_mesh> vertex_separator;
	dijkstra_create_separator(mesh, vertex_separator, vert_descript, edges);
	std::cout << "Our vertex_separator:" << std::endl;
	for (std::set<vertex_descriptor_mesh>::iterator tit = vertex_separator.begin(); tit != vertex_separator.end(); ++tit){
		std::cout << (*tit)->id() << "->";
	}
	std::cout << std::endl;

	// create a component map that labels each vertex by which component it is in, -1 if in the vertex separator
	std::map<vertex_descriptor_mesh, int> component_map;
	std::set<int> values;
	label_components(mesh, vertex_separator, component_map, values);

	

	


	std::vector<Polyhedron> new_meshes;
	std::vector<Parameterization_polyhedron_adaptor> mesh_adaptors;

	for (std::set<int>::iterator val_it = values.begin(); val_it != values.end(); ++val_it){
		new_meshes.push_back(partial_mesh_builder(mesh, component_map, (*val_it)));
	}

	for (std::vector<Polyhedron>::iterator p_it = new_meshes.begin(); p_it != new_meshes.end(); ++p_it){
		if (!(*p_it).is_valid()){
			std::cout << "cutting the mesh into pieces did not work" << std::endl;
		}
		else{

			//Parameterize our new meshes
			mesh_adaptors.push_back(parameterize_mesh(*p_it));
		}
	}

	char buffer[32]; // The filename buffer.
	int index_f = 0;
	for (std::vector<Parameterization_polyhedron_adaptor>::iterator ad_it = mesh_adaptors.begin(); ad_it != mesh_adaptors.end(); ++ad_it){
		// Put "file" then k then ".txt" in to filename.
		sprintf(buffer,  "../file%i.eps", index_f);
		++index_f;
		save_mesh((*ad_it), buffer);
	}

	return 0;
}


Polyhedron partial_mesh_builder(Polyhedron &mesh, std::map<vertex_descriptor_mesh, int> &component_map, int k){
	std::cout << "begin mesh builder" << std::endl;


	std::map<int, int> relative_index_map;
	std::set<int> keys;
	Polyhedron new_mesh;
	vertex_iterator_mesh vit, ve;
	v_handle n_v;
	int av, bv, cv;



	CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(new_mesh.hds(), true);
	B.begin_surface(100, 100, 300);


	//Add the vertices new this new mesh and create map between the two.
	int v_index = 0;
	for (boost::tie(vit, ve) = vertices(mesh); vit != ve; ++vit){
		if (component_map[(*vit)] == k || component_map[(*vit)] == -1){
			//std::cout << "placed at: " << v_index << " new id is: " << (*vit)->id() << std::endl;
			n_v = B.add_vertex((*vit)->point());
			relative_index_map[(*vit)->id()] = v_index;
			n_v->id() = (*vit)->id();
			++v_index;
		}
	}



	//create an easy accessor for the keys of the vertices we're using

	for (std::map<int, int>::iterator k = relative_index_map.begin(); k != relative_index_map.end(); ++k){
		keys.insert(k->first);
	}



	//Go through each triangle on our original mesh and if its vertices are in our new mesh, add the triangle
	for (Polyhedron::Facet_handle facet = mesh.facets_begin(); facet != mesh.facets_end(); ++facet){
		av = facet->halfedge()->vertex()->id();
		bv = facet->halfedge()->next()->vertex()->id();
		cv = facet->halfedge()->next()->next()->vertex()->id();

		//if this triangle in completely inside our new mesh
		if (isin(keys, av) && isin(keys, bv) && isin(keys, cv))
		{

			//std::cout << "add triangle: " << av << ", " << bv << ", " << cv << std::endl;
			//add this triangle
			B.begin_facet();
			B.add_vertex_to_facet(relative_index_map[av]);
			B.add_vertex_to_facet(relative_index_map[bv]);
			B.add_vertex_to_facet(relative_index_map[cv]);
			B.end_facet();
		}

	}



	B.end_surface();


	//remove all the extra vertices from vertex_separator that are not attached to our new surface.
	int deleted = new_mesh.keep_largest_connected_components(1);
	
	return new_mesh;

}





int test_real_images(){
	char* filename = "../../surface/surface.off";
	Polyhedron mesh;
	int fail = read_mesh(filename, mesh);
	std::cout << "fail if 1, sucess if 0: " << fail << std::endl;
	std::cout << "is valid: " << mesh.is_valid() << std::endl;
	std::cout << "number of vertices: " << num_vertices(mesh) << std::endl;
	return 0;
}


int main(int argc, char* argv[])
{

	std::vector<Point> points;
	std::vector<E> edges;
	create_frame(points, edges);

	char* filename = "../../surface/sphere.off";
	Polyhedron mesh;
	read_mesh(filename, mesh);





	split_and_parameterize_mesh(mesh, points, edges);
	//dijkstra_with_graph();

	//test_real_images();
	//take the shortest path found by dijkstra and find the connected components of mesh\separator


	//Parameterization_polyhedron_adaptor new_mesh = parameterize_mesh(mesh);
	//save_mesh(new_mesh);

	int h;
	std::cin >> h;
	return 0;
}

