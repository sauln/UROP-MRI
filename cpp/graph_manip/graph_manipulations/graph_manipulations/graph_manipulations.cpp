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
#include <climits>

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
//#include "Polyhedron_3_corners.h"
#include "Parameterization_polyhedron_adaptor_3_corners.h"
#include "Square_border_parameterizer_3_corners.h"
#include "graphs.h"
#include "tools.h"


using namespace boost;

//typedef std::pair<int, int> Ed;


// read the mesh from our file 
int read_mesh(char* filename, Polyhedron & mesh);
int print_mesh(const Polyhedron & mesh, const Parameterization_polyhedron_adaptor & mesh_adaptor);

int dijkstra_create_separator(Polyhedron &mesh,
std::set<vertex_descriptor_mesh> &vertex_separator,
std::vector<vertex_descriptor_mesh>& vert_descript,
std::vector<E> &edges);

int label_components(const Polyhedron & mesh,
	const std::set<vertex_descriptor_mesh> & vertex_separator,
	std::map<vertex_descriptor_mesh, int> &component,
	std::set<int> &values);

Polyhedron partial_mesh_builder(Polyhedron &mesh, std::map<vertex_descriptor_mesh, int> &component_map, int k);




//the frame defines our edges
struct Frame{
	std::vector<Point> points;
	std::vector<E> edges;
};

int read_mesh(char* filename, Polyhedron &mesh);


typedef std::map<vertex_iterator_m, Point3> p_map;


class meshes{
public:


	Frame												frame;
	Polyhedron											mesh;
	Polyhedron											standard_mesh;
	std::vector<vertex_descriptor_mesh>					corners;
	std::vector<Polyhedron>								original_meshes;
	std::vector<Parameterization_polyhedron_adaptor>	parameterized_meshes;

	


	std::vector<p_map>									p_maps;
	///these 3 should be deprecated
	//Polyhedron orig_meshes[6];
	//struct mesh_pair{
	//	Polyhedron *m_orig;
	//	Parameterization_polyhedron_adaptor  *m_param;
	//};
	//std::vector<mesh_pair> mesh_pairs;



	int find_opposites();
	meshes(char *f1, char *f2);
	int create_frame(char *filename);
	int split_and_parameterize_mesh();
	int test_corners();
	int find_corners();
	int set_4_corners(Polyhedron & o);
	Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron & m);
	int create_standard_mesh();
	p_map map_to_standard(Parameterization_polyhedron_adaptor& p_m, Polyhedron& o_m);
	Point3 find_closest_point(vertex_iterator_m & v, Parameterization_polyhedron_adaptor &p_m, Polyhedron& o_m);
	int raw_dump(Parameterization_polyhedron_adaptor& p_m, Polyhedron& o_m);
	double two_d_dist(double x1, double y1, double x2, double y2);


};

double meshes::two_d_dist(double x1, double y1, double x2, double y2){
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}


Point3 meshes::find_closest_point(vertex_iterator_m & v_c, Parameterization_polyhedron_adaptor &p_m, Polyhedron& o_m){

	//this function can be rewritten to either return the closest vertex, or find the barycentric coordinates of the point

	Point3 p_min;
	double min = DBL_MAX;


	Polyhedron::Vertex_const_iterator pVertex;
	for (pVertex = o_m.vertices_begin();
		pVertex != o_m.vertices_end();
		pVertex++)
	{


		double u = p_m.info(pVertex->halfedge())->uv().x();
		double v = p_m.info(pVertex->halfedge())->uv().y();

		double d = two_d_dist(v_c->point().x(), v_c->point().y(), u, v);
		if (d < min){
			min = d;
			p_min = pVertex->point();
		}


		//std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
	}


	//std::cout << "Closest point at " << p_min << " is " << min << " far away." << std::endl;

	return p_min;

}


int meshes::raw_dump(Parameterization_polyhedron_adaptor& p_m, Polyhedron& o_m){
	// Raw output: dump (u,v) pairs
	Polyhedron::Vertex_const_iterator pVertex;
	for (pVertex = o_m.vertices_begin();
		pVertex != o_m.vertices_end();
		pVertex++)
	{
		// (u,v) pair is stored in any halfedge
		double u = p_m.info(pVertex->halfedge())->uv().x();
		double v = p_m.info(pVertex->halfedge())->uv().y();
		std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
	}
	return 0;
}

int meshes::find_opposites(){
	//this is just for myself so I can see which to deal with.  after I figure this out i'll probably hardcode it. 

	//we have the 8 corners of our original mesh in  std::vector<vertex_descriptor_mesh>::iterator it = corners.begin()






	for (std::vector<vertex_descriptor_mesh>::iterator it = corners.begin(); it != corners.end(); ++it){
		std::cout << (*it)->id() << " at " <<(*it)->point() << std::endl;
	}
	std::vector<Polyhedron>::iterator tmp = original_meshes.begin();
	std::vector<int>::iterator tmp_c = tmp->corners.begin();
	std::cout << std::endl << std::endl;
	

		
		

	std::set<int> c_set;
	while (tmp_c != tmp->corners.end()){
		c_set.insert(*tmp_c);
		tmp_c++;
	}

	std::set<int>::iterator set_it;
	for (set_it = c_set.begin(); set_it != c_set.end(); ++set_it){
		std::cout << *set_it << ", ";
 	}


	tmp++;
	Polyhedron *opposite = &(*original_meshes.begin()) ; 
	std::cout << std::endl;
	while (tmp != original_meshes.end()){ 
		std::cout << "Corners in this mesh: ";
		for (tmp_c = tmp->corners.begin(); tmp_c != tmp->corners.end(); ++tmp_c){
			std::cout << *tmp_c << ", ";
		}
		std::cout << std::endl;


		for (tmp_c = tmp->corners.begin(); tmp_c != tmp->corners.end(); ++tmp_c){
			if (c_set.end() != c_set.find(*tmp_c)){
				std::cout << "We think that " << *tmp_c << " is in the set" << std::endl;
				std::cout << "Break and move onto next slab" << std::endl;
				break;
			}
			std::vector<int>::iterator peak = tmp_c;
			peak++;
			if (peak == tmp->corners.end()){//then we make it through all corners
				std::cout << "Found opposite slab" << std::endl;
				opposite = &(*tmp);
			}
		}
		tmp++;//cycle through each mesh.
	}
	
	std::cout << "The slab opposite is defined by the corners: ";
	for (tmp_c = opposite->corners.begin(); tmp_c != opposite->corners.end(); ++tmp_c){
		std::cout << *tmp_c << ", " ;

	}
		


	return 0;
}


p_map meshes::map_to_standard(Parameterization_polyhedron_adaptor& p_m, Polyhedron& o_m){
	//check to see if we have created_standard_mesh()
	std::cout << "Write a check to see if we have already ran create_standard_mesh()" << std::endl;

	//std::map<vertex_iterator_m, Point3> p_map;

	

	p_map tmp;

	for (vertex_iterator_m vit = standard_mesh.vertices_begin(); vit != standard_mesh.vertices_end(); ++vit){
		//find closest vertex in parameterized_meshes.front() to vit
		tmp[vit] = find_closest_point(vit, p_m, o_m);
	}


	return tmp;



}

int meshes::find_corners(){
	//find our vertices that are closest to the frame.
	for (std::vector<Point>::iterator pit = frame.points.begin(); pit != frame.points.end(); ++pit){
		//std::cout << "finding closest point to (" << (*pit).x() << ", " << (*pit).y() << ", " << (*pit).x() << ")" << std::endl;
		naive_closest_vertex(mesh, (*pit), corners);
	}
	return 0;
}

int meshes::split_and_parameterize_mesh(){
	//std::vector<Point> &points, std::vector<E> &edges,

	vertex_iterator_mesh						vit, ve;
	vertex_descriptor_mesh						st, en, it;

	std::set<vertex_descriptor_mesh>			vertex_separator;
	std::map<vertex_descriptor_mesh, int>		component_map;
	std::set<int>								values;

	//create vertex separator
	dijkstra_create_separator(mesh, vertex_separator, corners, frame.edges);


	// create a component map that labels each vertex by which component it is in, -1 if in the vertex separator
	label_components(mesh, vertex_separator, component_map, values);


	std::vector<Polyhedron> original_m_tmp;
	for (std::set<int>::iterator val_it = values.begin(); val_it != values.end(); ++val_it){

		//create new meshes from the connected components
		original_meshes.push_back(partial_mesh_builder(mesh, component_map, (*val_it)) );




		set_4_corners(original_meshes.back());



		if (original_meshes.back().is_valid()){
			parameterized_meshes.push_back(parameterize_mesh(original_meshes.back()));			
		}
	}

	return 0;
}

int meshes::set_4_corners(Polyhedron & o){
	//This function traverses the new mesh and finds the id of the corners 
	for (vertex_iterator_m it = o.vertices_begin(); it != o.vertices_end(); ++it){
		for (std::vector<vertex_descriptor_mesh>::iterator vit = corners.begin();
			vit != corners.end(); ++vit){
			if ((*it).point() == (*vit)->point()){
				o.corners.push_back((*it).id());
			}
		}
	}

	if (o.corners.size() != 4){ std::cout << "THERE ARE NOT FOR CORNERS ON OUR SQUARE MESH!!" << std::endl; }
	return 0;
}

int meshes::create_frame(char* filename){
	//this will simulate the act of reading the landmarks from a graph
	//take a filename, return a list of landmarks
	//std::vector<Point> &points, std::vector<E> &edges


	//				/--a
	//			d--/   |   b
	//			|	c  |   |
	//			|   |  |   |
	//			|	|  e   |
	//			h	|	   f
	//				g

	// figure out how to read this in from a file
	std::cout << "Creating frame" << std::endl;

	double a, b, c, d;
	std::ifstream infile(filename);

	while (infile >> a >> b >> c >> d){
		//std::cout << a << b << c << d << std::endl;
		frame.points.push_back(Point(b, c, d));
	}

	frame.edges.push_back(E(0, 1));
	frame.edges.push_back(E(1, 2));
	frame.edges.push_back(E(2, 3));
	frame.edges.push_back(E(3, 0));
	frame.edges.push_back(E(4, 5));
	frame.edges.push_back(E(5, 6));
	frame.edges.push_back(E(6, 7));
	frame.edges.push_back(E(7, 4));
	frame.edges.push_back(E(0, 4));
	frame.edges.push_back(E(1, 5));
	frame.edges.push_back(E(2, 6));
	frame.edges.push_back(E(3, 7));




	return 0;
}

meshes::meshes(char* f1, char* f2)
{
	parameterized_meshes.reserve(6);
	original_meshes.reserve(6);
	p_maps.reserve(6);

	create_frame(f2);
	read_mesh(f1, mesh);
	find_corners();
	split_and_parameterize_mesh();
	create_standard_mesh();


}

int meshes::test_corners(){
	//I need to figure out some way to test if the parameterizations are correct...

	std::cout << "Do the meshes work? and if so, how do they work?" << std::endl;
	vertex_iterator_param vit;

	std::cout << "Each corner should be one of 3 faces" << std::endl;


	std::cout << original_meshes[0].vertices_begin()->point() << std::endl
		<< original_meshes[1].vertices_begin()->point() << std::endl
		<< original_meshes[2].vertices_begin()->point() << std::endl
		<< original_meshes[3].vertices_begin()->point() << std::endl
		<< original_meshes[4].vertices_begin()->point() << std::endl
		<< original_meshes[5].vertices_begin()->point() << std::endl;


	return 0;


}


//Parameterize the map and return the mesh adaptor
Parameterization_polyhedron_adaptor meshes::parameterize_mesh(Polyhedron & m ){
	//std::cout << "parameterize mesh and return parameterization" << std::endl;
	// Create Polyhedron adaptor mesh
	// Assume mesh is a topological disk
	
	//int* foo = new int[10];
	Parameterization_polyhedron_adaptor mesh_adaptor(m);
	//Parameterization_polyhedron_adaptor mesh_adaptor(m);

	//add the 4 corners to this structure
	for (std::vector<int>::iterator ait = m.corners.begin(); ait != m.corners.end(); ++ait)
		mesh_adaptor.corners.push_back(*ait);
	
	if (mesh_adaptor.corners.size() != 4)
		std::cout << "DIFFERENT THAN 4 CORNERS SHOWED UP!: " << mesh_adaptor.corners.size() << std::endl;
	

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


	//std::cout << "The output works inside of parameterize_mesh() " << std::endl;







	return mesh_adaptor;
}










//Read in a mesh from a file.
int read_mesh(char* filename, Polyhedron &mesh){

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


Polyhedron partial_mesh_builder(Polyhedron &mesh, std::map<vertex_descriptor_mesh, int> &component_map, int k);

static bool write_file_eps(const Parameterization_polyhedron_adaptor& mesh_adaptor,
	const char *pFilename,
	double scale);
int read_mesh(char* filename, Polyhedron & mesh);


int save_mesh(Parameterization_polyhedron_adaptor & mesh_adaptor, char * filename);
int save_parameterized_meshes(meshes mesh_set);

int print_mesh(const Polyhedron & mesh, const Parameterization_polyhedron_adaptor & mesh_adaptor);
int print_vertices_of_mesh(Polyhedron& mesh);

bool isin(int item, std::set<int> keys);
bool isin(vertex_descriptor_mesh vit, std::set<vertex_descriptor_mesh> sep);

int create_frame(Frame &frame);

int label_deeper(const Polyhedron & mesh,
	const std::set<vertex_descriptor_mesh> & vertex_separator,
	std::map<vertex_descriptor_mesh, int>& component,
	int & component_number, vertex_descriptor_mesh vit);








int test_real_images();









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


int save_parameterized_meshes(meshes mesh_set){
	std::cout << "Try saving the meshes" << std::endl;

	char buffer[64]; // The filename buffer.
	int index_f = 0;
	for (std::vector<Parameterization_polyhedron_adaptor>::iterator ad_it = mesh_set.parameterized_meshes.begin(); ad_it != mesh_set.parameterized_meshes.end(); ++ad_it){
		// Put "file" then k then ".txt" in to filename.
		sprintf(buffer, "../../surface/output_surfaces/file%i.eps", index_f);
		++index_f;
		save_mesh((*ad_it), buffer);
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
bool isin(int item, std::set<int> keys){
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

int label_deeper(		const Polyhedron & mesh, 
					const std::set<vertex_descriptor_mesh> & vertex_separator, 
					std::map<vertex_descriptor_mesh, int>& component, 
					int & component_number, vertex_descriptor_mesh vit){

	//this is the poorly named recursive part of the breadth first search.
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
				label_deeper(mesh, vertex_separator, component, component_number, fish);
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
		check = label_deeper(mesh, vertex_separator, component, component_number, *vit);
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


int dijkstra_create_separator(Polyhedron &mesh, 
	std::set<vertex_descriptor_mesh> &vertex_separator, 
	std::vector<vertex_descriptor_mesh>& vert_descript, 
	std::vector<E> &edges){

	
	for (std::vector<E>::iterator it = edges.begin(); it != edges.end(); ++it){


		//std::cout << "finding separator for edge (" << (*it).first << ", " << (*it).second << ")" << std::endl;
		//std::cout << "        -> that is (" << vert_descript[(*it).first]->id() <<", " << 
		//		vert_descript[(*it).second]->id()<< ")" << std::endl;

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

		while (x != vert_descript[(*it).first]){
			tmp_sep.push_back(x);
			x = predecessorMap[x];
		}
		
		tmp_sep.push_back(vert_descript[(*it).first]);

		//std::cout << "Temporary vertex separator:" << std::endl;
		for (std::vector<vertex_descriptor_mesh>::iterator tm_it = tmp_sep.begin(); tm_it != tmp_sep.end(); ++tm_it){
			vertex_separator.insert((*tm_it));
			//std::cout << (*tm_it)->id() << "->";
		}
		//std::cout << std::endl;
		

	}//end for

	return 0;
} // dijkstra_create_separator




Polyhedron partial_mesh_builder(Polyhedron &mesh, std::map<vertex_descriptor_mesh, int> &component_map, int k){
	/// 
	/// Create a new mesh for the component number k according to the component map
	///
	///



	std::cout << "create new mesh for face " << k  << std::endl;
	std::map<int, int> relative_index_map;
	std::set<int> keys;
	Polyhedron new_mesh;
	vertex_iterator_mesh vit, ve;
	v_handle n_v;
	int av, bv, cv;

	CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(new_mesh.hds(), true);
	//This is just big enough now
	//TODO errors can be here
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
		if (isin(av, keys) && isin(bv, keys) && isin(cv, keys))
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



int meshes::create_standard_mesh(){
	//create a polyhedron:
	//we want vertices on [0,1]x[0,1]

	std::cout << "Create a stnadardized mesh" << std::endl;

	CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(standard_mesh.hds(), true);
	//This is just big enough now
	//TODO errors can be here

	double step = 0.1;


	B.begin_surface(100, 100, 300);

	for (double i = 0; i <= 1; i += step){
		std::cout << "Add row " << i << std::endl;
		for (double j = 0; j <= 1; j += step){
			B.add_vertex(Point(i, j, 0.0));
		}

	}
	
	B.end_surface();

	return 0;

}

int main(int argc, char* argv[])
{
	char* f1 = "../../surface/sphere.off";
	char* f2 = "../../surface/landmarks.txt";

	meshes mesh_set(f1, f2);
	mesh_set.test_corners();


	

	//for (vertex_iterator_m it = mesh_set.standard_mesh.vertices_begin(); it != mesh_set.standard_mesh.vertices_end(); ++it)
	//	std::cout << it->point() << std::endl;

	for (int i = 0; i < 6; i++){
		std::cout << "creating map for slab " << i << std::endl;
		mesh_set.p_maps.push_back(mesh_set.map_to_standard(mesh_set.parameterized_meshes[i], mesh_set.original_meshes[i]));
	}







	mesh_set.find_opposites();










	std::cout << "End" << std::endl;


	int h;
	std::cin >> h;
	return 0;
}

