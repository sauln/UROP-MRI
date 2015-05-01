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

class meshes{
public:
	//meshes();
	meshes(char *f1, char *f2);
	int create_frame(char *filename);
	int split_and_parameterize_mesh();
	int test_corners();
	int find_corners();
	int set_4_corners(Polyhedron & o);
	//int read_mesh(char* filename);





	//Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron &m);
	//int parameterize_mesh(Polyhedron & m, Parameterization_polyhedron_adaptor &mesh_adaptor);
	Parameterization_polyhedron_adaptor parameterize_mesh(Polyhedron & m);
	Frame												frame;
	Polyhedron											mesh;

	///
	//std::vector<cuts>  mesh_cuts;

	std::vector<vertex_descriptor_mesh>					corners;
	std::vector<Polyhedron>								original_meshes;
	std::vector<Parameterization_polyhedron_adaptor>	parameterized_meshes;


	//Parameterization_polyhedron_adaptor param_meshes[6];
	Polyhedron orig_meshes[6];

	struct mesh_pair{
		Polyhedron *m_orig;
		Parameterization_polyhedron_adaptor  *m_param;
	};
	std::vector<mesh_pair> mesh_pairs;


};

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

		//Pointer to the mesh
		//Polyhedron& o_tmp = original_meshes.back();
		set_4_corners(original_meshes.back());
		//mesh_pair tmp_pair;
		if (original_meshes.back().is_valid()){

			
			//mesh_adaptor is built off a mesh
			//we need the mesh to persist in the same place
			 // new Parameterization_polyhedron_adaptor(original_meshes.back());


			//TODO start here
			//These still aren't copying




			//only the most recent mesh_adaptor exists?!
			//Parameterization_polyhedron_adaptor mesh_adaptor = ;
			parameterized_meshes.push_back(parameterize_mesh(original_meshes.back()));
			
		}

		std::cout << "But the data cannot be accessed on the outside?  Why is the mesh not coming through?" << std::endl;
		std::cout << "Size of our mesh_pairs:" << mesh_pairs.size() << std::endl;


		///Issues:
		//  Parameterization_polyhedron_adaptor has no default constructor
		//  how do I get it out of the parameterize_mesh where it is initiated?




		//errors occur when we try to access the mesh_adaptor->halfedge()

		//Parameterization_polyhedron_adaptor *p_mesh = mesh_pairs[0].m_param;//mesh_pairs[0]->m_param;
		//Polyhedron *p_orig = mesh_pairs[0].m_orig;
		std::cout << "Try the most recent one" << std::endl;
		// Raw output: dump (u,v) pairs
		Polyhedron::Vertex_const_iterator pVertex;
		for (pVertex = original_meshes.back().vertices_begin();
			pVertex != original_meshes.back().vertices_end();
			pVertex++)
		{
			//double u, v;
			//u = 0.0; v = 1.0;
			// (u,v) pair is stored in any halfedge
			//parameterized_meshes[0].info(pVertex->halfedge())->uv().x();
			double u = parameterized_meshes.back().info(pVertex->halfedge())->uv().x();
			double v = parameterized_meshes.back().info(pVertex->halfedge())->uv().y();
			//std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
		}


		std::cout << "Try the first one" << std::endl;
		// Raw output: dump (u,v) pairs
		std::cout << "Only the most recent one works" << std::endl;
		for (pVertex = original_meshes.front().vertices_begin();
			pVertex != original_meshes.front().vertices_end();
			pVertex++)
		{
			//double u, v;
			//u = 0.0; v = 1.0;
			// (u,v) pair is stored in any halfedge
			std::cout << pVertex->halfedge()->vertex()->point() << std::endl;
			//parameterized_meshes[0].info(pVertex->halfedge())->uv().x();
			double u = parameterized_meshes.front().info(pVertex->halfedge())->uv().x();
			double v = parameterized_meshes.front().info(pVertex->halfedge())->uv().y();
			//std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
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
	create_frame(f2);
	read_mesh(f1, mesh);
	find_corners();
	split_and_parameterize_mesh();

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


	// Raw output: dump (u,v) pairs
	Polyhedron::Vertex_const_iterator pVertex;
	for (pVertex = original_meshes[0].vertices_begin();
		pVertex != original_meshes[0].vertices_end();
		pVertex++)
	{
		// (u,v) pair is stored in any halfedge
		//parameterized_meshes[0].info(pVertex->halfedge())->uv().x();
		double u = parameterized_meshes[0].info(pVertex->halfedge())->uv().x();
		double v = parameterized_meshes[0].info(pVertex->halfedge())->uv().y();
		std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
	}



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
	


	//for each vertex in our original mesh
	//for (vertex_iterator_m vit = m.vertices_begin(); vit != m.vertices_end(); ++vit){
		//and for each corner of this original mesh
		//std::cout << vit->id() << std::endl;
	//	for (std::vector<int>::iterator it = m.corners.begin(); it != m.corners.end(); ++it){
			//if this vertex is one of the corners (check by id)
	//		if (vit->id() == m.corners[0]){
				//then for each vertex in the adaptor
				//std::cout << "found corner" << std::endl;
	//			for (vertex_iterator_m mit = mesh_adaptor.mesh_vertices_begin(); mit != mesh_adaptor.mesh_vertices_end(); ++mit){
					//if the vertex in our original mesh and the vertex in the adaptor are the same point
	//				if (mit->point() == vit->point()){
						//then this vertex in our adaptor must be a corner!
	//					mesh_adaptor.corners.push_back(mit->id());
					//}}}}}






	//std::cout << std::endl;
	//std::cout << "New mesh vertex #" << m.corners[0] << std::endl;
	//std::cout << "Total vertex count: " << mesh_adaptor.count_mesh_vertices() << std::endl;
	//std::cout << "Mesh corners: ";
	//for (std::vector<int>::iterator ait = mesh_adaptor.corners.begin(); ait != mesh_adaptor.corners.end(); ++ait)
	//	std::cout << *ait << ", ";
	//std::cout << std::endl;
	//mesh_adaptor.corners_id.push_back(3);
	// we need to set the corners

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


	std::cout << "The output works inside of parameterize_mesh() " << std::endl;
	// Raw output: dump (u,v) pairs
	Polyhedron::Vertex_const_iterator pVertex;
	for (pVertex = m.vertices_begin();
		pVertex != m.vertices_end();
		pVertex++)
	{
		// (u,v) pair is stored in any halfedge
		double u = mesh_adaptor.info(pVertex->halfedge())->uv().x();
		double v = mesh_adaptor.info(pVertex->halfedge())->uv().y();
		//std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
	}






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



int main(int argc, char* argv[])
{
	char* f1 = "../../surface/sphere.off";
	char* f2 = "../../surface/landmarks.txt";

	meshes mesh_set(f1, f2);
	mesh_set.test_corners();

	
	std::cout << "finish parametrization" << std::endl;



	std::cout << "Test parameterizations" << std::endl;
	//save_parameterized_meshes(mesh_set);

	std::cout << "End" << std::endl;




	int h;
	std::cin >> h;
	return 0;
}

