


#ifndef GRAPHS
#define GRAPHS


#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/filtered_graph.hpp>



using namespace boost;



typedef double Weight;

////////////////////////////////////////////////////////////////////////////////////////

//Polyhedron types

typedef CGAL::Simple_cartesian<double>					Kernel;
typedef CGAL::Vector_3<Kernel>							Vector;
typedef Kernel::Point_3									Point;


typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;


//typedef Polyhedron::Vertex_iterator									vertex_iterator_mesh_2;

//typedef boost::graph_traits<Polyhedron>::vertex_descriptor			vertex_descriptor_mesh;
typedef boost::graph_traits<Polyhedron>::vertex_iterator			vertex_iterator_mesh;


//typedef boost::graph_traits<Polyhedron>::edge_descriptor			edge_descriptor_mesh;
typedef Polyhedron::Halfedge_handle edge_descriptor_mesh;


//typedef boost::graph_traits<Polyhedron>::edge_iterator edge_iterator_mesh;
typedef Polyhedron::Edge_iterator									edge_iterator_mesh;
typedef Polyhedron::Facet_iterator									facet_iterator_mesh;

typedef Polyhedron::Facet_handle									facet_descriptor_mesh;
typedef Polyhedron::Halfedge_handle									halfedge_handle_mesh;
typedef Polyhedron::Point_3											point_mesh;
typedef Polyhedron::HalfedgeDS										HalfedgeDS;
typedef HalfedgeDS::Vertex											Vertex;
typedef Vertex::Point												Point3;

typedef Polyhedron::Vertex_iterator									vertex_iterator_m;
typedef HalfedgeDS::Vertex_handle									v_handle;
typedef HalfedgeDS::Vertex_handle									vertex_descriptor_mesh;
//typedef Polyhedron::Vertex_iterator vertex_iterator_mesh;
typedef CGAL::Out_edge_iterator < Polyhedron >						edge_it_mesh;


// Floater Mean Value Coordinates parameterization with square border
// Parameterization Mesh
typedef CGAL::Parameterization_polyhedron_adaptor_3_corners<Polyhedron>		Parameterization_polyhedron_adaptor;
typedef CGAL::Square_border_arc_length_parameterizer_3_corner<Parameterization_polyhedron_adaptor>  
																	Border_parameterizer;
typedef CGAL::Eigen_solver_traits<>									Solver;


typedef CGAL::Mean_value_coordinates_parameterizer_3<Parameterization_polyhedron_adaptor, Border_parameterizer, Solver> 
																	Parameterizer;

typedef Parameterization_polyhedron_adaptor::Vertex_iterator		vertex_iterator_param;



/* Property maps for Dijkstra's shortest path */
typedef property_map < Polyhedron, vertex_index_t >::type			IndexMap_mesh;
typedef iterator_property_map < vertex_descriptor_mesh*, IndexMap_mesh, vertex_descriptor_mesh, vertex_descriptor_mesh& > 
																	PredecessorMap_mesh;
typedef iterator_property_map < Weight*, IndexMap_mesh, Weight, Weight& > 
																	DistanceMap_mesh;

//compare and combine used for the dijkstra's shortest path
struct dist_combine{
	double operator()(const double & a, const double &b) const {
		return sqrt(b) + a;
	}
};

struct dist_compare{
	bool operator()(const double & a, const double &b) const{
		return sqrt(a) < sqrt(b);
	}
};

template <typename TGraph>
struct vertex_id_filter_mesh
{
	//add an attribute for the vertex_separator
	std::vector<vertex_descriptor_mesh> vertex_separator;
	int size_of_vs;

	//predicate to check if the vertex is in our vertex_separator.  Intended to remove all 
	//vertices in separator.
	bool operator()(const typename boost::graph_traits<TGraph>::vertex_descriptor& v) const
	{
		//for unknown reason iterators are not working- resort to this ghetto way of checking 
		for (int i = 0; i < size_of_vs; i++)
			if (vertex_separator[i] == v)
				return false;

		return true; //default to returning true if v is not in the separator. (keep those not in separator
	}
};

typedef filtered_graph<Polyhedron, keep_all, vertex_id_filter_mesh<Polyhedron> > FilteredMeshType;



			//typedef Is_finite<Polyhedron> Filter;
			//typedef boost::filtered_graph<Polyhedron, Filter, Filter> Finite_Polyhedron;
			//typedef boost::graph_traits<Finite_Polyhedron>::vertex_descriptor vertex_descriptor_fin;
			//typedef boost::graph_traits<Finite_Polyhedron>::vertex_iterator vertex_iterator_fin;









/* Predicate and filter for the separator filter */

















//  Deprecated 



template <typename T>
struct Is_finite {

	const T* t_;

	Is_finite()
		: t_(NULL)
	{}

	Is_finite(const T& t)
		: t_(&t)
	{ }

	template <typename VertexOrEdge>
	bool operator()(const VertexOrEdge& voe) const {
		return !t_->is_infinite(voe);
	}
};

/* Properties for default graph */
struct V {
	double x;
	double y;
	double z;

	//bool operator=(const my_data_type& other) const { return my_i < other.my_i; }

};
typedef std::pair<std::size_t, std::size_t> E;

typedef int Name;
typedef property<edge_weight_t, Weight> WeightProperty;

typedef adjacency_list < setS, vecS, undirectedS,
	V, WeightProperty, no_property, vecS > default_Graph;

/* Descriptors and iterators */
typedef graph_traits < default_Graph >::vertex_descriptor vertex_descriptor;
typedef graph_traits < default_Graph >::vertex_iterator vertex_iterator;
typedef graph_traits <default_Graph>::edge_iterator EdgeIterator;
typedef std::pair<EdgeIterator, EdgeIterator> EdgePair;

/* Property maps for Dijkstra's shortest path */
typedef property_map < default_Graph, vertex_index_t >::type IndexMap;
typedef iterator_property_map < vertex_descriptor*, IndexMap, vertex_descriptor, vertex_descriptor& > PredecessorMap;
typedef iterator_property_map < Weight*, IndexMap, Weight, Weight& > DistanceMap;

/* Predicate and filter for the separator filter */
template <typename TGraph>
struct vertex_id_filter
{
	//add an attribute for the vertex_separator
	std::vector<vertex_descriptor> vertex_separator;
	int size_of_vs;

	//predicate to check if the vertex is in our vertex_separator.  Intended to remove all 
	//vertices in separator.
	bool operator()(const typename boost::graph_traits<TGraph>::vertex_descriptor& v) const
	{
		//for unknown reason iterators are not working- resort to this ghetto way of checking 
		for (int i = 0; i < size_of_vs; i++)
			if (vertex_separator[i] == v)
				return false;

		return true; //default to returning true if v is not in the separator. (keep those not in separator
	}
};

typedef filtered_graph<default_Graph, keep_all, vertex_id_filter<default_Graph> > FilteredGraphType;



/* Interface for graph functions */
/* One day these might be turned into a class */
int read_vertex(default_Graph g, int position);
int write_vertex(default_Graph &g, int position, double x, double y, double z);
int write_point(V &p, double x, double y, double z);
default_Graph read_graph(char* filename);
double diff_V(V p, V q);
V diff_direction(V a, V b);
int print_all_edges(default_Graph g);
int add_edge_N(int in, int out, default_Graph &g);
int print_point(V &p);






V project_ish(double scale, double portion, V direction, double dir_scale, V a);
V project_ish(double scale, double portion, V direction, double dir_scale, V a)
{

	std::cout << "Scale: " << scale << std::endl
		<< "portion: " << portion << std::endl
		<< "dir_scale: " << dir_scale << std::endl;
	std::cout << "Direction: " << std::endl;
	print_point(direction);
	std::cout << "Start: " <<  std::endl;
	print_point(a);
		 
	V n;
	double scaling =  portion / scale ;
	n.x = scaling * (direction.x / dir_scale) + a.x;
	n.y = scaling * (direction.y / dir_scale) + a.y;
	n.z = scaling * (direction.z / dir_scale) + a.z;

	return n;

}

default_Graph read_graph(char* filename)
{
	int		num_verts;
	double	a, b, c, d;

	std::ifstream infile(filename);
	
	std::string line; getline(infile, line);


	infile >> num_verts >> b >> c;
	default_Graph g(num_verts);

	//read and save vertex properties
	for (int i = 0; i < num_verts; i++){
		infile >> a >> b >> c;
		//write_vertex(graph, index, x_coord, y_coord, z_coord);
		write_vertex(g, i, a, b, c);
	}

	//read and save each triangle as 3 edges
	int q, w, e, r;
	while (infile >> q >> w >> e >> r){
		add_edge_N(w, e, g);
		add_edge_N(e, r, g);
		add_edge_N(r, w, g);
	}

	return g;
}
int write_point(V &p, double x, double y, double z){
	p.x = x; p.y = y; p.z = z;
	return 0;
}
int read_vertex(default_Graph g, int position){
	std::cout << "Vertex is located at: (" << g[position].x << "," << g[position].y << "," << g[position].z << ")\n";
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
V diff_direction(V a, V b){
	/*unfortunately, V is of length 3*/

	double norm = diff_V(a, b);
	V sum;
	sum.x = (a.x - b.x) / norm;
	sum.y = (a.y - b.y) / norm;
	sum.z = (a.z - b.z) / norm;


	return sum;
}
int add_edge_N(int in, int out, default_Graph &g){
	//Adds edge with ordered vertices and euclidean distance weight
	double w = diff_V(g[in], g[out]);

	if (in <= out)
		add_edge(in, out, w, g);
	else
		add_edge(out, in, w, g);
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
int print_point(V &p){
	std::cout << "Point is (" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
	return 0;
}
int print_all_vertices(default_Graph g){
	
	for (unsigned int i = 0; i < num_vertices(g); i++)
		print_point(g[i]);

	return 0;

}



#endif