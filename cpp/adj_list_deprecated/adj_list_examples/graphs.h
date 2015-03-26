


#ifndef GRAPHS
#define GRAPHS

#include <boost/array.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>





using namespace boost;








/* Properties for default graph */
struct V { double x; double y; double z; };
typedef std::pair<std::size_t, std::size_t> E;
typedef double Weight;
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
default_Graph create_graph(char* filename);
double diff_V(V p, V q);
int print_all_edges(default_Graph g);
int add_edge_N(int in, int out, default_Graph &g);






default_Graph create_graph(char* filename)
{
	int		num_verts;
	double	a, b, c, d;

	for (int i = 0; i < 9; i++)
	{
		std::cout << filename[i];
	}

	//std::cout << "filename is " << filename << "\n";

	std::ifstream infile(filename);
	std::string line; getline(infile, line);


	infile >> num_verts >> b >> c;
	default_Graph g(num_verts);

	//read and save vertex properties
	for (int i = 0; i < num_verts; i++){
		infile >> a >> b >> c;
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


#endif