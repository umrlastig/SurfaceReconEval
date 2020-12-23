#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Point_set_3<Point> Point_set;

// shortest path
typedef CGAL::Surface_mesh_shortest_path_traits<K, Mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;

// iterators and descriptors
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

// accelerating data structure
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Tree_traits;
typedef CGAL::AABB_tree<Tree_traits> Tree;



Mesh read_OBJ_mesh(const char* fileName){
	std::cout << "---> READING: " << fileName << std::endl;
	Mesh mesh;
	// Variables to store polygon soup:
	std::ifstream is (fileName);
	std::vector<Point> points;
	std::vector<std::vector<std::size_t>> polygons;

	// Extract points and polygons from stream 'is':
	CGAL::read_OBJ(is, points, polygons);
	std::cout << "number of points: " << points.size() << std::endl;
	std::cout << "number of polygons: " << polygons.size() << std::endl;

	// Processing polygon soup:
	// CGAL::Polygon_mesh_processing::repair_polygon_soup 	(
	if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons))
	{
		// polygons do not define a valid mesh
		std::cout << "Warning: [" << fileName << "] does not define a valid polygon mesh " << std::endl;
		if ( !CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons) ){ // try to orient
			// orientation failed
			std::cout << "Warning: the orientation of the polygon soup ["
			<< fileName << "] failed. => Points were duplicated." << std::endl;
		} else {
			// orientation succeeded
			std::cout << "orientation of the polygon soup [" << fileName << "] was a success." << std::endl;
		}
	} else {
		// polygons define a valid mesh
		std::cout << "Success: [" << fileName << "] defines a valid polygon mesh " << std::endl;
	}

	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

	std::cout << "Number of points: " << mesh.number_of_vertices() << std::endl;
	std::cout << "Number of faces: " << mesh.number_of_faces() << std::endl;
	std::cout << "Done with reading" << std::endl << std::endl;

	return mesh;
}

Surface_mesh_shortest_path build_shortest_path(Mesh &mesh, Point_set &pcd){
	std::cout << "---> Building Shortest Path Sequence Tree" << std::endl;
	Surface_mesh_shortest_path shortest_paths(mesh);
	// Tree tree(faces(mesh).first, faces(mesh).second, mesh); // AABB_tree of the mesh
	Tree tree;
	shortest_paths.build_aabb_tree(tree); // Creates an AABB_tree suitable for use with locate
	// browse point set:
	for (Point_set::const_iterator pi = pcd.begin(); pi != pcd.end(); ++pi){
		// add point *pi to source:
		shortest_paths.add_source_point( shortest_paths.locate(pcd.point(*pi), tree) );
	}
	std::cout << " Number of source points: " << shortest_paths.number_of_source_points() << std::endl;
	shortest_paths.build_sequence_tree();
	std::cout << " Done with building sequence tree" << std::endl;
	return shortest_paths;
}

std::array<vertex_descriptor, 3> face_alpha(std::map<vertex_descriptor, vertex_descriptor> &mapRefAlpha,
											std::array<vertex_descriptor, 3> &face_ref)
{
	/*
		usage:
			- get the correspondance between vertex descriptors of two meshes
			  where one mesh is a sub-mesh of the other one.
		input:
			- mapRefAlpha: the map correspondance of vertex descriptors:
				o key: vertex descriptor of vertex in mesh_ref
				o value: vertex descriptor of vertex in sub-mesh mesh_alpha
			- face_ref: array of 3 vertex descriptors defining a face in
			  mesh_ref
		output:
			- face_a: array of corresponding 3 vertex descriptors defining
			  the same face but in sub-mesh mesh_alpha
	*/
	std::map<vertex_descriptor, vertex_descriptor>::iterator it_ref;
	std::array<vertex_descriptor, 3> face_a;
	for (int i=0; i != 3; i++){
		vertex_descriptor v_ref = face_ref[i];
		vertex_descriptor v_alpha = mapRefAlpha[v_ref];
		face_a[i] = v_alpha;
	}
	return face_a;
}

Mesh Mesh_alpha(Mesh &mesh, Point_set &pcd, double alpha){
	/*
		usage:
			- compute the "reconstructible" part of a mesh based on a
			  partial sampling of it: the set of triangles for which
			  at least one vertex is closer than a given geodesic
			  distance from at least one sample point.
		input:
			- mesh: the full mesh that has been sampled
			- pcd: the point cloud representing a sampling of in:mesh
			- alpha: the geodesic distance above which a triangle is
					 not appended to out:mesh_alpha
		output:
			- mesh_alpha: the subset of in:mesh for which triangles are
						  closer than in:alpha to in:pcd
	*/
	std::cout << "---> Computing mesh_alpha (reconstructible part of mesh)"
			  << std::endl;
	std::cout << "  max geodesic distance: " << alpha << std::endl;
	Mesh mesh_alpha;

	// vertices and faces of mesh already appended to mesh_alpha:
	std::vector<vertex_descriptor> mAlpha_v;
	std::vector<Mesh::Face_index> mAlpha_f;

	// correspondance map between vertex indices of mesh and mesh_alpha:
	// Vertex_map_meshes mapRefAlpha;
	std::map<vertex_descriptor, vertex_descriptor> mapRefAlpha;

	Surface_mesh_shortest_path shortest_paths = build_shortest_path(mesh, pcd);

	CGAL::Face_around_target_iterator<Mesh> face_b, face_e;

	for (vertex_descriptor vd : vertices(mesh)){
		Point p = mesh.point(vd);
		// std::cout << std::endl << vd << " : " << p << std::endl;

		// calculate face and barycentric coordinate of vd
		// and set as target
		double d = get<0>(shortest_paths.shortest_distance_to_source_points(vd));
		// std::cout << "distance to source points: " << d << std::endl;

		if (d < alpha && d >= 0){
		// if (d < alpha){
			// browse all faces incident to vertex vd:
			for (boost::tie(face_b, face_e) = CGAL::faces_around_target(mesh.halfedge(vd),mesh);
			face_b != face_e;
			++face_b)
			{
				if (face_b->is_valid()) // faces_around_target returns invalid faces for border vertices
				{
					// std::cout << " > " << *face_b << std::endl;
					// check if face_b is in mAlpha_f:
					if (std::find(mAlpha_f.begin(), mAlpha_f.end(), *face_b) == mAlpha_f.end()){
						// not already in mAlpha_f
						std::array<vertex_descriptor, 3> face_ref; // to store indices of face in mesh_ref
						int i = 0;
						// std::cout << "  - ";
						// browse all vertices adjacent to the face face_b:
						CGAL::Vertex_around_face_iterator<Mesh> vb,ve;
						for (boost::tie(vb,ve)=CGAL::vertices_around_face(mesh.halfedge(*face_b),mesh);
						vb != ve;
						++vb)
						{
							// std::cout << *vb << " ";
							// add vertex *vb to current face:
							face_ref[i] = *vb; i++;

							// check if vertex *vb is in mAlpha_v:
							if (std::find(mAlpha_v.begin(), mAlpha_v.end(), *vb) == mAlpha_v.end()){ // if not already in mAlpha_v
								mAlpha_v.push_back(*vb); // ... add it

								// add vertex *vb to mesh alpha:
								vertex_descriptor vb_alpha = mesh_alpha.add_vertex(mesh.point(*vb));

								// update map:
								mapRefAlpha.insert( std::pair<vertex_descriptor, vertex_descriptor>(*vb, vb_alpha) );
							}
						}
						// std::cout << std::endl;
						mAlpha_f.push_back(*face_b); // ... add face
						// Now, we are sure every vertices of *face_b are in mesh alpha
						// so we can add *face_b to M_alpha
						std::array<vertex_descriptor, 3> face_a = face_alpha(mapRefAlpha, face_ref);
						mesh_alpha.add_face(face_a);
					}
				}
			}
		}
	}
	std::cout << " Done with mesh_alpha computation" << std::endl << std::endl;
	return mesh_alpha;
}


int main(int argc, char** argv)
{
	std::cout << "--> program started <--" << std::endl << std::endl;

	// Mesh mesh;
	// std::ifstream is ("input_data/light/custom_mesh.off");
	// Mesh mesh; is >> mesh;
	// Mesh mesh = read_OBJ_mesh("input_data/heavy/open_data_strasbourg/PC3E45/PC3E45_3.obj");
	Mesh mesh = read_OBJ_mesh("input_data/heavy/cow.obj");

	std::cout << "test validity: " << mesh.is_valid(true) << std::endl;
	std::cout << "test whether it is closed: " << CGAL::is_closed(mesh) << std::endl;

	Point_set pcd;
	std::ifstream is_pcd ("output_data/out_cow.off");
	// CGAL::read_ply_point_set(is_pcd, pcd);
	// bool readSuccess = CGAL::read_off_points(is_pcd, )
	// CGAL::read_off_points(is_pcd, pcd.index_back_inserter());
	is_pcd >> pcd;


	int nbRmVert = CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
	std::cout << nbRmVert << " isolated vertices removed" << std::endl;

	Mesh mesh_alpha = Mesh_alpha(mesh, pcd, 0.65);

	std::ofstream of("output_data/mesh_alpha_cow.off");
	write_off(of,mesh_alpha);

	std::cout << "--> program ended <--" << std::endl;
	return 0;
}