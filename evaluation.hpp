#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>

// Input / Output
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/boost/graph/io.h>

#include <CGAL/boost/graph/iterator.h>

// Polygon_mesh_processing
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

// neighbor searching
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <list>
#include <cmath>


// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Point_set_3<Point> Point_set;

// shortest path
typedef CGAL::Surface_mesh_shortest_path_traits<K, Mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> SM_shortest_path;

// iterators and descriptors
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

// accelerating data structure
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Tree_traits;
typedef CGAL::AABB_tree<Tree_traits> Tree;

// neighbor searching
typedef CGAL::Search_traits_3<K> TreeTraitsNS;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraitsNS> Neighbor_search;
typedef Neighbor_search::Tree TreeNS;

namespace PMP = CGAL::Polygon_mesh_processing;

std::string getFileExt(const char* fileName) {
	const std::string s(fileName);
	size_t i = s.rfind('.', s.length());
	if (i != std::string::npos) {
		return(s.substr(i+1, s.length() - i));
	}
	return("");
}

Mesh read_OBJ_mesh(const char* fileName, bool verbose){
	Mesh mesh;

	// Variables to store polygon soup:
	std::ifstream is (fileName);
	std::vector<Point> points;
	std::vector<std::vector<std::size_t>> polygons;

	if (!CGAL::read_OBJ(is, points, polygons)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
	if (verbose) std::cout << "number of points: " << points.size() << std::endl;
	if (verbose) std::cout << "number of polygons: " << polygons.size() << std::endl;

	// Processing polygon soup:
	if (!PMP::is_polygon_soup_a_polygon_mesh(polygons))
	{
		// polygons do not define a valid mesh
		if (verbose) std::cout << "Warning: [" << fileName << "] does not define a valid polygon mesh " << std::endl;
		if ( !PMP::orient_polygon_soup(points, polygons) ){ // try to orient
			// orientation failed
			if (verbose) std::cout << "Warning: the orientation of the polygon soup [" << fileName << "] failed. => Points were duplicated." << std::endl;
			if (verbose) std::cout << "Number of points: " << mesh.number_of_vertices() << std::endl;
			if (verbose) std::cout << "Number of faces: " << mesh.number_of_faces() << std::endl;
		} else {
			// orientation succeeded
			if (verbose) std::cout << "orientation of the polygon soup [" << fileName << "] was a success." << std::endl;
		}
	} else {
		// polygons define a valid mesh
		if (verbose) std::cout << "Success: [" << fileName << "] defines a valid polygon mesh " << std::endl;
	}

	PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

	return mesh;
}

Mesh readMesh(const char* fileName, bool verbose){
	if (verbose) std::cout << "---> READING MESH: " << fileName << std::endl;
	Mesh mesh; std::ifstream is (fileName);

	// Determine extension then extract points and polygons
	std::string extension = getFileExt(fileName);
	if (extension == "obj"){
		mesh = read_OBJ_mesh(fileName, verbose);
	} else if (extension == "ply") {
		if (!read_ply(is, mesh)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
		// std::vector<Point> points;
		// std::vector<std::vector<std::size_t>> polygons;
		// std::vector<CGAL::Color> fcolors;
		// std::vector<CGAL::Color> vcolors;
		// bool success = CGAL::read_PLY (is, points, polygons, fcolors, vcolors);
		// if (success) std::cout << "success" << std::endl;
		// std::cout << "number of points: " << points.size() << std::endl;
		// std::cout << "number of polygons: " << polygons.size() << std::endl;
		// PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
	} else if (extension == "off") {
		if (!read_off(is, mesh)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
	} else {
		std::cout << "\nERROR: extension [." << extension << "] not readable.\n" << std::endl;
	}
	if (verbose) std::cout << "Number of points: " << mesh.number_of_vertices() << std::endl;
	if (verbose) std::cout << "Number of faces: " << mesh.number_of_faces() << std::endl;
	return mesh;
}

Point_set readPointSet(const char* fileName, bool verbose){
	if (verbose) std::cout << "---> READING POINT SET: " << fileName << std::endl;
	Point_set pcd; std::ifstream is (fileName);

	std::string extension = getFileExt(fileName);
	if (extension == "ply") {
		if (!CGAL::read_ply_point_set(is, pcd)) std::cerr
			<< "\nERROR: impossible to read point set\n" << std::endl;
	} else if (extension == "off") {
		if (!CGAL::read_off_point_set(is, pcd)) std::cerr
			<< "\nERROR: impossible to read mesh\n" << std::endl;
	} else {
		std::cerr << "\nERROR: extension [." << extension << "] not readable.\n" << std::endl;
	}
	if (verbose) std::cout << "Number of points: " << pcd.number_of_points() << std::endl;
	return pcd;
}

void writeMesh(const char* fileName, Mesh &mesh, bool verbose){
	std::ofstream of(fileName);
	bool success = false;

	// Determine extension then use appropriate function
	std::string extension = getFileExt(fileName);
	if (extension == "ply"){
		success = write_ply(of, mesh);
	} else if (extension == "off") {
		success = write_off(of, mesh);
	} else {
		std::cerr << "ERROR: extension [." << extension << "] not writable." << std::endl;
	}
	// Outcome:
	if (success){
		if (verbose) std::cout << "wrote: " << fileName << std::endl;
	} else {
		std::cerr << "ERROR: impossible to write '" << fileName << "'" << std::endl;
	}
}

void writePointSet(const char* fileName, Point_set &pcd, bool verbose){
	std::ofstream of(fileName);
	bool success = false;

	// Determine extension then use appropriate function
	std::string extension = getFileExt(fileName);
	if (extension == "ply"){
		success = CGAL::write_ply_point_set(of, pcd);
	} else if (extension == "off") {
		success = CGAL::write_off_point_set(of, pcd);
	} else {
		std::cerr << "ERROR: extension [." << extension << "] not writable." << std::endl;
	}
	// Outcome:
	if (success){
		if (verbose) std::cout << "wrote: " << fileName << std::endl;
	} else {
		std::cerr << "ERROR: impossible to write '" << fileName << "'" << std::endl;
	}
}

void build_shortest_path(Point_set &pcd, SM_shortest_path &shortest_paths, bool verbose){
	if (verbose) std::cout << "---> Building Shortest Path Sequence Tree" << std::endl;
	// Tree tree(faces(mesh).first, faces(mesh).second, mesh); // AABB_tree of the mesh
	Tree tree;
	shortest_paths.build_aabb_tree(tree); // Creates an AABB_tree suitable for use with locate
	// browse point set:
	for (Point_set::const_iterator pi = pcd.begin(); pi != pcd.end(); ++pi){
		// add point *pi to source:
		shortest_paths.add_source_point( shortest_paths.locate(pcd.point(*pi), tree) );
	}
	if (verbose) std::cout << " Number of source points: " <<
		shortest_paths.number_of_source_points() << std::endl;

	shortest_paths.build_sequence_tree();

	if (verbose) std::cout << " Done with building sequence tree" << std::endl;
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

void build_search_structure(bool geodesic,
							Point_set &pcd,
							SM_shortest_path &shortest_paths,
							TreeNS &tree,
							bool verbose ){	
	if (geodesic) // using geodesic distance
	{
		build_shortest_path(pcd, shortest_paths, verbose);
	} else { // using euclidean distance
		for (Point_set::iterator it = pcd.begin(); it != pcd.end(); ++it){
			tree.insert(pcd.point(*it));
		}
	}
}

double geod_or_eucli_distance(bool geodesic,
							  SM_shortest_path &shortest_paths,
							  vertex_descriptor &vd,
							  Point &p,
							  TreeNS &tree
							  ){
	double dist = 0;
	if (geodesic) // geodesic distance
	{
		dist = get<0>(shortest_paths.shortest_distance_to_source_points(vd));
	} else {  // euclidean distance
		Neighbor_search search(tree, p, 1);
		Neighbor_search::iterator itS = search.begin(); // only one neighbor
		dist = std::sqrt(itS->second);
	}
	return dist;
}

Mesh Mesh_alpha(Mesh &mesh, Point_set &pcd, double alpha, bool geodesic, bool verbose, bool debug){
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
	if (verbose) std::cout << "\n---> Computing mesh_alpha (reconstructible part of mesh)"
			  << std::endl;
	if (geodesic) {
		if (verbose) std::cout << "     max geodesic distance: " << alpha << std::endl;
	} else {
		if (verbose) std::cout << "     max euclidean distance: " << alpha << std::endl;
	}
	
	Mesh meshAlpha;

	// vertices and faces of mesh already appended to mesh_alpha:
	std::vector<vertex_descriptor> mAlpha_v;
	std::vector<Mesh::Face_index> mAlpha_f;

	// correspondance map between vertex indices of mesh and mesh_alpha:
	std::map<vertex_descriptor, vertex_descriptor> mapRefAlpha;

	// search structure:
	TreeNS tree; // euclidean distance
	SM_shortest_path shortest_paths(mesh); // geodesic distance
	build_search_structure(geodesic, pcd, shortest_paths, tree, verbose);

	for (vertex_descriptor vd : vertices(mesh)){
		Point p = mesh.point(vd);
		if (debug) std::cout << std::endl << vd << ": [" << p << "]" << std::endl;
		double d = geod_or_eucli_distance(geodesic, shortest_paths, vd, p, tree);
		if (debug) std::cout << "distance to source points: " << d << std::endl;

		if (d < alpha && d >= 0){
		// if (d < alpha){
			// browse all faces incident to vertex vd:
			CGAL::Face_around_target_iterator<Mesh> face_b, face_e;
			for (boost::tie(face_b, face_e) = CGAL::faces_around_target(mesh.halfedge(vd),mesh);
			face_b != face_e;
			++face_b)
			{
				if (face_b->is_valid()) // faces_around_target returns invalid faces for border vertices
				{
					if (debug) std::cout << " > " << *face_b << std::endl;
					// check if face_b is in mAlpha_f:
					if (std::find(mAlpha_f.begin(), mAlpha_f.end(), *face_b) == mAlpha_f.end()){
						// not already in mAlpha_f
						std::array<vertex_descriptor, 3> face_ref; // to store indices of face in mesh_ref
						int i = 0;
						if (debug) std::cout << "  - ";
						// browse all vertices adjacent to the face face_b:
						CGAL::Vertex_around_face_iterator<Mesh> vb,ve;
						for (boost::tie(vb,ve)=CGAL::vertices_around_face(mesh.halfedge(*face_b),mesh);
						vb != ve;
						++vb)
						{
							if (debug) std::cout << *vb << " ";
							// add vertex *vb to current face:
							face_ref[i] = *vb; i++;

							// check if vertex *vb is in mAlpha_v:
							if (std::find(mAlpha_v.begin(), mAlpha_v.end(), *vb) == mAlpha_v.end()){ // if not already in mAlpha_v
								mAlpha_v.push_back(*vb); // ... add it

								// add vertex *vb to mesh alpha:
								vertex_descriptor vb_alpha = meshAlpha.add_vertex(mesh.point(*vb));

								// update map:
								mapRefAlpha.insert( std::pair<vertex_descriptor, vertex_descriptor>(*vb, vb_alpha) );
							}
						}
						if (debug) std::cout << std::endl;
						mAlpha_f.push_back(*face_b); // ... add face
						// Now, we are sure every vertices of *face_b are in mesh alpha
						// so we can add *face_b to M_alpha
						std::array<vertex_descriptor, 3> face_a = face_alpha(mapRefAlpha, face_ref);
						meshAlpha.add_face(face_a);
					}
				}
			}
		}
	}
	if (verbose) std::cout << "     Done with mesh_alpha computation" << std::endl;
	return meshAlpha;
}

void remove_points_too_far_from_P(Point_set &pcd, Point_set &P, double alpha, bool verbose){
	/*
		usage:
			- Remove points from in:pcd that are farther from P than in:alpha
		input:
			- pcd: a point cloud representing a sampling of a surface
			- P: the input point cloud used for a reconstruction
			- alpha: the distance above which a point is deleted
	*/

	if (verbose) std::cout << "Started removing points farther than alpha = "
		<< alpha << " from P" << std::endl;

	TreeNS tree(P.points().begin(), P.points().end());
	std::vector<Point_set::const_iterator> ptsToRemove;

	for (Point_set::const_iterator pi = pcd.begin(); pi != pcd.end(); ++pi){
		Point query = pcd.point(*pi);
		Neighbor_search search(tree, query, 1);
		Neighbor_search::iterator itS = search.begin(); // only one neighbor

		// to display point and its distance to nearest neighbor:
		// std::cout << itS->first << " "<< std::sqrt(itS->second) << std::endl;

		double dist = std::sqrt(itS->second);
		if (dist > alpha)
		{
			ptsToRemove.push_back(pi); // store point to remove
		}
	}
	// remove all points too far:
	for (int i = 0; i != ptsToRemove.size(); ++i){
		pcd.remove(*ptsToRemove[i]);
	}
	if (verbose) std::cout << "Done: " << ptsToRemove.size() << " point(s) removed" << std::endl;
}

void mean_and_max_distance_from_P_to_mesh(Point_set &pcd, Mesh &mesh){
	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	double mean_distance = 0;
	double max_distance = 0;
	for (Point_set::const_iterator pi = pcd.begin(); pi != pcd.end(); ++pi){
		Point query = pcd.point(*pi);
		// Point closest = tree.closest_point(query); // optional
		double distance = std::sqrt( tree.squared_distance(query) );

		// std::cout << "query: " << query << std::endl;
		// std::cout << "closest: " << closest << std::endl;
		// std::cout << "distance: " << distance << std::endl << std::endl;

		// update distances
		mean_distance += distance;
		max_distance = (distance > max_distance) ? distance : max_distance;
	}
	int nbPts = pcd.end() - pcd.begin();
	mean_distance /= nbPts;

	// CGAL built-in function for max distance computation:
	// double max = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(pcd.points(), mesh);

	std::cout << "mean distance: " << mean_distance << std::endl;
	std::cout << "max distance: " << max_distance << std::endl;
}

