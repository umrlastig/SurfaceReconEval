#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Point_set_3<Point> Point_set;

// shortest path
typedef CGAL::Surface_mesh_shortest_path_traits<K, Mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;

// iterators
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;


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

void test(){

	std::ifstream is ("output_data/out_origin.ply");
	Point_set pcd;
	// X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	// bool success = false;
	// boost::tie (x_origin, success) = pcd.add_property_map<double>("x_origin", 0);
	// boost::tie (y_origin, success) = pcd.add_property_map<double>("y_origin", 0);
	// boost::tie (z_origin, success) = pcd.add_property_map<double>("z_origin", 0);
	// assert(success);
	is >> pcd;
	Mesh mesh = read_OBJ_mesh("input_data/light/cube.obj");

	// vertices and faces of mesh already appended to mesh_alpha:
	std::vector<vertex_descriptor> mAlpha_v;
	std::vector<Mesh::Face_index> mAlpha_f;

	Surface_mesh_shortest_path shortest_paths(mesh);
	// write a for loop which browse Point Cloud pcd
	// for each point, compute its face and barycentric coordinates
	// using Surface_mesh_shortest_path::locate() (see doc to avoid
	// rebuilding the tree each time). By the way it is probably
	// useful to write a function which takes as input pcd and mesh
	// and returns the set of pairs (face, barycentric coord) for
	// each point. It does that once for all (so that we don't have
	// to recompute everything after changing lambda).

	CGAL::Face_around_target_iterator<Mesh> face_b, face_e;

	for (vertex_descriptor vd : vertices(mesh)){
		Point p = mesh.point(vd);
		std::cout << std::endl << vd << " : " << p << std::endl;

		// calculate face and barycentric coordinate of vd
		// and set as target

		//check if p is in mAlpha_v:
		if (std::find(mAlpha_v.begin(), mAlpha_v.end(), vd) == mAlpha_v.end()){ // if not already in mAlpha_v
			mAlpha_v.push_back(vd); // ... add it
		}
		// browse all faces incident to vertex vd:
		for (boost::tie( face_b, face_e) = CGAL::faces_around_target(mesh.halfedge(vd),mesh);
		face_b != face_e;
		++face_b)
		{
			std::cout << " > " << *face_b << std::endl;

			// check if face_b is in mAlpha_f:
			if (std::find(mAlpha_f.begin(), mAlpha_f.end(), *face_b) == mAlpha_f.end()){ // if not already in mAlpha_f
				std::cout << "  - ";
				// browse all vertices adjacent to the face face_b:
				CGAL::Vertex_around_face_iterator<Mesh> vb,ve;
				for (boost::tie(vb,ve)=CGAL::vertices_around_face(mesh.halfedge(*face_b),mesh);
				vb != ve;
				++vb)
				{
					std::cout << *vb << " ";
					// check if vertex *vb is in mAlpha_v:
					if (std::find(mAlpha_v.begin(), mAlpha_v.end(), *vb) == mAlpha_v.end()){ // if not already in mAlpha_v
						mAlpha_v.push_back(*vb); // ... add it
						// add vertex *vb to mesh alpha
					}
				}
				std::cout << std::endl;
				mAlpha_f.push_back(*face_b); // ... add face
				// Now, we are sure every vertices of *face_b are in mesh alpha
				// so we can add *face_b to M_alpha
			}

			


		}

	}
	std::cout << std::endl;
	// display added vertices:
	for (std::vector<vertex_descriptor>::const_iterator iii=mAlpha_v.begin(); iii != mAlpha_v.end(); ++iii){
		std::cout << *iii << " - ";
	}
	std::cout << std::endl;
	// display added faces:
	for (std::vector<Mesh::Face_index>::const_iterator iii=mAlpha_f.begin(); iii != mAlpha_f.end(); ++iii){
		std::cout << *iii << " - ";
	}
	std::cout << std::endl;




	// std::ofstream of("output_data/noisy_pcd.ply");

	// CGAL::write_ply_point_set(of,pcd);

}



int main(int argc, char** argv)
{

	test();

	std::cout << "done" << std::endl;
	return 0;
}