#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/boost/graph/io.h>
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Mesh::Face_range Face_range;

struct Skip {
  face_descriptor fd;
  Skip(const face_descriptor fd)
    : fd(fd)
  {}
  bool operator()(const face_descriptor& t) const
  { if(t == fd){
      std::cerr << "ignore " << t  <<std::endl;
    };
    return(t == fd);
  }
};


double angle(int k, double start, double stop, int n){
	return start + (2*k + 1)*(stop - start)/(2*n);
}

void print_point_set(const Point_set& point_set){
	std::cerr << "Content of point set :" << std::endl;
	for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); it++){
		Point P = point_set.point(*it);
		std::cerr << "Point " << *it << " : " << P << std::endl;
	}
}

void write_point_set_to_file(const Point_set& point_set, std::string file_name){
	std::ofstream outputFile("output_data/" + file_name);
	outputFile << "OFF" << std::endl; // header
	std::size_t nPts = point_set.size();
	outputFile << nPts << " 0 0" << std::endl; // nb of vertices, faces and edges
	for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); it++){
		Point P = point_set.point(*it);
		outputFile << P << std::endl;
	}
	std::cout << "wrote : " << file_name << std::endl << std::endl;
}

Point centroid(const Mesh mesh){
	double areaM = 0; // total area of mesh
	Point cM = Point(0,0,0); // centroid of mesh
	for( face_descriptor fd : faces(mesh) ){
		// area of face fd :
		double areaT = CGAL::Polygon_mesh_processing::face_area(fd,mesh);

		// centroid of face fd :
		halfedge_descriptor hd = halfedge(fd,mesh);
		Point cT = CGAL::centroid(mesh.point(source(hd,mesh)),
								 mesh.point(target(hd,mesh)),
								 mesh.point(target(next(hd,mesh),mesh)));

		// update :
		areaM += areaT;
		Vector v = Vector(Point(0,0,0), cT).operator*=(areaT);
		cM.operator+=(v);
	}
	return Point( cM.x(), cM.y(), cM.z(), areaM );
}

Point_set spherical_ray_shooting(Mesh mesh, Point origin, int nTheta, int nPhi){
	std::cout << "Spherical Ray Shooting from position (" << origin << ")" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	Point_set point_set; // for rayshooting output storage
	double thetaMin=0; double thetaMax=M_PI; // polar angle theta in [0,pi]
	double phiMin=0; double phiMax=2*M_PI; // azimuthal angle phi in [0,2pi]
	for (int kTheta=0; kTheta<nTheta; kTheta++){
		double theta = angle(kTheta, thetaMin, thetaMax, nTheta);
		for (int kPhi=0; kPhi<nPhi; kPhi++){
			double phi = angle(kPhi, phiMin, phiMax, nPhi);
			double x = sin(theta) * cos(phi);
			double y = sin(theta) * sin(phi);
			double z = cos(theta);
			Vector dir = Vector(x,y,z);
			Ray ray(origin, dir); // ray shooted

			Ray_intersection intersection = tree.first_intersection(ray);
			if(intersection){
				if(boost::get<Point>(&(intersection->first))){
					const Point* p =  boost::get<Point>(&(intersection->first) );
					point_set.insert (*p);
				}
			}
		}
	}
	std::cout << "Done with Spherical Ray Shooting" << std::endl << std::endl;
	return point_set;
}

Mesh read_OFF_mesh(const char* fileName){
	std::ifstream is (fileName, std::ifstream::in);
	Mesh mesh;
	is >> mesh;
	// read_off(is, mesh);
	return mesh;
}

Mesh read_OBJ_mesh(const char* fileName){
	std::cout << "READING : " << fileName << std::endl;
	Mesh mesh;
	// Variables to store polygon soup :
	std::ifstream is (fileName, std::ifstream::in);
	std::vector<Point> points;
	std::vector<std::vector<std::size_t>> polygons;

	// Extract points and polygons from stream 'is' :
	CGAL::read_OBJ(is, points, polygons);	
	std::cout << "number of points : " << points.size() << std::endl;
	std::cout << "number of polygons : " << polygons.size() << std::endl;

	// Processing polygon soup :
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

	std::cout << "Number of points : " << mesh.number_of_vertices() << std::endl;
	std::cout << "Number of faces : " << mesh.number_of_faces() << std::endl;
	std::cout << "Done with reading" << std::endl << std::endl;

	return mesh;
}

Mesh read_input_file(const char* fileName){
	/*
	Use the proper function absed on the extension of the input file
	*/
	Mesh mesh;
	std::string extension = std::string(fileName).substr(std::string(fileName).find_last_of(".") + 1);
	if (extension == "obj"){
		mesh = read_OBJ_mesh(fileName);
	} else if (extension == "off") {
		mesh = read_OFF_mesh(fileName);
	} else {
		std::cout << "Error: extension [." << extension << "] not readable." << std::endl;
	}
	return mesh;
}

int main(int argc, char* argv[])
{
  // const char* filename = (argc > 1) ? argv[1] : "input_data/tetrahedron.off";
  // const char* filename = (argc > 1) ? argv[1] : "input_data/crocodile_statue.ply";
  // const char* filename = (argc > 1) ? argv[1] : "input_data/open_data_strasbourg/PC3E45/PC3E45_3.obj";
  //const char* filename = (argc > 1) ? argv[1] : "input_data/cube.obj";
  const char* filename = (argc > 1) ? argv[1] : "input_data/cow.obj";
  Mesh mesh = read_input_file(filename);
  CGAL::write_off("input_data/heavy/open_data_strasbourg/PC3E45/without_texture.off", mesh);
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  double x_origin = ( bbox.xmin() + bbox.xmax() ) / 2;
  double y_origin = ( bbox.ymin() + bbox.ymax() ) / 2;
  double z_origin = bbox.zmin() + 3*( bbox.zmax() - bbox.zmin() );
  Point origin_lidar(x_origin, y_origin, z_origin);

  // Point origin_lidar(centroid(mesh));
  // Point origin_lidar(-2, 2, 2);
  Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, 150, 300);
  write_point_set_to_file(point_set, "out_pcd.off");

  std::cerr << "done" << std::endl;
  return 0;
}