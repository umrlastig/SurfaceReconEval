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
#include <CGAL/Point_set_3.h>
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
}

void test(){
	std::cerr << "beginning of test : " << std::endl;
	double d(4.5);
	Point P1(4,10,10);
	Point P2(5,18,-4);
	Point P3(0,14,5);
	// Point c = CGAL::centroid(P1,P2,P3);
	// std::cerr << "centroid : " << P << std::endl;
	std::cerr << "end of test : " << std::endl;
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
	return point_set;
}

Mesh read_OFF_mesh(const char* fileName){
	std::ifstream input(fileName);
	Mesh mesh;
	input >> mesh;
	return mesh;
}

Mesh read_OFF_mesh2(const char* fileName){
	Mesh mesh;
	std::filebuf fb; // Stream buffer to read from and write to files
	if (fb.open (fileName, std::ios::in)){
		std::istream is(&fb);
		// is >> mesh;
		read_off(is, mesh);
	} else {
		std::cout << "Error: impossible to read " << fileName << std::endl;
	}
	return mesh;
}

// Mesh read_PLY_mesh(const char* fileName){
// 	Mesh mesh;
// 	std::filebuf fb; // Stream buffer to read from and write to files
// 	if (fb.open (fileName, std::ios::in)){
// 		std::istream is(&fb);
// 		// is >> mesh;
// 		// CGAL::read_ply(is, mesh);
// 	} else {
// 		std::cout << "Error: impossible to read " << fileName << std::endl;
// 	}
// 	return mesh;
// }

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "input_data/tetrahedron.off";
  // const char* filename = (argc > 1) ? argv[1] : "input_data/crocodile_statue.ply";
  Mesh mesh = read_OFF_mesh(filename);
  Point cent(0.25, 0.25, 0.25);
  Point_set point_set = spherical_ray_shooting(mesh, cent, 100, 200);
  write_point_set_to_file(point_set, "out_pcd.off");

  Point A(centroid(mesh));
  std::cout << "centroid A : " << A << std::endl;
  
  // Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  // double d = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)?-1:1;
  // for(face_descriptor fd : faces(mesh)){
  //   halfedge_descriptor hd = halfedge(fd,mesh);
  //   Point p = CGAL::centroid(mesh.point(source(hd,mesh)),
  //                            mesh.point(target(hd,mesh)),
  //                            mesh.point(target(next(hd,mesh),mesh)));
  //   Vector v = CGAL::Polygon_mesh_processing::compute_face_normal(fd,mesh);
  //   Ray ray(p,d * v);
  //   Skip skip(fd);
  //   Ray_intersection intersection = tree.first_intersection(ray, skip);
  //   if(intersection){
  //     if(boost::get<Point>(&(intersection->first))){
  //       const Point* p =  boost::get<Point>(&(intersection->first) );
  //       std::cout <<  *p << std::endl;
  //     }
  //   }
  // }

  // test();

  std::cerr << "done" << std::endl;
  return 0;
}