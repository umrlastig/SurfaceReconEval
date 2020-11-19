#include <iostream>
#include <fstream>
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

// double* angles(double ang[], double start, double stop, int n){
//   for (int i=0; i<n; i++){
//     ang[i] = 23;
//   }
//   return ang;
// }
// void display(double tab[], int size){
//   std::cout << "elements of tab :" << std::endl;
//   for (int i=0; i<size; i++){
//     std::cout << tab[i] << std::endl;
//   }
//   std::cout << "=================" << std::endl;
// }

double angle(int k, double start, double stop, int n){
	return start + (2*k + 1)*(stop - start)/(2*n);
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tetrahedron.off";
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  Point_set point_set;
  std::ofstream outputFile("data/output_point_set.off");
  outputFile << "OFF" << std::endl;



  // std::cout << &mesh << std::endl;
  // int nTheta=10;
  // double a[nTheta];
  // a = angles(a,0,10,nTheta);
  // display(a, nTheta);
  // std::cout << a << std::endl;

  // Point b=mesh.vertices().begin();
  // for (Mesh::Vertex_index a : mesh.vertices() ){
  //   std::cout << "test" << std::endl;
  // }

  // Point A=CGAL::Polygon_mesh_processing::centroid(mesh);
  Point cent(0.25, 0.25, 0.25);
  std::cout << cent << std::endl;
  double radius = 10;

  double thetaMin=0; double thetaMax=M_PI; int nTheta = 100;
  double phiMin=0; double phiMax=2*M_PI; int nPhi = 200;
  int nPts = nPhi * nTheta;
  outputFile << nPts << " 0 0" << std::endl; // nb of vertices, faces and edges

  for (int kTheta=0; kTheta<nTheta; kTheta++){
    double theta = angle(kTheta, thetaMin, thetaMax, nTheta);
    for (int kPhi=0; kPhi<nPhi; kPhi++){
    	double phi = angle(kPhi, phiMin, phiMax, nPhi);
    	double x = sin(theta) * cos(phi);
    	double y = sin(theta) * sin(phi);
    	double z = cos(theta);
    	Vector dir = Vector(x,y,z);
    	Ray ray(cent, dir);
    	Ray_intersection intersection = tree.first_intersection(ray);

    	if(intersection){
    		if(boost::get<Point>(&(intersection->first))){
    			const Point* p =  boost::get<Point>(&(intersection->first) );
    			point_set.insert (*p);
    			outputFile << *p << std::endl;
    			// std::cout <<  *p << std::endl;
    		}
    	}
    }
  }


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



  std::cerr << "done" << std::endl;
  return 0;
}