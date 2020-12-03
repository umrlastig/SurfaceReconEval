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
#include <CGAL/Point_set_3/IO.h>
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
// typedef boost::optional< Primitive_id > Primitive_id;
typedef Tree::Primitive_id Primitive_id;


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
	/*
		usage:
			Shoot rays in all directions from a fixed position in space
			and compute intersections with a mesh
		input:
			Mesh mesh: mesh to be raycasted/sampled
			Point origin: The position from which the rays must be shooted
			int nTheta: nb of polar angles in [0,pi] interval to generate
			int nPhi: nb of azimuthal angles in [0,2pi] interval to generate
		output:
			Point_set point_set: intersections of the rays with the mesh

	*/
	std::cout << "Spherical Ray Shooting from position (" << origin << ")" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	Point_set point_set; // for rayshooting output storage
	point_set.add_normal_map(); // add normal property
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
					const Point& p = boost::get<Point>(intersection->first);					
					const face_descriptor f = boost::get<face_descriptor>(intersection->second);
					const Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f,mesh);
					// point_set.insert (p);
					point_set.insert(p,n); // A normal property must have been added to the point set before using this method
				}
			}
		}
	}
	std::cout << "Done with Spherical Ray Shooting" << std::endl << std::endl;
	return point_set;
}

Vector normalize_vector_3(Vector &v){
	double norm_v = sqrt(v.squared_length());
	return v.operator/(norm_v);
}

void compute_local_frame(Vector &AB, Vector &vec_i, Vector &vec_j, Vector &vec_k){
	std::cout << "Computing local frame" << std::endl;
	vec_k = normalize_vector_3(AB);

	vec_j = CGAL::cross_product(Vector(0,0,1), vec_k);
	vec_j = normalize_vector_3(vec_j);

	vec_i = CGAL::cross_product(vec_j, vec_k);

	std::cout << "vec_k : [" << vec_k << "]" << std::endl;
	std::cout << "vec_i : [" << vec_i << "]" << std::endl;
	std::cout << "vec_j : [" << vec_j << "]" << std::endl << std::endl;
}

void compute_and_add_intersection(Mesh &mesh, Tree &tree, Ray &ray, Point_set &point_set){
	Ray_intersection intersection = tree.first_intersection(ray);
	if(intersection){
		if(boost::get<Point>(&(intersection->first))){
			const Point& p = boost::get<Point>(intersection->first);					
			const face_descriptor f = boost::get<face_descriptor>(intersection->second);
			const Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f,mesh);
			// point_set.insert (p);
			point_set.insert(p,n); // A normal property must have been added to the point set before using this method
		}
	}
}

Point_set aerial_lidar(Mesh &mesh, Point &A, Point &B, double v0, double omega, double theta_0, double freq){
	/*
		usage:
			Simulate an aerial LiDAR acquisition
		input:
			Mesh mesh: mesh to be raycasted/sampled
			Point A: starting point of the vehicle
			Point B: ending point of the vehicle
			double v0: speed of the vehicle
			double omega: angular speed of the laser pointer
			double theta_0: initial angle of the laser pointer
			double freq: frequency (time resolution) at which
						 the analysis should be carried out
		output:
			Point_set point_set: intersections of the rays with the mesh
	*/
	std::cout << "---> Aerial LiDAR acquisition from A(" << A << ") to B(" << B << ")" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	Point_set point_set; // for rayshooting output storage
	point_set.add_normal_map(); // add normal property

	Vector AB = Vector(A,B);
	Vector vec_i; Vector vec_j; Vector vec_k;
 	compute_local_frame(AB, vec_i, vec_j, vec_k);

 	Point M; Vector di; Vector dj; Vector ML; double theta; // Some necesary variables

 	double alpha = 90 * M_PI / 180; // max semi-angle for a ray being actually shooted
 	double tB = sqrt(AB.squared_length()) / v0; // duration of simulation
 	int nT = (int) (freq * tB); // number of time steps
 	double dt = tB / (nT-1); // infinitesimal unit of time
 	double t = 0; // initial time 	

	for (int ti=0; ti<nT; ti++){
		theta = omega*t + theta_0; // angular position
		theta = std::fmod(theta, 2*M_PI); // euclidean division
		if (theta < alpha || theta > (2*M_PI - alpha)){ // theta is eligible
			M = operator+(A, operator*(v0*t, vec_k)); // position of M
			di = operator*(cos(theta), vec_i);
			dj = operator*(sin(theta), vec_j);
			ML = di.operator+(dj); // laser pointer
			Ray ray(M, ML);
			compute_and_add_intersection(mesh, tree, ray, point_set);
		} else {
			// do nothing
		}
		t += dt; // update time
	}
	std::cout << "Done with Aerial LiDAR acquisition" << std::endl << std::endl;
	return point_set;
}

Mesh read_OFF_mesh(const char* fileName){
	std::ifstream is (fileName);
	Mesh mesh;
	is >> mesh;
	// read_off(is, mesh);
	return mesh;
}

Mesh read_OBJ_mesh(const char* fileName){
	std::cout << "---> READING : " << fileName << std::endl;
	Mesh mesh;
	// Variables to store polygon soup :
	std::ifstream is (fileName);
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

void Strasbourg_Scene_Raytracing(int nTheta, int nPhi){
	/*
		Routine to raytrace the model PC3E45_3.obj from central position in x/y
		and z = 3*bbox (above the model) and write the results to a ply pcd
		with normals
	*/
	const char* filename = "input_data/heavy/open_data_strasbourg/PC3E45/PC3E45_3.obj";
	Mesh mesh = read_input_file(filename);
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	double x_origin = ( bbox.xmin() + bbox.xmax() ) / 2;
	double y_origin = ( bbox.ymin() + bbox.ymax() ) / 2;
	double z_origin = bbox.zmin() + 3*( bbox.zmax() - bbox.zmin() );
	Point origin_lidar(x_origin, y_origin, z_origin);
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi);
	std::ofstream of("output_data/out_pcd_normal_PC3E45_3.ply");
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/out_pcd_normal_PC3E45_3.ply" << std::endl;
}

void Strasbourg_Scene_Aerial_Lidar(double v0, double omega, double theta_0, double freq){
	/*
		Routine to simulate aerial LiDAR on the model PC3E45_3.obj from central 
		position in x, min and max y and z = 3*bbox (above the model) and write
		the results to a ply pcd with normals
	*/
	const char* filename = "input_data/heavy/open_data_strasbourg/PC3E45/PC3E45_3.obj";
	Mesh mesh = read_input_file(filename);
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	double x_mid = ( bbox.xmin() + bbox.xmax() ) / 2;
	double y_A = bbox.ymin();
	double y_B =  bbox.ymax();
	double z_above = bbox.zmin() + 3*( bbox.zmax() - bbox.zmin() );

	Point A(x_mid, y_A, z_above);
	Point B(x_mid, y_B, z_above);

	Point_set point_set = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq);
	std::ofstream of("output_data/out_pcd_normal_AERIAL_PC3E45_3.ply");
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/out_pcd_normal_AERIAL_PC3E45_3.ply" << std::endl;
}

void object_raytracing_from_centroid(const char* filename, int nTheta, int nPhi){
	/*
		Routine to raytrace any closed object from its centroid and
		write the results to a ply pcd with normals
	*/
	Mesh mesh = read_input_file(filename);
	Point origin_lidar(centroid(mesh));
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi);
	std::ofstream of("output_data/out_pcd_normal_00.off");
	CGAL::write_off_point_set(of,point_set);
	std::cerr << "wrote output_data/out_pcd_normal_00.off" << std::endl;

	// just change the name of the file and user this function for a PLY pcd :
	// CGAL::write_ply_point_set(ofPLY,point_set);
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "input_data/light/cube.obj"; 
  // object_raytracing_from_centroid(filename, 150, 300);

  // Strasbourg_Scene_Raytracing(150, 300);
  double v0 = 50;
  double omega = 200*M_PI;
  double freq = 300000;
  double theta_0 = 0;
  Strasbourg_Scene_Aerial_Lidar(v0, omega, theta_0, freq);

  std::cerr << "done" << std::endl;
  return 0;
}