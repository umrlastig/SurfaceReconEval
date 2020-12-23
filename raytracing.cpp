#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <random>
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
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/tags.h>

#include <CGAL/jet_estimate_normals.h>

// kernel
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;

// geometrical objects
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Point_set_3<Point> Point_set;

// iterators and descriptors
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

// accelerating data structure
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

// optical centers
typedef Point_set::Property_map<double> X_Origin_Map;
typedef Point_set::Property_map<double> Y_Origin_Map;
typedef Point_set::Property_map<double> Z_Origin_Map;



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
	std::cerr << "Content of point set:" << std::endl;
	for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); it++){
		Point P = point_set.point(*it);
		std::cerr << "Point " << *it << " : " << P << std::endl;
	}
}

Point centroid(const Mesh mesh){
	double areaM = 0; // total area of mesh
	Point cM = Point(0,0,0); // centroid of mesh
	for( face_descriptor fd : faces(mesh) ){
		// area of face fd:
		double areaT = CGAL::Polygon_mesh_processing::face_area(fd,mesh);

		// centroid of face fd:
		halfedge_descriptor hd = halfedge(fd,mesh);
		Point cT = CGAL::centroid(mesh.point(source(hd,mesh)),
								 mesh.point(target(hd,mesh)),
								 mesh.point(target(next(hd,mesh),mesh)));

		// update:
		areaM += areaT;
		Vector v = Vector(Point(0,0,0), cT).operator*=(areaT);
		cM.operator+=(v);
	}
	return Point( cM.x(), cM.y(), cM.z(), areaM );
}

Vector normalize_vector_3(Vector &v){
	double norm_v = sqrt(v.squared_length());
	return v.operator/(norm_v);
}

void compute_local_frame(Vector &AB, Vector &vec_i, Vector &vec_j, Vector &vec_k){
	std::cout << " > Computing local frame" << std::endl;
	vec_k = normalize_vector_3(AB);

	vec_j = CGAL::cross_product(Vector(0,0,1), vec_k);
	vec_j = normalize_vector_3(vec_j);

	vec_i = CGAL::cross_product(vec_j, vec_k);

	std::cout << "  vec_k: [" << vec_k << "]" << std::endl;
	std::cout << "  vec_i: [" << vec_i << "]" << std::endl;
	std::cout << "  vec_j: [" << vec_j << "]" << std::endl;
}

void compute_and_orient_normals(Point_set &pcd, unsigned int k){
	/*
		Usage:
			- compute and orient normals of a point set containing optical centers
				o normal estimation: PCA-based
				o normal orientation: based on:
					Hugues Hoppe, Tony DeRose, Tom Duchamp, John McDonald, and
					Werner Stuetzle. Surface reconstruction from unorganized
					points.	In Computer Graphics (Proc. SIGGRAPH '90),
					volume 26, pages 71–77, 1992.
		Input:
			- pcd: input point set containing optical centers
			- k: number of neighbors for normal computation
	*/
	std::cout << "---> Point set normal estimation and orientation" << std::endl;
	pcd.add_normal_map();

	// Estimate normals:
	std::cout << " normal estimation" << std::endl;
	CGAL::pca_estimate_normals<CGAL::Sequential_tag>(
									pcd,
									k,
									CGAL::parameters::
										point_map (pcd.point_map()).
										normal_map (pcd.normal_map())
									);

	// Try to consistently orient normals:
	// CGAL::mst_orient_normals() modifies the order of
	// input points so as to pack all sucessfully oriented
	// points first, and returns an iterator over the
	// first point with an unoriented normal
	std::cout << " normal orientation" << std::endl;
	Point_set::iterator unoriented_points_begin = CGAL::mst_orient_normals(
											pcd,
											k,
											CGAL::parameters::
												point_map (pcd.point_map()).
												normal_map (pcd.normal_map())
											);

	// Remove points which orientation failed:
	int nbRmPoints = pcd.end()-unoriented_points_begin;
	pcd.remove(unoriented_points_begin, pcd.end());
	std::cout << " Done with point set normal estimation and orientation: "
				<< nbRmPoints << " points removed" << std::endl << std::endl;
}

Point_set compute_and_orient_normals_based_on_origin(Point_set &pcd, unsigned int k){
	/*
		Usage:
			- compute and orient normals of a point set containing optical centers
				o normal estimation: PCA-based
				o normal orientation:
					normal n_i and vec(PiOi) must have positive dot product
		Input:
			- pcd: input point set containing optical centers
			- k: number of neighbors for normal computation
		Output:
			- out_pcd: point set containing the same points as in:pcd but with
					   oriented normals
	*/
	std::cout << "---> Point set normal estimation and orientation" << std::endl;
	pcd.add_normal_map();

	// Estimate normals:
	std::cout << " normal estimation" << std::endl;
	CGAL::pca_estimate_normals<CGAL::Sequential_tag>(
									pcd,
									k,
									CGAL::parameters::
										point_map (pcd.point_map()).
										normal_map (pcd.normal_map())
									);

	// Read optical centers and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcd.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcd.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcd.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cout << "ERROR: optical centers absent from in:pcd";
	}

	// New point set for output with normals but without optical centers:
	Point_set out_pcd; out_pcd.add_normal_map();

	std::cout << " orientation based on optical centers" << std::endl;
	for (Point_set::iterator p = pcd.begin(); p != pcd.end(); ++p)
	{
		Point Pi = pcd.point(*p);
		Vector ni = pcd.normal(*p);
		Point Oi(x_origin[*p], y_origin[*p], z_origin[*p]);
		Vector PiOi(Pi,Oi);

		if (CGAL::scalar_product(ni, PiOi) < 0) // if normal is mis-oriented
		{
			ni.operator*=(-1); // flip normal
		}
		out_pcd.insert(Pi, ni);
	}
	std::cout << " Done with normal estimation and orientation" << std::endl << std::endl;
	return out_pcd;
}

Point_set initialize_point_set(int &outProperty, X_Origin_Map &x_origin,
												 Y_Origin_Map &y_origin,
												 Z_Origin_Map &z_origin){
	/*
		Usage: Builds a Point_set object and add properties based on the value
		of 'outProperty':
			- 0: vertex position ONLY
			- 1: vertex + NORMAL
			- 2: vertex + OPTICAL CENTER
	*/
	Point_set point_set; // for rayshooting output storage
	std::cout << "Output Point Cloud initialized with: ";
	if (outProperty == 0){
	std::cout << "Vertex position ONLY" << std::endl;

	} else if (outProperty == 1) {
		std::cout << "Vertex position + NORMAL" << std::endl;
		point_set.add_normal_map(); // add normal property

	} else if (outProperty == 2) {
		std::cout << "Vertex position + OPTICAL CENTER" << std::endl;
		bool success_x = false; bool success_y = false; bool success_z = false;
		boost::tie (x_origin, success_x) = point_set.add_property_map<double>("x_origin", 0);
		boost::tie (y_origin, success_y) = point_set.add_property_map<double>("y_origin", 0);
		boost::tie (z_origin, success_z) = point_set.add_property_map<double>("z_origin", 0);
		if (!success_x || !success_y || !success_z){
			std::cout << "ERROR: impossible to add optical center properties";
		}
	} else {
		std::cerr << "ERROR: 'outProperty' argument not valid. Must be 0, 1 or 2" << std::endl;
		std::cerr << " 'outProperty' changed to 0 value" << std::endl;
		outProperty = 0;
	}
	return point_set;
}

void add_desired_output_to_pcd(Mesh &mesh, Tree &tree, Point &M, Ray &ray, Point_set &point_set, int &outProperty, X_Origin_Map &x_origin,
																												   Y_Origin_Map &y_origin,
																												   Z_Origin_Map &z_origin){
	/*
		Usage: Computes the intersection of the ray with the mesh and adds to the point cloud:
			   (based on the value of 'outProperty')
			   		- 0: vertex position ONLY
			 	 	- 1: vertex + NORMAL
			 	 	- 2: vertex + OPTICAL CENTER
	*/
	Ray_intersection intersection = tree.first_intersection(ray);
	if(intersection){
		if(boost::get<Point>(&(intersection->first))){
			const Point& p = boost::get<Point>(intersection->first);
			if (outProperty == 0){
				// only add point location:
				point_set.insert (p);

			} else if (outProperty == 1) {
				// compute and add corresponding face normal:
				const face_descriptor f = boost::get<face_descriptor>(intersection->second);
				const Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f,mesh);
				point_set.insert(p,n);

			} else if (outProperty == 2) {
				// add position of optical center:
				Point_set::iterator it = point_set.insert(p);
				x_origin[*it] = M.x();
				y_origin[*it] = M.y();
				z_origin[*it] = M.z();
			}
		}
	}
}

Point_set aerial_lidar(Mesh &mesh, Point &A, Point &B, double v0, double omega, double theta_0, double freq, int outProperty){
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
			int outProperty: desired property of output point set:
			 	 	- 0: vertex position ONLY
			 	 	- 1: vertex + NORMAL
			 	 	- 2: vertex + OPTICAL CENTER
		output:
			Point_set point_set: intersections of the rays with the mesh
	*/
	std::cout << "---> Aerial LiDAR acquisition from A(" << A << ") to B(" << B << ")" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);

	// Point Set initialization:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	Point_set point_set = initialize_point_set(outProperty, x_origin, y_origin, z_origin);
	Point_set positionsM;

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
			positionsM.insert(M);
			positionsM.insert(operator+(M, 50*ML));
			Ray ray(M, ML);
			add_desired_output_to_pcd(mesh, tree, M, ray, point_set, outProperty, x_origin,
																				  y_origin,
																				  z_origin);
		} else {
			// do nothing
		}
		t += dt; // update time
	}
	std::ofstream ofM("output_data/positionsM.ply");
	CGAL::write_ply_point_set(ofM,positionsM);
	std::cout << " Done with Aerial LiDAR acquisition" << std::endl << std::endl;
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
	std::cout << "---> READING: " << fileName << std::endl;
	Mesh mesh;
	// Variables to store polygon soup:
	std::ifstream is (fileName);
	std::vector<Point> points;
	std::vector<std::vector<std::size_t>> polygons;

	// Extract points and polygons from stream 'is':
	CGAL::read_OBJ(is, points, polygons);	
	std::cout << " number of points: " << points.size() << std::endl;
	std::cout << " number of polygons: " << polygons.size() << std::endl;

	// Processing polygon soup:
	if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons))
	{
		// polygons do not define a valid mesh
		std::cout << " Warning: [" << fileName << "] does not define a valid polygon mesh. "
		<< "Trying to orient..." << std::endl;
		if ( !CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons) ){ // try to orient
			// orientation failed
			std::cout << "  Warning: the orientation of the polygon soup ["
			<< fileName << "] failed. => Points were duplicated." << std::endl;
		} else {
			// orientation succeeded
			std::cout << "  orientation of the polygon soup [" << fileName << "] was a success." << std::endl;
		}
	} else {
		// polygons define a valid mesh
		std::cout << " Success: [" << fileName << "] defines a valid polygon mesh " << std::endl;
	}

	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

	std::cout << " Number of points: " << mesh.number_of_vertices() << std::endl;
	std::cout << " Number of faces: " << mesh.number_of_faces() << std::endl;
	std::cout << " Done with reading" << std::endl << std::endl;

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

Point_set spherical_ray_shooting(Mesh mesh, Point origin, int nTheta, int nPhi, int outProperty){
	/*
		usage:
			Shoot rays in all directions from a fixed position in space
			and compute intersections with a mesh
		input:
			Mesh mesh: mesh to be raycasted/sampled
			Point origin: The position from which the rays must be shooted
			int nTheta: nb of polar angles in [0,pi] interval to generate
			int nPhi: nb of azimuthal angles in [0,2pi] interval to generate
			int outProperty: desired property of output point set:
						 	 	- 0: vertex position ONLY
						 	 	- 1: vertex + NORMAL
						 	 	- 2: vertex + OPTICAL CENTER
		output:
			Point_set point_set: intersections of the rays with the mesh

	*/
	std::cout << "Spherical Ray Shooting from position (" << origin << ")" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);

	// Point Set initialization:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	Point_set point_set = initialize_point_set(outProperty, x_origin, y_origin, z_origin);

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
			add_desired_output_to_pcd(mesh, tree, origin, ray, point_set, outProperty, x_origin,
																					   y_origin,
																					   z_origin);
		}
	}
	std::cout << " Done with Spherical Ray Shooting" << std::endl << std::endl;
	return point_set;
}

void Strasbourg_Scene_Raytracing(int nTheta, int nPhi, int outProperty){
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
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi, outProperty);
	std::ofstream of("output_data/out_pcd_normal_PC3E45_3.ply");
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/out_pcd_normal_PC3E45_3.ply" << std::endl << std::endl;
}

void Strasbourg_Scene_Aerial_Lidar(double v0, double omega, double theta_0, double freq, int outProperty){
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
	// double z_above = bbox.zmin() + 100*( bbox.zmax() - bbox.zmin() );
	double z_above = 1000;

	// Point A(x_mid, y_A, z_above);
	// Point B(x_mid, y_B, z_above);

	Point A(14174.3, 20500.2, 147);
	Point B(14148.6, 20658.6, 147);

	Point_set point_set = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, outProperty);
	std::ofstream of("output_data/out_origin_street_n.ply");
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/out_origin_street_n.ply" << std::endl;
}

void object_raytracing_from_centroid(const char* filename, int nTheta, int nPhi){
	/*
		Routine to raytrace any closed object from its centroid and
		write the results to a ply pcd with normals
	*/
	Mesh mesh = read_input_file(filename);
	Point origin_lidar(centroid(mesh));
	int outProperty = 2;
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi, outProperty);
	std::ofstream of("output_data/tetra_oc.ply");
	// CGAL::write_off_point_set(of,point_set);
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/tetra_oc.ply" << std::endl;

	// just change the name of the file and user this function for a PLY pcd:
	// CGAL::write_ply_point_set(ofPLY,point_set);
}

void add_normal_noise(Point_set &pcd, double mu, double sigma){
	/*
		usage:
			Add depth-noise (along the ray) folowwing a normal distribution
			of parameters (mu, sigma) to a given point set.
		input:
			Point_set pcd: noise-free point set to which adding the noise
								 Must contain XYZ Origin Maps with the
								 corresponding optical centers.
			double mu, double sigma: parameters of the normal distribution
	*/
	std::cout << "---> Adding normal noise (" << mu <<", " << sigma << ")"
				<< " to point set" << std::endl;

	// Normal distribution (mu, sigma):
	std::default_random_engine generator;
	std::normal_distribution<double> normDistri(mu, sigma);

	// Read optical centers and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcd.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcd.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcd.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cout << "ERROR: optical centers absent from in:pcd";
	}

	// Browse point set:
	for (Point_set::const_iterator it = pcd.begin(); it != pcd.end(); it++){
		Point M = Point(x_origin[*it], y_origin[*it], z_origin[*it]); // origin
		Point Pi = pcd.point(*it); // sample point corresponding to origin
		Vector MPi = Vector(M,Pi);
		Vector d = normalize_vector_3(MPi); // normalized ray
		double k = normDistri(generator); // scalar following normal distribution

		// translates sample point along the ray by a factor of 'k':
		pcd.point(*it).operator+=( d.operator*=(k) );
	}
}

int main(int argc, char* argv[])
{
	std::cout << "--> program started <--" << std::endl << std::endl;
	const char* filename = (argc > 1) ? argv[1] : "input_data/light/cube.obj"; 
	// object_raytracing_from_centroid(filename, 150, 300);

  
	double v0 = 50;
	double omega = 200*M_PI;
	double freq = 300000;
	double theta_0 = 0;
	int outProperty = 1;
	/*
		- 0: vertex position ONLY
		- 1: vertex + NORMAL
		- 2: vertex + OPTICAL CENTER
	*/
	// Strasbourg_Scene_Aerial_Lidar(v0, omega, theta_0, freq, outProperty);
	// object_raytracing_from_centroid("input_data/light/tetrahedron.off", 5, 10);

	Point_set pcd;
	std::ifstream is_pcd ("output_data/tetra_oc.ply");
	// std::ifstream is_pcd ("output_data/out_cow.off");
	is_pcd >> pcd;
	// compute_and_orient_normals(pcd, 18);
	Point_set out_pcd = compute_and_orient_normals_based_on_origin(pcd, 5);
	std::ofstream of("output_data/tetra_normal.ply");
	CGAL::write_ply_point_set(of,out_pcd);


  std::cout << "--> program ended <--" << std::endl;
  return 0;
}