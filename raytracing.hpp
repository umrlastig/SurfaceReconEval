#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>

// Accelerating data structure
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

// Geometrical objects
#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>

// Polygon mesh processing
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

// Normal computation
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

// Kernel
#include <CGAL/Simple_cartesian.h>

// Input / Output
#include "IO.hpp"

#include <CGAL/tags.h>


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

namespace PMP = CGAL::Polygon_mesh_processing;

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
		double areaT = PMP::face_area(fd,mesh);

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

void compute_local_frame(Vector &AB, Vector &vec_i, Vector &vec_j, Vector &vec_k, bool verbose){
	if (verbose) std::cout << " > Computing local frame" << std::endl;
	vec_k = normalize_vector_3(AB);

	vec_j = CGAL::cross_product(Vector(0,0,1), vec_k);
	vec_j = normalize_vector_3(vec_j);

	vec_i = CGAL::cross_product(vec_j, vec_k);
	if (verbose) {
		std::cout << "  vec_k: [" << vec_k << "]" << std::endl;
		std::cout << "  vec_i: [" << vec_i << "]" << std::endl;
		std::cout << "  vec_j: [" << vec_j << "]" << std::endl;
	}
}

void compute_and_orient_normals(Point_set &pcd, unsigned int k, bool verbose){
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
	if (verbose) std::cout << "\n---> Point set normal estimation and orientation" << std::endl;
	pcd.add_normal_map();

	// Estimate normals:
	if (verbose) std::cout << "     > normal estimation" << std::endl;
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
	if (verbose) std::cout << "     > normal orientation" << std::endl;
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
	if (verbose) std::cout << "     Done with point set normal estimation and orientation: "
				<< nbRmPoints << " points removed" << std::endl;
}

Point_set compute_and_orient_normals_based_on_origin(Point_set &pcd, unsigned int k, bool verbose){
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
			- outPcd: point set containing the same points as in:pcd but with
					   oriented normals
	*/
	if (verbose) std::cout << "\n---> Point set normal estimation and orientation" << std::endl;
	pcd.add_normal_map();

	// Estimate normals:
	if (verbose) std::cout << "     > normal estimation" << std::endl;
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
		std::cerr << "ERROR: optical centers absent from in:pcd";
	}

	// New point set for output with normals but without optical centers:
	Point_set outPcd; outPcd.add_normal_map();

	// Orient normals:
	if (verbose) std::cout << "     > orientation based on optical centers" << std::endl;
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
		outPcd.insert(Pi, ni);
	}

	pcd.remove_normal_map();
	if (verbose) std::cout << "Done with normal estimation and orientation" << std::endl;
	return outPcd;
}

void remove_optical_centers(Point_set &pcd){

	// Read optical center property maps and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcd.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcd.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcd.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: optical centers absent from in:pcd";
	}
	
	// Remove property maps:
	success_x = pcd.remove_property_map(x_origin);
	success_y = pcd.remove_property_map(y_origin);
	success_z = pcd.remove_property_map(z_origin);
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: impossible to remove optical centers from in:pcd";
	}
}

Point_set initialize_point_set(int &outProperty, X_Origin_Map &x_origin,
												 Y_Origin_Map &y_origin,
												 Z_Origin_Map &z_origin,
												 bool verbose){
	/*
		Usage: Builds a Point_set object and add properties based on the value
		of 'outProperty':
			- 0: vertex position ONLY
			- 1: vertex + NORMAL
			- 2: vertex + OPTICAL CENTER
	*/
	Point_set point_set; // for rayshooting output storage
	if (verbose) std::cout << "     Output Point Cloud initialized with: ";
	if (outProperty == 0){
	if (verbose) std::cout << "Vertex position ONLY" << std::endl;

	} else if (outProperty == 1) {
		if (verbose) std::cout << "Vertex position + NORMAL" << std::endl;
		point_set.add_normal_map(); // add normal property

	} else if (outProperty == 2) {
		if (verbose) std::cout << "Vertex position + OPTICAL CENTER" << std::endl;
		bool success_x = false; bool success_y = false; bool success_z = false;
		boost::tie (x_origin, success_x) = point_set.add_property_map<double>("x_origin", 0);
		boost::tie (y_origin, success_y) = point_set.add_property_map<double>("y_origin", 0);
		boost::tie (z_origin, success_z) = point_set.add_property_map<double>("z_origin", 0);
		if (!success_x || !success_y || !success_z){
			std::cerr << "ERROR: impossible to add optical center properties";
		}
	} else {
		std::cerr << "ERROR: 'outProperty' argument not valid. Must be 0, 1 or 2" << std::endl;
		std::cerr << " 'outProperty' changed to 0 value" << std::endl;
		outProperty = 0;
	}
	return point_set;
}

void add_desired_output_to_pcd(Mesh &mesh, Tree &tree, Point &M, Ray &ray, Point_set &point_set,
																		   std::vector<face_descriptor> &hitFaces,
																		   int &outProperty,
																		   X_Origin_Map &x_origin,
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
				point_set.insert(p);

			} else if (outProperty == 1) {
				// compute and add corresponding face normal:
				const face_descriptor f = boost::get<face_descriptor>(intersection->second);
				const Vector n = PMP::compute_face_normal(f,mesh);
				point_set.insert(p,n);

			} else if (outProperty == 2) {
				// add position of optical center:
				Point_set::iterator it = point_set.insert(p);
				x_origin[*it] = M.x();
				y_origin[*it] = M.y();
				z_origin[*it] = M.z();
			}
			const face_descriptor fd = boost::get<face_descriptor>(intersection->second);
			if (std::find(hitFaces.begin(), hitFaces.end(), fd) == hitFaces.end()){
				// fd is not already in hitFaces so add it
				hitFaces.push_back(fd);
			}
		}
	}
}

Point_set parallel_line_aerial_lidar(Mesh &mesh, Point &A, Point &B, double v0,
													   double omega,
													   double fov,
													   double theta_0,
													   double freq,
													   int outProperty,
													   bool verbose){
	/*
		usage:
			Simulate an aerial LiDAR acquisition
			following a parallel-line pattern
		input:
			Mesh mesh: mesh to be raycasted/sampled
			Point A: starting point of the vehicle
			Point B: ending point of the vehicle
			double v0: speed of the vehicle [m/s]
			double omega: angular speed of the laser pointer [rad/s]
			double fov: field of view of LiDAR scanner [rad]
			double theta_0: initial angle of the laser pointer [rad]
			double freq: frequency (time resolution) at which [Hz]
						 the analysis should be carried out
			int outProperty: desired property of output point set:
			 	 	- 0: vertex position ONLY
			 	 	- 1: vertex + NORMAL
			 	 	- 2: vertex + OPTICAL CENTER
		output:
			Point_set point_set: intersections of the rays with the mesh
	*/
	if (verbose) std::cout << "\n---> Aerial LiDAR acquisition from A[" << A << "] to B[" << B << "]"
							<< "\nScanning pattern: parallel-line" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	std::vector<face_descriptor> hitFaces;

	// Point Set initialization:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	Point_set point_set = initialize_point_set(outProperty, x_origin, y_origin, z_origin, false);
	Point_set positionsM;

	Vector AB = Vector(A,B);
	Vector vec_i; Vector vec_j; Vector vec_k;
	compute_local_frame(AB, vec_i, vec_j, vec_k, false);

	Point M; Vector di; Vector dj; Vector ML; double theta; // Some necesary variables

	double beta = fov / 2; // max semi-angle for a ray being actually shooted
	double tB = sqrt(AB.squared_length()) / v0; // duration of simulation
	int nT = (int) (freq * tB); // number of time steps
	double dt = tB / (nT-1); // infinitesimal unit of time
	double t = 0; // initial time

	for (int ti=0; ti<nT; ti++){
		theta = omega*t + theta_0; // angular position
		theta = std::fmod(theta, 2*M_PI); // euclidean division
		if (theta <= beta || theta >= (2*M_PI - beta)){ // theta is eligible
			M = operator+(A, operator*(v0*t, vec_k)); // position of M
			di = operator*(cos(theta), vec_i);
			dj = operator*(sin(theta), vec_j);
			ML = di.operator+(dj); // laser pointer
			positionsM.insert(M);
			positionsM.insert(operator+(M, 50*ML));
			Ray ray(M, ML);
			add_desired_output_to_pcd(mesh, tree, M, ray, point_set, hitFaces, outProperty,
																			   x_origin,
																			   y_origin,
																			   z_origin);
		} else {
			// do nothing
		}
		t += dt; // update time
	}
	write_point_set("output_data/positionsM.ply", positionsM, false);
	double aeraHit = PMP::area(hitFaces, mesh);
	double aeraTot = PMP::area(mesh);
	if (verbose) std::cout << "Total area hit during acquisition: "
		<< aeraHit << " (" << aeraHit / aeraTot * 100 <<' '<<'%'<<" of total)" << std::endl;
	if (verbose) std::cout << "Done with Aerial LiDAR acquisition" << std::endl;
	return point_set;
}

Point_set elliptical_aerial_lidar(Mesh &mesh, Point &A, Point &B, double v0,
														double omega,
														double theta,
														double phi_0,
														double freq,
														int outProperty,
														bool verbose){
	/*
		usage:
			Simulate an aerial LiDAR acquisition
			following an elliptical pattern
		input:
			Mesh mesh: mesh to be raycasted/sampled
			Point A: starting point of the vehicle
			Point B: ending point of the vehicle
			double v0: speed of the vehicle [m/s]
			double omega: angular speed of the laser pointer [rad/s]
			double theta: azimuthal angle of the laser pointer [rad]
			double phi_0: initial polar angle of the laser pointer [rad]
			double freq: frequency (time resolution) at which [Hz]
						 the analysis should be carried out
			int outProperty: desired property of output point set:
					- 0: vertex position ONLY
					- 1: vertex + NORMAL
					- 2: vertex + OPTICAL CENTER
		output:
			Point_set point_set: intersections of the rays with the mesh
	*/
	if (verbose) std::cout << "\n---> Aerial LiDAR acquisition from A[" << A << "] to B[" << B << "]"
							<< "\nScanning pattern: elliptical" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	std::vector<face_descriptor> hitFaces;

	// Point Set initialization:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	Point_set point_set = initialize_point_set(outProperty, x_origin, y_origin, z_origin, false);
	Point_set positionsM;

	Vector AB = Vector(A,B);
	Vector vec_i; Vector vec_j; Vector vec_k;
	compute_local_frame(AB, vec_i, vec_j, vec_k, false);

	Point M; Vector di; Vector dj; Vector dk; Vector ML; double phi; // Some necesary variables

	// double beta = fov / 2; // max semi-angle for a ray being actually shooted
	double tB = sqrt(AB.squared_length()) / v0; // duration of simulation
	int nT = (int) (freq * tB); // number of time steps
	double dt = tB / (nT-1); // infinitesimal unit of time
	double t = 0; // initial time

	for (int ti=0; ti<nT; ti++){
		phi = omega*t + phi_0; // angular position
		phi = std::fmod(phi, 2*M_PI); // euclidean division
		M = operator+(A, operator*(v0*t, vec_k)); // position of M
		di = operator*(-cos(theta), vec_i);
		dj = operator*(-sin(theta)*cos(phi), vec_j);
		dk = operator*(sin(theta)*sin(phi), vec_k);
		ML = di.operator+(dj).operator+(dk); // laser pointer
		positionsM.insert(M);
		positionsM.insert(operator+(M, 50*ML));
		Ray ray(M, ML);
		add_desired_output_to_pcd(mesh, tree, M, ray, point_set, hitFaces, outProperty, x_origin,
																			  y_origin,
																			  z_origin);
		t += dt; // update time
	}
	write_point_set("output_data/positionsM.ply", positionsM, false);
	double aeraHit = PMP::area(hitFaces, mesh);
	double aeraTot = PMP::area(mesh);
	if (verbose) std::cout << "Total area hit during acquisition: "
		<< aeraHit << " (" << aeraHit / aeraTot * 100 <<' '<<'%'<<" of total)" << std::endl;
	if (verbose) std::cout << "Done with Aerial LiDAR acquisition" << std::endl;
	return point_set;
}

Point_set spherical_ray_shooting(Mesh mesh, Point origin, int nTheta, int nPhi, int outProperty, bool verbose){
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
	std::cout << "\nSpherical Ray Shooting from position [" << origin << "]" << std::endl;

	Tree tree(faces(mesh).first, faces(mesh).second, mesh);
	std::vector<face_descriptor> hitFaces;

	// Point Set initialization:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	Point_set point_set = initialize_point_set(outProperty, x_origin, y_origin, z_origin, verbose);

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
			add_desired_output_to_pcd(mesh, tree, origin, ray, point_set, hitFaces,
																		  outProperty,
																		  x_origin,
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
	const char* inMeshFile = "input_data/heavy/open_data_strasbourg/PC3E45/PC3E45_3.obj";
	bool verbose = true;
	Mesh mesh = read_mesh<Mesh, Point>(inMeshFile, verbose);

	CGAL::Bbox_3 bbox = PMP::bbox(mesh);
	double x_origin = ( bbox.xmin() + bbox.xmax() ) / 2;
	double y_origin = ( bbox.ymin() + bbox.ymax() ) / 2;
	double z_origin = bbox.zmin() + 3*( bbox.zmax() - bbox.zmin() );
	Point origin_lidar(x_origin, y_origin, z_origin);
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi, outProperty, verbose);

	const char* outFileName = "output_data/out_pcd_normal_PC3E45_3.ply";
	write_point_set(outFileName, point_set, verbose);
}

void object_raytracing_from_centroid(const char* fileName, int nTheta, int nPhi){
	/*
		Routine to raytrace any closed object from its centroid and
		write the results to a ply pcd with normals
	*/
	bool verbose = true;
	Mesh mesh = read_mesh<Mesh,Point>(fileName, verbose);

	Point origin_lidar(centroid(mesh));
	int outProperty = 2;
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi, outProperty, verbose);

	const char* outFile("output_data/tetra_oc.ply");
	write_point_set(outFile, point_set, verbose);
}

void add_depth_normal_noise(Point_set &pcd, double mu, double sigma, bool verbose){
	/*
		usage:
			Add depth-noise (along the ray) following a normal distribution
			of parameters (mu, sigma) to a given point set.
		input:
			Point_set pcd: noise-free point set to which adding the noise
								 Must contain XYZ Origin Maps with the
								 corresponding optical centers.
			double mu, double sigma: parameters of the normal distribution
	*/
	if (verbose) std::cout << "---> Adding normal noise (" << mu <<", " << sigma << ")"
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
		std::cerr << "ERROR: optical centers absent from in:pcd";
	}

	// Browse point set:
	for (Point_set::const_iterator it = pcd.begin(); it != pcd.end(); it++){
		Point M = Point(x_origin[*it], y_origin[*it], z_origin[*it]); // origin
		Point Pi = pcd.point(*it); // sample point corresponding to origin
		Vector MPi = Vector(M,Pi);
		Vector d = normalize_vector_3(MPi); // normalized ray
		double k = normDistri(generator); // scalar following normal distribution

		// Translates sample point along the ray by a factor of 'k':
		pcd.point(*it).operator+=( d.operator*=(k) );
	}
}

void add_normal_noise(Point_set &pcd, double muXY, double sigmaXY,
									  double muZ, double sigmaZ,
									  bool verbose){
	/*
		usage:
			Add altimetric and planimetric noise folowwing normal
			distributions to a given point set.
		input:
			Point_set pcd: noise-free point set to which adding the noise
			double muZ, double sigmaZ: parameters of altimetric error
			double muXY, double sigmaXY: parameters of planimetric error
	*/
	if (verbose) std::cout << "\n---> Adding planimetric noise (" << muXY <<", " << sigmaXY << ")"
				<< " and altimetric noise (" << muZ <<", " << sigmaZ << ")"
				<< " to point set" << std::endl;

	// Normal distributions (mu, sigma):
	std::default_random_engine generator;
	std::normal_distribution<double> planimetricError(muXY, sigmaXY);
	std::normal_distribution<double> altimetricError(muZ, sigmaZ);

	// Browse point set:
	for (Point_set::const_iterator it = pcd.begin(); it != pcd.end(); it++){
		Point Pi = pcd.point(*it); // sample point

		// Scalar errors:
		double dX = planimetricError(generator);
		double dY = planimetricError(generator);
		double dZ = altimetricError(generator);
		Vector d = Vector(dX,dY,dZ);
		// std::cout <<"["<<dX<<", "<<dY<<", "<<dZ<<"]," << std::endl;

		// Translates sample point along the ray by a factor of 'k':
		pcd.point(*it).operator+=(d);
	}
}

Point_set optical_centers_2_rays(Point_set pcdOC){
	/*
		usage:
			Look for sensor positions in in:pcdOC, computes optical rays
			and appends them into the normal field of a new point cloud
			so that optical rays can be visualized as normals in a mesh
			viewer like Meshlab
		input:
			Point_set pcdOC: point cloud which must contain optical centers
		output:
			Point_set pcdRays: point cloud with optical rays as normal field
	*/
	Point_set pcdRays; pcdRays.add_normal_map();

	// Read optical centers and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcdOC.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcdOC.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcdOC.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: optical centers absent from in:pcdOC";
	}

	// Browse pcdOC:
	for (Point_set::iterator p = pcdOC.begin(); p != pcdOC.end(); ++p)
	{
		Point Pi = pcdOC.point(*p);
		Point Oi(x_origin[*p], y_origin[*p], z_origin[*p]);
		Vector PiOi(Pi,Oi);
		pcdRays.insert(Pi, PiOi);
	}
	return pcdRays;
}

Point_set getOrigins(Point_set pcdOC){
	/*
		usage:
			Look for sensor positions in in:pcdOC and appends them into
			a new point cloud so that sensor positions can be visualized
			in a mesh viewer like Meshlab
		input:
			Point_set pcdOC: point cloud which must contain optical centers
		output:
			Point_set pcdRays: point cloud with optical centers only
	*/
	Point_set pcdRays;

	// Read optical centers and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcdOC.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcdOC.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcdOC.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: optical centers absent from in:pcdOC";
	}

	// Browse pcdOC:
	for (Point_set::iterator p = pcdOC.begin(); p != pcdOC.end(); ++p)
	{
		// Point Pi = pcdOC.point(*p);
		Point Oi(x_origin[*p], y_origin[*p], z_origin[*p]);

		pcdRays.insert(Oi);
	}
	return pcdRays;
}
