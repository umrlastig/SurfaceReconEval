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
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Mesh::Face_range Face_range;
typedef Tree::Primitive_id Primitive_id;

// Optical Centers
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
	std::cout << "Computing local frame" << std::endl;
	vec_k = normalize_vector_3(AB);

	vec_j = CGAL::cross_product(Vector(0,0,1), vec_k);
	vec_j = normalize_vector_3(vec_j);

	vec_i = CGAL::cross_product(vec_j, vec_k);

	std::cout << "vec_k: [" << vec_k << "]" << std::endl;
	std::cout << "vec_i: [" << vec_i << "]" << std::endl;
	std::cout << "vec_j: [" << vec_j << "]" << std::endl;
}

void compute_and_orient_normals(Point_set &point_set, int &k){
	std::cout << "---> Point set normal estimation" << std::endl;
	point_set.add_normal_map();

	// Estimate normals:
	CGAL::pca_estimate_normals<CGAL::Sequential_tag>(
								point_set.begin(),
								point_set.end(),
								point_set.point_map(),
								point_set.normal_map(),
								k
								);

	// Try to consistently orient normals:
	std::cout << " normal orientation" << std::endl;
	Point_set::iterator unoriented_points_begin = CGAL::mst_orient_normals(
														point_set.begin(),
														point_set.end(),
														point_set.point_map(),
														point_set.normal_map(),
														k
														);

	// Remove points which orientation failed:
	int nbRmPoints = point_set.end()-unoriented_points_begin;
	point_set.remove(unoriented_points_begin, point_set.end());
	std::cout << " Done with point set normal estimation and orientation: "
				<< nbRmPoints << " points removed" << std::endl << std::endl;
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
		bool success = false;
		boost::tie (x_origin, success) = point_set.add_property_map<double>("x_origin", 0);
		boost::tie (y_origin, success) = point_set.add_property_map<double>("y_origin", 0);
		boost::tie (z_origin, success) = point_set.add_property_map<double>("z_origin", 0);
		assert(success);
	} else {
		std::cerr << "Error: 'outProperty' argument not valid. Must be 0, 1 or 2" << std::endl;
		std::cerr << "'outProperty' changed to 0 value" << std::endl;
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
	std::cout << "Done with Spherical Ray Shooting" << std::endl << std::endl;
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

	Point A(x_mid, y_A, z_above);
	Point B(x_mid, y_B, z_above);

	Point_set point_set = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, outProperty);
	std::ofstream of("output_data/out_origin.ply");
	CGAL::write_ply_point_set(of,point_set);
	std::cerr << "wrote output_data/out_origin.ply" << std::endl;
}

void object_raytracing_from_centroid(const char* filename, int nTheta, int nPhi){
	/*
		Routine to raytrace any closed object from its centroid and
		write the results to a ply pcd with normals
	*/
	Mesh mesh = read_input_file(filename);
	Point origin_lidar(centroid(mesh));
	Point_set point_set = spherical_ray_shooting(mesh, origin_lidar, nTheta, nPhi, 0);
	std::ofstream of("output_data/out_pcd_normal_00AA.off");
	CGAL::write_off_point_set(of,point_set);
	std::cerr << "wrote output_data/out_pcd_normal_00AA.off" << std::endl;

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

	// extract property maps:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success = false;
	boost::tie (x_origin, success) = pcd.property_map<FT>("x_origin"); assert(success);
	boost::tie (y_origin, success) = pcd.property_map<FT>("y_origin"); assert(success);
	boost::tie (z_origin, success) = pcd.property_map<FT>("z_origin"); assert(success);

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

	CGAL::Face_around_target_iterator<Mesh> face_b, face_e;

	for (vertex_descriptor vd : vertices(mesh)){
		Point p = mesh.point(vd);
		std::cout << vd << " : " << p << std::endl;

		//check if p is in mAlpha_v:
		if (std::find(mAlpha_v.begin(), mAlpha_v.end(), vd) == mAlpha_v.end()){ // if not already in mAlpha_v
			mAlpha_v.push_back(vd); // ... add it
		}
		// browse all faces incident to vertex vd:
		for (boost::tie( face_b, face_e) = CGAL::faces_around_target(mesh.halfedge(vd),mesh);
		face_b != face_e;
		++face_b)
		{
			std::cout << " : " << *face_b << std::endl;

			// check if face_b is in mAlpha_f:
			if (std::find(mAlpha_f.begin(), mAlpha_f.end(), *face_b) == mAlpha_f.end()){ // if not already in mAlpha_f
				mAlpha_f.push_back(*face_b); // ... add it

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




	// std::ofstream of("output_data/noisy_pcd.ply");

	// CGAL::write_ply_point_set(of,pcd);

}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "input_data/light/cube.obj"; 
  // object_raytracing_from_centroid(filename, 150, 300);

  
  // double v0 = 50;
  // double omega = 200*M_PI;
  // double freq = 300000;
  // double theta_0 = 0;
  // int outProperty = 2;
  	/*
 	- 0: vertex position ONLY
	- 1: vertex + NORMAL
	- 2: vertex + OPTICAL CENTER
	*/
  // Strasbourg_Scene_Aerial_Lidar(v0, omega, theta_0, freq, outProperty);

  // Strasbourg_Scene_Raytracing(150, 300, outProperty);

  test();

  std::cerr << "done" << std::endl;
  return 0;
}