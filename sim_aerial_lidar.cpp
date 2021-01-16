#include "raytracing.hpp"


int main(int argc, char* argv[])
{
	if (argc != 6)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	const char* inFileName = argv[1];
	const std::string outBaseName = argv[2];
	const std::string ext = argv[3];
	double sigma = std::atof(argv[4]);
	const std::string verb = argv[5];

	// Set verbose flag
	bool verbose = (verb == "1") ? true : false;
	if (verbose) std::cout << "verbose flag active" << std::endl;
  
	double v0 = 50;
	double omega = 200*M_PI;
	double freq = 300000;
	double theta_0 = 0;
	/* outProperty:
		- 0: vertex position ONLY
		- 1: vertex + NORMAL
		- 2: vertex + OPTICAL CENTER
	*/

	Mesh mesh = read_mesh<Mesh,Point>(inFileName, verbose);

	// Starting and ending positions
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	double xMid = ( bbox.xmin() + bbox.xmax() ) / 2;
	double yA = bbox.ymin();
	double yB =  bbox.ymax();
	double altitude = 500; // 150 --> 900 for a low-altitude LiDAR

	Point A(xMid, yA, altitude);
	Point B(xMid, yB, altitude);

	// Street positions:
	// Point A(14174.3, 20500.2, 147);
	// Point B(14148.6, 20658.6, 147);

	// Output file names
	std::string outFileOC = outBaseName + std::string("_OptCtr") + ext;
	std::string outFileN = outBaseName + std::string("_normals") + ext;
	std::string outFileP = outBaseName + std::string("_pts") + ext;

	if (sigma == 0){
		if (verbose) std::cout << "\nPerfect aerial LiDAR scan" << std::endl;

		// points only
		Point_set pcdPts = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, 0);
		write_point_set(outFileP.c_str(), pcdPts, verbose);

		// normals
		Point_set pcdN = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, 1);
		write_point_set(outFileN.c_str(), pcdN, verbose);

		// optical centers
		Point_set pcdOC = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, 2);
		write_point_set(outFileOC.c_str(), pcdOC, verbose);
	} else {
		if (verbose) std::cout << "\nRealistic aerial LiDAR scan" << std::endl;

		int outProperty = 2;
		Point_set pcdOC = aerial_lidar(mesh, A, B, v0, omega, theta_0, freq, outProperty);

		// normal noise
		double mu = 0;
		if (sigma != 0) add_normal_noise(pcdOC, mu, sigma);

		// optical centers		
		write_point_set(outFileOC.c_str(), pcdOC, verbose);

		// normals
		int kNN = 20;
		Point_set pcdN = compute_and_orient_normals_based_on_origin(pcdOC, kNN);		
		write_point_set(outFileN.c_str(), pcdN, verbose);

		// points only
		remove_optical_centers(pcdOC);		
		write_point_set(outFileP.c_str(), pcdOC, verbose);
	}

	return 0;
}
