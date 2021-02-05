#include "raytracing.hpp"

void display_help(char* argv[]){
	std::cout << "\n========= LIDAR SCAN SIMULATOR =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input-file, -i              |    File containing the virtual environment you wish to scan\n"
		<< " --out-base, -o                |    Base name of the output files where to store point clouds (with points only, normals, and optical centers"
		<< std::endl << std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --extension, -e               |    The extension output files should be exported with ('.ply' or '.off')\n"
		<< " --perfect-scan, -p            |    No noise is added to point locations and normals are perfect\n"
		<< " --verbose, -v                 |    Display information throughout execution\n"
		<< " FLIGHT parameters:\n"
		<< "  --flying-altitude, -z        |    [m]\n"
		<< "  --flying-speed, -v0          |    [m/s]\n"
		<< " LiDAR parameters:\n"
		<< "  --angular-speed, -omega      |    [rot/s]\n"
		<< "  --field-of-view, -fov        |    [deg]\n"
		<< "  --pulse-frequency, -fp       |    [Hz]\n"
		<< " NOISE parameters:\n"
		<< "  --planimetric-mean, -muXY    |    [m]\n"
		<< "  --planimetric-std, -sigmaXY  |    [m]\n"
		<< "  --altimetric-mean, -muZ      |    [m]\n"
		<< "  --altimetric-std, -sigmaZ    |    [m]\n"
		<< std::endl;
}

int main(int argc, char* argv[])
{
	// Default values for cmd-line parameters:
	std::string inFileName = "";
	std::string outBaseName = "";
	std::string ext = ".ply";
	double muXY = 0, sigmaXY = 0.13; // planimetric error
	double muZ = 0, sigmaZ = 0.05; // altimetric error
	bool perfectScan = false;
	bool verbose = false;

	// Flight parameters:
	double altitude = 1000; // [m] 150 --> 1000 for a low-altitude LiDAR
	double v0 = 60; // [m/s]
	double omega = 150; // [rotations/s]
	double fov = 40; // [deg] full angle range
	double freq = 400000; // [Hz]
	double theta_0 = 0;

	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help(argv);
			return 0;
		} else if (std::string(argv[i]) == "--input-file" || std::string(argv[i]) == "-i"){
			inFileName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--out-base" || std::string(argv[i]) == "-o"){
			outBaseName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--extension" || std::string(argv[i]) == "-e"){
			ext = argv[++i];
		} else if (std::string(argv[i]) == "--planimetric-mean" || std::string(argv[i]) == "-muXY"){
			muXY = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--planimetric-std" || std::string(argv[i]) == "-sigmaXY"){
			sigmaXY = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--altimetric-mean" || std::string(argv[i]) == "-muZ"){
			muZ = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--altimetric-std" || std::string(argv[i]) == "-sigmaZ"){
			sigmaZ = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--flying-altitude" || std::string(argv[i]) == "-z"){
			altitude = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--flying-speed" || std::string(argv[i]) == "-v0"){
			v0 = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--angular-speed" || std::string(argv[i]) == "-omega"){
			omega = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--field-of-view" || std::string(argv[i]) == "-fov"){
			fov = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--pulse-frequency" || std::string(argv[i]) == "-fp"){
			freq = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--perfect-scan" || std::string(argv[i]) == "-p"){
			perfectScan = true;
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help(argv);
			return 1;
		}
	}

	/* outProperty:
		- 0: vertex position ONLY
		- 1: vertex + NORMAL
		- 2: vertex + OPTICAL CENTER
	*/

	// Check if in/out files were set
	if (inFileName == "" || outBaseName == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help(argv);
		return 1;
	}

	if (verbose){
		std::cout << "verbose flag active\n" << std::endl;

		std::cout << "Simulating LiDAR scan with parameters:" << std::endl;
		std::cout << " - Input mesh: '" << inFileName << "'" << std::endl;
		std::cout << " - Output files: '" << outBaseName << "[...]" << ext << "'" << std::endl;
		if (perfectScan){
			std::cout << " - Perfect scan" << std::endl;
		} else {
			std::cout << " - Noise parameters:\n"
					<<"    > Planimetric error ~ N("<<muXY<<", "<<sigmaXY<<"^2)\n"
					<<"    > Altimetric error ~ N("<<muZ<<", "<<sigmaZ<<"^2)" << std::endl;
		}
		std::cout << " - Flying altitude: " << altitude << " m" << std::endl;
		std::cout << " - Flying speed: " << v0 << " m/s" << std::endl;
		std::cout << " - Angular speed: " << omega << " rotations/s ("
										  << omega/M_PI << " Hz)" << std::endl;
		std::cout << " - Field of view: " << fov << " deg" << std::endl;
		std::cout << " - Pulse frequency: " << freq << " Hz" << std::endl;
	}

	Mesh mesh = read_mesh<Mesh,Point>(inFileName.c_str(), verbose);

	// Starting and ending positions
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	double xMid = ( bbox.xmin() + bbox.xmax() ) / 2;
	double yA = bbox.ymin();
	double yB =  bbox.ymax();

	Point A(xMid, yA, altitude);
	Point B(xMid, yB, altitude);

	// Output file names
	std::string outFileOC = outBaseName + std::string("_OptCtr") + ext;
	std::string outFileN = outBaseName + std::string("_normals") + ext;
	std::string outFileP = outBaseName + std::string("_pts") + ext;

	// Unit conversions:
	omega *= 2*M_PI; // from [rot/s] to [rad/s]
	fov *= M_PI / 180; // from [deg] to [rad]

	if (perfectScan){
		if (verbose) std::cout << "\nPerfect aerial LiDAR scan" << std::endl;

		// points only
		Point_set pcdPts = aerial_lidar(mesh, A, B, v0, omega, fov, theta_0, freq, 0, verbose);
		write_point_set(outFileP.c_str(), pcdPts, verbose);

		// normals
		Point_set pcdN = aerial_lidar(mesh, A, B, v0, omega, fov, theta_0, freq, 1, verbose);
		write_point_set(outFileN.c_str(), pcdN, verbose);

		// optical centers
		Point_set pcdOC = aerial_lidar(mesh, A, B, v0, omega, fov, theta_0, freq, 2, verbose);
		write_point_set(outFileOC.c_str(), pcdOC, verbose);
	} else {
		if (verbose) std::cout << "\nRealistic aerial LiDAR scan" << std::endl;

		int outProperty = 2;
		Point_set pcdOC = aerial_lidar(mesh, A, B, v0, omega, fov, theta_0, freq, outProperty, verbose);

		// normal noise
		add_normal_noise(pcdOC, muXY, sigmaXY, muZ, sigmaZ);

		// optical centers		
		write_point_set(outFileOC.c_str(), pcdOC, verbose);

		// normals
		int k = 20;
		Point_set pcdN = compute_and_orient_normals_based_on_origin(pcdOC, k, verbose);
		write_point_set(outFileN.c_str(), pcdN, verbose);

		// points only
		remove_optical_centers(pcdOC);		
		write_point_set(outFileP.c_str(), pcdOC, verbose);
	}

	return 0;
}
