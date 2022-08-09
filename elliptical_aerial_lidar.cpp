#include "raytracing.hpp"

void display_help(char* argv[]){
	std::cout << "\n========= ELLIPTICAL LIDAR SCAN SIMULATOR =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input-file , -i              |    File containing the virtual environment you wish to scan\n"
		<< " --out-base , -o                |    Base name of the output files where to store point clouds (with points only, normals, and optical centers"
		<< std::endl << std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --extension , -e               |    The extension output files should be exported with ('.ply' or '.off')\n"
		<< " --noise-free , -nf             |    No noise is added to point locations\n"
		<< " --perfect-normals , -pn        |    Normals are the ones of the facet each sampled point belongs to\n"
		<< " --use-opt-ctrs , -useOC        |    Normals are estimated and their orientation is choosen according to the sensor position\n"
		<< " --verbose , -v                 |    Display information throughout execution\n"
		<< " FLIGHT parameters:\n"
		<< "  --alpha-X , -alphaX           |    in [0,1]: fraction of X amplitude along which to fly\n"
		<< "  --flying-altitude , -z        |    [m]\n"
		<< "  --flying-speed , -v0          |    [m/s]\n"
		<< " LiDAR parameters:\n"
		<< "  --angular-speed , -omega      |    [rot/s]\n"
		<< "  --angle-from-nadir , -theta   |    [deg]\n"
		<< "  --pulse-frequency , -fp       |    [Hz]\n"
		<< " NOISE parameters:\n"
		<< "  --planimetric-mean , -muXY    |    [m]\n"
		<< "  --planimetric-std , -sigmaXY  |    [m]\n"
		<< "  --altimetric-mean , -muZ      |    [m]\n"
		<< "  --altimetric-std , -sigmaZ    |    [m]\n"
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
	bool noiseFree = false;
	bool perfectNormals = false;
	bool useOC = false;
	bool verbose = false;

	// Flight parameters:
	double alphaX = 0.5; // by default, fly along xMid
	double altitude = 1000; // [m] 150 --> 1000 for a low-altitude LiDAR
	double v0 = 60; // [m/s]
	double omega = 150; // [rotations/s]
	double theta = 20; // [deg] angle from nadir
	double freq = 400000; // [Hz]

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
		} else if (std::string(argv[i]) == "--alpha-X" || std::string(argv[i]) == "-alphaX"){
			alphaX = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--flying-altitude" || std::string(argv[i]) == "-z"){
			altitude = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--flying-speed" || std::string(argv[i]) == "-v0"){
			v0 = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--angular-speed" || std::string(argv[i]) == "-omega"){
			omega = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--angle-from-nadir" || std::string(argv[i]) == "-theta"){
			theta = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--pulse-frequency" || std::string(argv[i]) == "-fp"){
			freq = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--noise-free" || std::string(argv[i]) == "-nf"){
			noiseFree = true;
		} else if (std::string(argv[i]) == "--perfect-normals" || std::string(argv[i]) == "-pn"){
			perfectNormals = true;
		} else if (std::string(argv[i]) == "--use-opt-ctrs" || std::string(argv[i]) == "-useOC"){
			useOC = true;
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
		if (noiseFree){
			std::cout << " - Noise-free scan" << std::endl;
		} else {
			std::cout << " - Noise parameters:\n"
					<<"    > Planimetric error ~ N("<<muXY<<", "<<sigmaXY<<"^2)\n"
					<<"    > Altimetric error ~ N("<<muZ<<", "<<sigmaZ<<"^2)" << std::endl;
		}
		std::cout << " - alphaX: " << alphaX << std::endl;
		std::cout << " - Flying altitude: " << altitude << " m" << std::endl;
		std::cout << " - Flying speed: " << v0 << " m/s" << std::endl;
		std::cout << " - Angular speed (azimuthal): " << omega << " rotations/s ("
										  << omega * 360 << " deg/s)" << std::endl;
		std::cout << " - Angle from nadir (theta): " << theta << " deg" << std::endl;
		std::cout << " - Pulse frequency: " << freq << " Hz" << std::endl;
	}

	Mesh mesh = read_mesh<Mesh,Point>(inFileName.c_str(), verbose);

	// Starting and ending positions
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	// double xMid = ( bbox.xmin() + bbox.xmax() ) / 2;
	double xAB = (1-alphaX)*bbox.xmin() + alphaX*bbox.xmax();
	double deltaY = bbox.ymax() - bbox.ymin();
	double yA = bbox.ymin() - 2 * deltaY;
	double yB = bbox.ymax() + 2 * deltaY;

	Point A(xAB, yA, altitude);
	Point B(xAB, yB, altitude);

	// Output file names
	std::string outFileOC = outBaseName + std::string("_OptCtr") + ext;
	std::string outFileN = outBaseName + std::string("_normals") + ext;
	std::string outFileP = outBaseName + std::string("_pts") + ext;

	// Unit conversions:
	omega *= 2*M_PI; // from [rot/s] to [rad/s]
	theta *= M_PI / 180; // from [deg] to [rad]

	// Change convention for theta: becomes the polar angle
	theta = M_PI - theta;
	double phi_0 = 0;

	// Acquisition with Optical Centers + Normals
	int outProperty = 12;
	Point_set pcdAll = elliptical_aerial_lidar(mesh, A, B, v0, omega, theta, phi_0, freq, outProperty, verbose);

	// normal noise
	if (!noiseFree) add_normal_noise(pcdAll, muXY, sigmaXY, muZ, sigmaZ, verbose);

	// points only
	Point_set pcdPts = pcdAll;
	remove_optical_centers(pcdPts);
	pcdPts.remove_normal_map();
	write_point_set(outFileP.c_str(), pcdPts, verbose);

	// normals		
	Point_set pcdN;
	if (perfectNormals){ // perfect normals
		pcdN = pcdAll;
		remove_optical_centers(pcdN);
	} else { // estimated normals
		int k = 20;
		if (useOC){ // use optical centres to orient normals
			pcdN = compute_and_orient_normals_based_on_origin(pcdAll, k, verbose);
		} else { // don't use OC for orientation
			pcdN = pcdPts;
			compute_and_orient_normals(pcdN, k, verbose);
		}
		
	}
	write_point_set(outFileN.c_str(), pcdN, verbose);

	// optical centers
	Point_set pcdOC = pcdAll;
	pcdOC.remove_normal_map();
	write_point_set(outFileOC.c_str(), pcdOC, verbose);

	return 0;
}
