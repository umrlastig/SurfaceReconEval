#include "evaluation.hpp"

void display_help(){
	std::cout << "\n========= COMPUTE MESH_ALPHA =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input-mesh, -iM       |    File containing the mesh you wish to process\n"
		<< " --point-cloud, -P       |    File containing the point cloud P to which distances must be computed\n"
		<< " --output-mesh, -o       |    File where to write mesh_alpha\n"
		<< " --alpha, -a             |    Value of alpha parameter\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --verbose, -v           |    Display information throughout execution\n"
		<< " --debug, -d             |    Display additional information: index of each input sample, nearest neighbor, distance to it and whether it is removed\n"
		<< " --help, -h              |    Display this information\n"
		<< std::endl;
}

int main(int argc, char** argv)
{
	// Default values for cmd-line parameters:
	std::string meshFile = "";
	std::string pcdFile = "";
	std::string outFile = "";
	double alpha = -1;
	bool verbose = false;
	bool debug = false;

	// Extraction of command-line parameters
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help();
			return 0;
		} else if (std::string(argv[i]) == "--input-mesh" || std::string(argv[i]) == "-iM"){
			meshFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--point-cloud" || std::string(argv[i]) == "-P"){
			pcdFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output-mesh" || std::string(argv[i]) == "-o"){
			outFile = argv[++i];
		} else if (std::string(argv[i]) == "--alpha" || std::string(argv[i]) == "-a"){
			alpha = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else if (std::string(argv[i]) == "--debug" || std::string(argv[i]) == "-d"){
			debug = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help();
			return 1;
		}
	}

	// Check if mandatory parameters were set
	if (meshFile == "" || pcdFile == "" || outFile == "" || alpha == -1){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}
	if (verbose) std::cout << "verbose flag active" << std::endl;
	if (debug) std::cout << "debug mode active" << std::endl;

	// Read data
	Mesh mesh = read_mesh<Mesh,Point>(meshFile.c_str(), verbose);
	Point_set pcd = read_point_set<Point_set>(pcdFile.c_str(), verbose);

	// Compute mesh_alpha
	bool geodesic = false;
	Mesh meshAlpha = mesh_alpha(mesh, pcd, alpha, geodesic, verbose, debug);

	// Write mesh_alpha to file
	write_mesh(outFile.c_str(), meshAlpha, verbose);

	return 0;
}
