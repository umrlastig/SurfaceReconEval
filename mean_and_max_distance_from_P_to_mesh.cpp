#include "evaluation.hpp"

void display_help(){
	std::cout << "\n========= MEAN AND MAX DISTANCE FROM P TO MESH =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --point-cloud, -P     |     Input mesh\n"
		<< " --mesh, -M            |     Output mesh\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --verbose, -v         |     Display information throughout execution\n"
		<< " --help, -h :          |     Display this information"
		<< std::endl << std::endl;
}

int main(int argc, char** argv)
{
	// Default values for cmd-line parameters:
	std::string pcdFile = "";
	std::string meshFile = "";
	bool verbose = false;

	// Extraction of command-line parameters
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help();
			return 0;
		} else if (std::string(argv[i]) == "--point-cloud" || std::string(argv[i]) == "-P"){
			pcdFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--mesh" || std::string(argv[i]) == "-M"){
			meshFile = argv[++i];
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help();
			return 1;
		}
	}

	// Check if mandatory parameters were set
	if (pcdFile == "" || meshFile == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}

	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Point_set pcd = read_point_set<Point_set>(pcdFile.c_str(), verbose);
	Mesh mesh = read_mesh<Mesh,Point>(meshFile.c_str(), verbose);

	// Compute distances
	mean_and_max_distance_from_P_to_mesh(pcd, mesh);

	return 0;
}
