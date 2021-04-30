#include "raytracing.hpp"

void display_help(char* argv[]){
	std::cout << "\n========= NORMALS ESTIMATION USING OPTICAL CENTERS =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input-file, -i              |    File containing the point cloud of which you wish to estimate normals\n"
		<< " --output-file, -o             |    Name of the output file where to store point cloud"
		<< std::endl << std::endl;

	std::cout << "OPTIONAL parameter:\n"
			  << "--------------------\n"
		<< " --nb-neighbors, -k            |    Number of neighbors to use for normal estimation"
		<< " --verbose, -v                 |    Display information throughout execution"
		<< std::endl << std::endl;
}

int main(int argc, char* argv[])
{
	// Default values for cmd-line parameters:
	std::string inFileName = "";
	std::string outFileName = "";
	int k = 10;
	bool verbose = false;

	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help(argv);
			return 0;
		} else if (std::string(argv[i]) == "--input-file" || std::string(argv[i]) == "-i"){
			inFileName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output-file" || std::string(argv[i]) == "-o"){
			outFileName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--nb-neighbors" || std::string(argv[i]) == "-k"){
			k = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help(argv);
			return 1;
		}
	}

	// Check if in/out files were set
	if (inFileName == "" || outFileName == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help(argv);
		return 1;
	}

	if (verbose){
		std::cout << "verbose flag active\n" << std::endl;
		std::cout << " - Input point cloud: '" << inFileName << "'" << std::endl;
		std::cout << " - Output point cloud: '" << outFileName << "'" << std::endl;
		std::cout << " - Number of neighbors: " << k << std::endl << std::endl;
	}

	Point_set pcd = read_point_set<Point_set>(inFileName.c_str(), verbose);
	// compute_and_orient_normals(pcd, k, verbose);
	Point_set outPcd = compute_and_orient_normals_based_on_origin(pcd, k, verbose);
	write_point_set(outFileName.c_str(), outPcd, verbose);

	return 0;
}
