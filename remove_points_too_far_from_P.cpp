#include "evaluation.hpp"

void display_help(){
	std::cout << "\n========= REMOVE POINTS TOO FAR FROM P =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input-sampling, -iS   |    File containing the sampling you wish to process\n"
		<< " --input-P, -iP          |    File containing the point cloud P to which distances must be computed\n"
		<< " --output-file, -o       |    File where to write processed point cloud (subset of the input sampling)\n"
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
	std::string samplingFile = "";
	std::string PFile = "";
	std::string outFile = "";
	double alpha = -1;
	bool verbose = false;
	bool debug = false;

	// Extraction of command-line parameters
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help();
			return 0;
		} else if (std::string(argv[i]) == "--input-sampling" || std::string(argv[i]) == "-iS"){
			samplingFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--input-P" || std::string(argv[i]) == "-iP"){
			PFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output-file" || std::string(argv[i]) == "-o"){
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
	if (samplingFile == "" || PFile == "" || outFile == "" || alpha == -1){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}
	if (verbose) std::cout << "verbose flag active" << std::endl;
	if (debug) std::cout << "debug mode active" << std::endl;

	// Read data
	Point_set pds = read_point_set<Point_set>(samplingFile.c_str(), verbose);
	Point_set P = read_point_set<Point_set>(PFile.c_str(), verbose);

	// Apply processing
	Point_set processedPcd = remove_points_too_far_from_P(pds, P, alpha, verbose, debug);

	// Write filtered point set to file
	write_point_set(outFile.c_str(), processedPcd, verbose);

	return 0;
}
