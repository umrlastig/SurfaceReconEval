#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 6)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	// Extraction of command-line parameters
	const char* samplingFile = argv[1];
	const char* PFile = argv[2];
	double alpha = std::atof(argv[3]);
	const char* outFile = argv[4];
	const std::string verb = argv[5];

	// Set verbose flag
	bool verbose = (verb == "1") ? true : false;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Point_set pds = read_point_set<Point_set>(samplingFile, verbose);
	Point_set P = read_point_set<Point_set>(PFile, verbose);

	// Apply processing
	remove_points_too_far_from_P(pds, P, alpha, verbose);

	// Write filtered point set to file
	write_point_set(outFile, pds, verbose);

	return 0;
}
