#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 6)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
	} else {

		// Extraction of command-line parameters
		const char* samplingFile = argv[1];
		const char* PFile = argv[2];
		double alpha = std::atof(argv[3]);
		const char* outFile = argv[4];
		const char* verb = argv[5];

		// Set verbose flag
		bool verbose = (verb == "1") ? true : false;

		// Read data
		Point_set pds = readPointSet(samplingFile, verbose);
		Point_set P = readPointSet(PFile, verbose);

		// Apply processing
		remove_points_too_far_from_P(pds, P, alpha, verbose);

		// Write filtered point set to file
		writePointSet(outFile, pds, verbose);
	}
	return 0;
}