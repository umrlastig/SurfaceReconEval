#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cout << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
	} else {

		// Extraction of command-line parameters
		const char* samplingFile = argv[1];
		const char* PFile = argv[2];
		double alpha = std::atof(argv[3]);
		const char* outFile = argv[4];

		// Read data
		Point_set pds = readPointSet(samplingFile);
		Point_set P = readPointSet(PFile);

		// Apply processing
		remove_points_too_far_from_P(pds, P, alpha);

		// Write filtered point set to file
		writePointSet(outFile, pds);
	}
	return 0;
}