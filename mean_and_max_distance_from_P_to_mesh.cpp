#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 4)
	{
		std::cout << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
	} else {

		// Extraction of command-line parameters
		const char* pcdFile = argv[1];
		const char* meshFile = argv[2];
		const char* outFile = argv[3];

		// Read data
		Point_set pcd = readPointSet(pcdFile);
		Mesh mesh = readMesh(meshFile);

		// Compute distances
		mean_and_max_distance_from_P_to_mesh(pcd, mesh);

		// Write distances to file
		// writePointSet(outFile, pds);
	}
	return 0;
}