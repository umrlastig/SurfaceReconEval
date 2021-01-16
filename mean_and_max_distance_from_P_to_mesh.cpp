#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	// Extraction of command-line parameters
	const char* pcdFile = argv[1];
	const char* meshFile = argv[2];
	const char* outFile = argv[3];
	const std::string verb = argv[4];

	// Set verbose flag
	bool verbose = (verb == "1") ? true : false;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Point_set pcd = read_point_set<Point_set>(pcdFile, verbose);
	Mesh mesh = read_mesh<Mesh,Point>(meshFile, verbose);

	// Compute distances
	mean_and_max_distance_from_P_to_mesh(pcd, mesh);

	return 0;
}
