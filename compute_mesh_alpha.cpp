#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 6)
	{
		std::cout << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	// Extraction of command-line parameters
	const char* meshFile = argv[1];
	const char* lidarFile = argv[2];
	double alpha = std::atof(argv[3]);
	const char* outFile = argv[4];
	const std::string verb = argv[5];

	// Set verbose flag
	bool verbose = (verb == "1") ? true : false;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Mesh mesh = read_mesh<Mesh,Point>(meshFile, verbose);
	Point_set pcd = read_point_set<Point_set>(lidarFile, verbose);

	// Compute mesh_alpha
	bool geodesic = false; bool debug=false;
	Mesh meshAlpha = mesh_alpha(mesh, pcd, alpha, geodesic, verbose, debug);

	// Write mesh_alpha to file
	write_mesh(outFile, meshAlpha, verbose);

	return 0;
}
