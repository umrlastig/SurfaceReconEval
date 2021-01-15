#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 6)
	{
		std::cout << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
	} else {
		// Extraction of command-line parameters
		const char* meshFile = argv[1];
		const char* lidarFile = argv[2];
		double alpha = std::atof(argv[3]);
		const char* outFile = argv[4];
		const char* verb = argv[5];

		// Set verbose flag
		bool verbose = (verb == "1") ? true : false;

		// Read data
		Mesh mesh = readMesh(meshFile, verbose);
		Point_set pcd = readPointSet(lidarFile, verbose);

		// Compute mesh_alpha
		bool geodesic = false; bool debug=false;
		Mesh meshAlpha = Mesh_alpha(mesh, pcd, alpha, geodesic, verbose, debug);

		// Write mesh_alpha to file
		writeMesh(outFile, meshAlpha, verbose);
	}
	return 0;
}
