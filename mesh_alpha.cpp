#include "evaluation.hpp"


int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cout << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
	} else {
		// Extraction of command-line parameters
		const char* meshFile = argv[1];
		const char* lidarFile = argv[2];
		double alpha = std::atof(argv[3]);
		const char* outFile = argv[4];

		// Read data
		Mesh mesh = readMesh(meshFile);
		std::ifstream is_lidar (lidarFile);
		Point_set pcd; CGAL::read_ply_point_set(is_lidar, pcd);

		// Compute mesh_alpha
		bool geodesic = false; bool verbose = false;
		Mesh mesh_alpha = Mesh_alpha(mesh, pcd, alpha, geodesic, verbose);

		// Write mesh_alpha to file
		writeMesh(outFile, mesh_alpha);
	}
	return 0;
}
