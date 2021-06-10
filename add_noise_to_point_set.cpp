#include "raytracing.hpp"

/*
	Commodity executable to add planimetric and altimetric normal noise to
	a given point cloud without having to re-run a LiDAR simulation.
*/

int main(int argc, char* argv[])
{
	if (argc != 8)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< "Arguments should be passed in order following:"
		<< "$1: input point cloud file"
		<< "$2: output point cloud file"
		<< "$3: muXY"
		<< "$4: sigmaXY"
		<< "$5: muZ"
		<< "$6: sigmaZ"
		<< "$7: verbose (\"1\" to activate verbose mode)"
		<< std::endl;
		return 1;
	}

	const char* inFileName = argv[1];
	const char* outFileName = argv[2];

	double muXY = std::atof(argv[3]);
	double sigmaXY = std::atof(argv[4]);
	double muZ = std::atof(argv[5]);
	double sigmaZ = std::atof(argv[6]);

	const std::string verb = argv[7];

	// Set verbose flag
	bool verbose = (verb == "1") ? true : false;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	Point_set pcd = read_point_set<Point_set>(inFileName, verbose);
	add_normal_noise(pcd, muXY, sigmaXY, muZ, sigmaZ, verbose);
	write_point_set(outFileName, pcd, verbose);

	return 0;
}
