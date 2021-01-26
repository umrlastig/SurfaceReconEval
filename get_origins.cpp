#include "raytracing.hpp"

/*
	Commodity executable to implement getOrigins() function.
	Reads the point cloud of input file and writes the resulting
	point cloud containing sensor positions.
*/

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	const char* inFileName = argv[1];
	const char* outFileName = argv[2];

	bool verbose = true;

	Point_set pcdOC = read_point_set<Point_set>(inFileName, verbose);

	Point_set pcdRays = getOrigins(pcdOC);

	write_point_set(outFileName, pcdRays, verbose);	

	return 0;
}
