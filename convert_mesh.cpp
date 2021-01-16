#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "IO.hpp"

// kernel
typedef CGAL::Simple_cartesian<double> K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;


int main(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cerr << "ERROR: illegal number of command line arguments to execute " << argv[0]
		<< std::endl;
		return 1;
	}

	// Extraction of command-line parameters
	const char* inFile = argv[1];
	const char* outFile = argv[2];

	// Set verbose flag
	bool verbose = true;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Mesh mesh = read_mesh<Mesh,Point>(inFile, verbose);

	// Write data
	write_mesh(outFile, mesh, verbose);

	return 0;
}
