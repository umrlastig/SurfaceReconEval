#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "IO.hpp"

// kernel
typedef CGAL::Simple_cartesian<double> K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

/*
	Commodity executable to convert a mesh from one extension
	to another (or just to check if a mesh can be read by CGAL
	before running the evaluation pipeline, for example)
*/

void display_help(){
	std::cout << "\n========= CONVERT MESH =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input, -i           |     Input mesh\n"
		<< " --output, -o          |     Output mesh\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --verbose, -v         |     Display information throughout execution\n"
		<< " --help, -h :          |     Display this information"
		<< std::endl << std::endl;
}

int main(int argc, char** argv)
{

	// Default values for cmd-line parameters:
	std::string inFile = "";
	std::string outFile = "";
	bool verbose = false;

	// Extraction of command-line parameters
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help();
			return 0;
		} else if (std::string(argv[i]) == "--input" || std::string(argv[i]) == "-i"){
			inFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output" || std::string(argv[i]) == "-o"){
			outFile = argv[++i];
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help();
			return 1;
		}
	}

	// Check if mandatory parameters were set
	if (inFile == "" || outFile == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}

	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Mesh mesh = read_mesh<Mesh,Point>(inFile.c_str(), verbose);

	// Write data
	write_mesh(outFile.c_str(), mesh, verbose);

	return 0;
}
