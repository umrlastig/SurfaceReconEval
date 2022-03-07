#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include "IO.hpp"

// kernel
typedef CGAL::Simple_cartesian<double> K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;

/*
	Commodity executable to merge several point clouds into one
*/

void display_help(){
	std::cout << "\n========= MERGE POINT CLOUDS =========\n\n";

	std::cout << "First string must be the name of the output file.\n"
			  << "All the followings must be the name of the files containing the point clouds you wish to merge"
	<< std::endl;
}

int main(int argc, char** argv)
{

	// Default values for cmd-line parameters:
	Point_set global_pcd;
	std::string inFile = "";
	std::string outFile = std::string(argv[1]);
	bool verbose = true;
	if (verbose) std::cout << "verbose flag active" << std::endl;

	if (outFile == "-h" || outFile == "--help"){
		display_help();
		return 0;
	}

	// Extraction of command-line parameters
	for (int i = 2; i < argc; ++i) {
		inFile = std::string(argv[i]);
		std::cout << inFile << std::endl;
		Point_set curr_pcd = read_point_set<Point_set>(inFile.c_str(), verbose);
		global_pcd.copy_properties(curr_pcd);
		if ( !global_pcd.join(curr_pcd) ){
			std::cerr << "ERROR while trying to join point set." << std::endl;
		}
	}

	// Write data
	write_point_set(outFile.c_str(), global_pcd, verbose);

	return 0;
}
