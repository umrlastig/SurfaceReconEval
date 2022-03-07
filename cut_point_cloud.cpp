#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
// #include "IO.hpp"
#include "raytracing.hpp"
#include <algorithm>

// kernel
typedef CGAL::Simple_cartesian<double> K;

// geometrical objects
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

/*
	Commodity executable to cut a given point
	cloud into several sub point clouds, specifying
	the maximum size of the resulting chunks
*/

void display_help(){
	std::cout << "\n========= CUT POINT CLOUD =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input , -i              |     Input point cloud\n"
		<< " --output , -o             |     Output base name for output files\n"
		<< " --number-of-points , -n   |     Maximum number of points per chunk\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --verbose , -v            |     Display information throughout execution\n"
		<< " --help , -h               |     Display this information"
		<< std::endl << std::endl;
}

int main(int argc, char** argv)
{

	// Default values for cmd-line parameters:
	std::string inFile = "";
	std::string outBase = "";
	int nbPtsPerChunk = -1;
	bool verbose = false;

	// Extraction of command-line parameters
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help();
			return 0;
		} else if (std::string(argv[i]) == "--input" || std::string(argv[i]) == "-i"){
			inFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output" || std::string(argv[i]) == "-o"){
			outBase = argv[++i];
		} else if (std::string(argv[i]) == "--number-of-points" || std::string(argv[i]) == "-n"){
			nbPtsPerChunk = std::atoi(argv[++i]);
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help();
			return 1;
		}
	}

	// Check if mandatory parameters were set
	if (inFile == "" || outBase == "" || nbPtsPerChunk == -1){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}

	if (verbose) std::cout << "verbose flag active" << std::endl;

	// Read data
	Point_set inPcd = read_point_set<Point_set>(inFile.c_str(), verbose);
	X_Origin_Map in_x_origin; Y_Origin_Map in_y_origin; Z_Origin_Map in_z_origin;
	read_optical_centers(inPcd, in_x_origin, in_y_origin, in_z_origin);

	int nbChunks = inPcd.number_of_points() / nbPtsPerChunk;

	// Generate all "full" (actually nbPtsPerChunk-sized) chunks
	for (int iChunk = 0; iChunk <= nbChunks; iChunk++){
		// Initialize new point set:
		X_Origin_Map out_x_origin; Y_Origin_Map out_y_origin; Z_Origin_Map out_z_origin;
		int OCProperty = 2;
		Point_set outPcd = initialize_point_set(OCProperty, out_x_origin, out_y_origin, out_z_origin, false);

		// Browse points to insert in current chunck:
		int indexBeginningChunk = iChunk * nbPtsPerChunk;
		int indexEndChunk       = std::min( (iChunk + 1) * nbPtsPerChunk, (int)inPcd.number_of_points() );
		for (Point_set::const_iterator 	pi =  inPcd.begin() + indexBeginningChunk;
										pi != inPcd.begin() + indexEndChunk;
										++pi)
		{
			// Insert point with OC:
			Point_set::iterator it = outPcd.insert(inPcd.point(*pi));
			out_x_origin[*it] = in_x_origin[*pi];
			out_y_origin[*it] = in_y_origin[*pi];
			out_z_origin[*it] = in_z_origin[*pi];
		}

		// Write current sub point cloud to file:
		std::string outFile = outBase + std::to_string(iChunk) + std::string(".ply");
		write_point_set(outFile.c_str(), outPcd, verbose);
	}

	return 0;
}
