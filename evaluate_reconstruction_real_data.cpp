#include "evaluation_real_data.hpp"

// Input / Output
#include "IO.hpp"

void display_help(char* argv[]){
	std::cout << "\n========= REAL DATA-BASED EVALUATION =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --in-gtPcd, -iGT              |    File containing the ground truth point cloud\n"
		<< " --in-recon, -iS               |    File containing reconstructed surface mesh to evaluate\n"
		<< " --out-pcd, -o                 |    File where to store the quality-based coloured point cloud\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< "  --max-distance, -dmax        |    [m] (Threshold distance at which to evaluate the reconstruction)\n\n"
		<< "  --verbose, -v                |    Flag to display information along the execution\n"
		<< "  --debug, -d                  |    Flag to display additional information along the execution\n"
		<< std::endl;
}

int main(int argc, char* argv[])
{
	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();

	// Default values for cmd-line parameters:
	std::string gtPcdFile = "";
	std::string reconMeshFile = "";
	std::string outPcdFile = "";
	double maxDistance = 0.2; // [m]
	bool verbose = false;
	bool debug = false;



	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help(argv);
			return 0;
		} else if (std::string(argv[i]) == "--in-gtPcd" || std::string(argv[i]) == "-iGT"){
			gtPcdFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--in-recon" || std::string(argv[i]) == "-iS"){
			reconMeshFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--out-pcd" || std::string(argv[i]) == "-o"){
			outPcdFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--max-distance" || std::string(argv[i]) == "-dmax"){
			maxDistance = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--debug" || std::string(argv[i]) == "-d"){
			debug = true;
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			display_help(argv);
			return 1;
		}
	}

	// Check if in/out files were set
	if (gtPcdFile == "" || reconMeshFile == "" || outPcdFile == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help(argv);
		return 1;
	}

	if (verbose){
		std::cout << "verbose flag active\n" << std::endl;

		std::cout << "Evaluating reconstruction with parameters:" << std::endl;
		std::cout << " - Input GT point cloud: '" << gtPcdFile << "'" << std::endl;
		std::cout << " - Evaluated reconstructed mesh: '" << reconMeshFile << "'" << std::endl;
		std::cout << " - Output point cloud (visualisation): '" << outPcdFile << "'" << std::endl;
		std::cout << " - Threshold distance: " << maxDistance << " m\n" << std::endl;
	}

	Point_set gtPcd = read_point_set<Point_set>(gtPcdFile.c_str(), verbose);
	Mesh reconMesh = read_mesh<Mesh,Point>(reconMeshFile.c_str(), verbose);

	Point_set outPcd = eval_raytracing(reconMesh, gtPcd, maxDistance, debug);

	write_point_set(outPcdFile.c_str(), outPcd, verbose);

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();
	std::cout << "\nTOTAL TIME ELAPSED = "
		<< std::chrono::duration_cast<std::chrono::milliseconds> (toc - tic).count() * 0.001
		<< "[s]" << std::endl;

	return 0;
}
