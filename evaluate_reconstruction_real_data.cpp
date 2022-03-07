#include "evaluation_real_data.hpp"

// Input / Output
#include "IO.hpp"

void display_help(char* argv[]){
	std::cout << "\n========= REAL DATA-BASED EVALUATION =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --in-gtPcd , -iGT              |    Files containing the ground truth point cloud (Separated by '***')\n"
		<< " --in-recon , -iS               |    File containing reconstructed surface mesh to evaluate\n"
		<< " --out-pcd , -o                 |    File where to store the quality-based coloured point cloud\n"		
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --max-distance , -dmax         |    [m] (Threshold distance at which to evaluate the reconstruction)\n"
		<< " --visual-eval , -visu          |    Flag to generate a quality-related coloured point cloud\n"
		<< " --verbose , -v                 |    Flag to display information along the execution\n"
		<< " --debug , -d                   |    Flag to display additional information along the execution\n"
		<< std::endl;
}

int main(int argc, char* argv[])
{
	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();

	// Default values for cmd-line parameters:
	std::string gtPcdFiles = "";
	std::string reconMeshFile = "";
	std::string outPcdFile = "";
	double maxDistance = 0.2; // [m]
	bool visual_eval = false;
	bool verbose = false;
	bool debug = false;



	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help(argv);
			return 0;
		} else if (std::string(argv[i]) == "--in-gtPcd" || std::string(argv[i]) == "-iGT"){
			gtPcdFiles = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--in-recon" || std::string(argv[i]) == "-iS"){
			reconMeshFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--out-pcd" || std::string(argv[i]) == "-o"){
			outPcdFile = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--max-distance" || std::string(argv[i]) == "-dmax"){
			maxDistance = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--visual-eval" || std::string(argv[i]) == "-visu"){
			visual_eval = true;
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
	if (gtPcdFiles == "" || reconMeshFile == "" || outPcdFile == ""){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help(argv);
		return 1;
	}

	if (verbose){
		std::cout << "verbose flag active\n" << std::endl;

		std::cout << "Evaluating reconstruction with parameters:" << std::endl;
		std::cout << " - Input GT point cloud: '" << gtPcdFiles << "'" << std::endl;
		std::cout << " - Evaluated reconstructed mesh: '" << reconMeshFile << "'" << std::endl;
		std::cout << " - Output point cloud (visualisation): '" << outPcdFile << "'" << std::endl;
		std::cout << " - Threshold distance: " << maxDistance << " m\n" << std::endl;
	}

	Mesh reconMesh = read_mesh<Mesh,Point>(reconMeshFile.c_str(), verbose);

	// Metrics:
	std::vector<double> distances; // distances of all nearest intersection points between GT ray and reconstructed surface
	int weakPrecision = 0; // total number of points from recon in front of the closest point from the GT point
	int recall = 0; // number of points of GT for which ray champfer is above a given threshold
	int oddIntersections = 0; // number of rays for which the closest intersection with the reconS is odd
	int intersectedRays = 0; // number of rays that have at least one intersection
	int nIntervals = 20; // for histogram
	int nbRays = 0; // for histogram

	// Assessment point clouds:
	Point_set outPcd; // Point cloud to visualise quality
	// Add color maps:
	color_map red, green, blue;
	add_color_maps(outPcd, red, green, blue);

	Point_set WP_pcd; // Point cloud containing all points that belong to 'WeakPrecision' metric

	std::string delimiter = "***";
	int ePos, iPos = 0;
	std::string gtFile;
	do {
		ePos = gtPcdFiles.find(delimiter, iPos);
		if (ePos == std::string::npos) ePos = gtPcdFiles.size();
		gtFile = gtPcdFiles.substr(iPos, ePos-iPos);
		if (gtFile.size() == 0) break;
		std::cout << "gtFile: " << gtFile << "\n";


		Point_set gtPcd = read_point_set<Point_set>(gtFile.c_str(), verbose);

		eval_raytracing_sequential(reconMesh, gtPcd, maxDistance,
			outPcd, red, green, blue,
			WP_pcd,
			distances,
			weakPrecision,
			recall,
			oddIntersections,
			intersectedRays,
			visual_eval,
			debug);

		nbRays+=gtPcd.number_of_points();

		std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();
		std::cout << "\nTOTAL TIME ELAPSED = "
			<< std::chrono::duration_cast<std::chrono::milliseconds> (toc - tic).count() * 0.001
			<< "[s]\n\n";


		iPos = ePos + delimiter.size();
	} while (ePos != gtPcdFiles.size());

	// Print metrics values:
	print_metrics_values(distances, weakPrecision, recall, oddIntersections, intersectedRays, nbRays, maxDistance, nIntervals);

	if (visual_eval) write_point_set(outPcdFile.c_str(), outPcd, verbose);

	// for (face_descriptor f : faces(reconMesh)){
	// 	if (PMP::is_degenerate_triangle_face(f, reconMesh)){
	// 		// std::cout << "Face: " << f << " is degenerate." << std::endl;
	// 		reconMesh.remove_face(f);
	// 		// n_degen++;
	// 	}
	// }

	return 0;
}
