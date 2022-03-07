#include "raytracing.hpp"
#include <algorithm>

// validity code: 0 --> invalid
typedef Point_set::Property_map<unsigned char> validityCode;


void sub_point_set_with_origin(std::string inFile, std::string outFile, int ratio, Point Ori){
	/*
		This function just keeps one point every [ratio] if it's valid
	*/
	Point_set inPcd = read_point_set<Point_set>(inFile.c_str(), true);
	Point_set outPcd;

	// Add origin properties:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie (x_origin, success_x) = outPcd.template add_property_map<double>("x_origin", 0);
	boost::tie (y_origin, success_y) = outPcd.template add_property_map<double>("y_origin", 0);
	boost::tie (z_origin, success_z) = outPcd.template add_property_map<double>("z_origin", 0);
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: impossible to add optical center properties";
	}

	// Read validity codes:
	validityCode valCodes; bool success_c = false;
	boost::tie (valCodes, success_c) = inPcd.property_map<unsigned char>("code");
	if (!success_c) std::cerr << "ERROR: impossible to read validity codes property";

	int k = 0; // to know if point is kept or not
	for (typename Point_set::const_iterator pi = inPcd.begin(); pi != inPcd.end(); ++pi){

		if (valCodes[*pi] != 0){ // if point is valid, otherwise just go to the next one
			if (k % ratio == 0){
				typename Point_set::iterator it = outPcd.insert(inPcd.point(*pi));
				x_origin[*it] = Ori.x();
				y_origin[*it] = Ori.y();
				z_origin[*it] = Ori.z();
			}
			k++;
		}
	}
	write_point_set(outFile.c_str(), outPcd, true);
}

void sub_point_set_with_origin_matrix_format(std::string inFile, std::string outFile, int height, int ratio, Point Ori){
	/*
		This function keeps one point every [ratio] in both dimensions of the matrix format, if that point happens to be valid
	*/
	Point_set inPcd = read_point_set<Point_set>(inFile.c_str(), true);
	Point_set outPcd;
	int nbInvalidPts = 0;

	// Add origin properties:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie (x_origin, success_x) = outPcd.template add_property_map<double>("x_origin", 0);
	boost::tie (y_origin, success_y) = outPcd.template add_property_map<double>("y_origin", 0);
	boost::tie (z_origin, success_z) = outPcd.template add_property_map<double>("z_origin", 0);
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: impossible to add optical center properties";
	}

	// Read validity codes:
	validityCode valCodes; bool success_c = false;
	boost::tie (valCodes, success_c) = inPcd.property_map<unsigned char>("code");
	if (!success_c) std::cerr << "ERROR: impossible to read validity codes property";

	int i_row = 0, i_col = 0; // to know if point is kept or not
	for (typename Point_set::const_iterator pi = inPcd.begin(); pi != inPcd.end(); ++pi){
		if (i_row == height){ // we reached the end of the current row
			i_col++;
			i_row=0;
		}
		if (i_col % ratio == 0){ // row is eligible
			if (i_row % ratio == 0){ // colomn is eligible
				if ((int)(valCodes[*pi]) != 0){ // if point is valid, otherwise just go to the next one*/
					typename Point_set::iterator it = outPcd.insert(inPcd.point(*pi));
					x_origin[*it] = Ori.x();
					y_origin[*it] = Ori.y();
					z_origin[*it] = Ori.z();
				} else {
					nbInvalidPts++;
				}
			}
		}
		i_row++;
	}
	std::cout << "Number of invalid points: " << nbInvalidPts << std::endl;
	write_point_set(outFile.c_str(), outPcd, true);
}

void display_help(){
	std::cout << "\n========= MATRIX SUB-SAMPLING =========\n\n";

	std::cout << "MANDATORY parameters:\n"
			  << "---------------------\n"
		<< " --input , -i              |     Input point cloud\n"
		<< " --output , -o             |     Output base name for output files\n"
		<< " --height , -h             |     Height of the matrix format\n"
		<< " --ratio , -r              |     Sampling ratio along both matrix dimensions: 2 => 1/2 row and 1/2 column\n"
		<< " --OriginXYZ , -oXYZ       |     Sensor position for all points: Ox Oy Oz\n"
		<< std::endl;

	std::cout << "OPTIONAL parameters:\n"
			  << "--------------------\n"
		<< " --verbose , -v            |     Display information throughout execution\n"
		<< " --help , -h               |     Display this information"
		<< std::endl << std::endl;
}

int main(int argc, char* argv[])
{
	// Default values for cmd-line parameters:
	std::string inFileName = "";
	std::string outFileName = "";
	int ratio = -1;
	int height = -1; // matrix format
	bool verbose = false;
	double Ox, Oy, Oz;

	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){
			display_help(); return 1;
		} else if (std::string(argv[i]) == "--input" || std::string(argv[i]) == "-i"){
			inFileName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--output" || std::string(argv[i]) == "-o"){
			outFileName = std::string(argv[++i]);
		} else if (std::string(argv[i]) == "--height" || std::string(argv[i]) == "-h"){
			height = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--ratio" || std::string(argv[i]) == "-r"){
			ratio = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--OriginXYZ" || std::string(argv[i]) == "-oXYZ"){
			Ox = std::atof(argv[++i]);
			Oy = std::atof(argv[++i]);
			Oz = std::atof(argv[++i]);
		} else if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v"){
			verbose = true;
		} else {
			std::cerr << "Invalid option: '" << argv[i] << "'" << std::endl;
			return 1;
		}
	}

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();

	// Check if mandatory parameters were set
	if (inFileName == "" || outFileName == "" || ratio == -1 || height == -1){
		std::cerr << "Error: one of mandatory parameters has not been set." << std::endl;
		display_help();
		return 1;
	}

	if (verbose) std::cout << "verbose flag active" << std::endl;

	Point Ori(Ox, Oy, Oz);
	sub_point_set_with_origin_matrix_format(inFileName, outFileName, height, ratio, Ori);

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();
	std::cout << "TOTAL TIME ELAPSED = "
		<< std::chrono::duration_cast<std::chrono::milliseconds> (toc - tic).count() * 0.001
		<< "[s]" << std::endl;

	return 0;
}
