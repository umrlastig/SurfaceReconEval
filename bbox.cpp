#include "raytracing.hpp"

int main(int argc, char* argv[])
{
	std::string inFileName = std::string(argv[1]);

	if (inFileName == ""){
		std::cerr << "Error: no input file provided" << std::endl;
		return 1;
	}

	bool verbose = true;
	Mesh mesh = read_mesh<Mesh,Point>(inFileName.c_str(), verbose);
	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	double xMin=bbox.xmin(), yMin=bbox.ymin(), zMin=bbox.zmin();
	double xMax=bbox.xmax(), yMax=bbox.ymax(), zMax=bbox.zmax();

	std::cout << "\n========== Bounding box dimensions ==========\n" << std::endl;
	std::cout << " Xmin: " << xMin << std::endl;
	std::cout << " Xmax: " << xMax << std::endl;
	std::cout << " Ymin: " << yMin << std::endl;
	std::cout << " Ymax: " << yMax << std::endl;
	std::cout << " Zmin: " << zMin << std::endl;
	std::cout << " Zmax: " << zMax << std::endl << std::endl;

	std::cout << " Delta_X: " << xMax - xMin << std::endl;
	std::cout << " Delta_Y: " << yMax - yMin << std::endl;
	std::cout << " Delta_Z: " << zMax - zMin << std::endl << std::endl;

	std::cout << " Scene area: " << (xMax - xMin) * (yMax - yMin) << std::endl;

	std::cout << "\n=============================================" << std::endl;
	return 0;
}
