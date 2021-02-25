#include <iostream>
#include <fstream>
#include <string>

// Input / Output
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/OBJ_reader.h>

// Polygon_mesh_processing
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;

std::string getFileExt(const char* fileName) {
	const std::string s(fileName);
	size_t i = s.rfind('.', s.length());
	if (i != std::string::npos) {
		return(s.substr(i+1, s.length() - i));
	}
	return("");
}


template <typename Mesh, typename Point>
Mesh read_OBJ_mesh(const char* fileName, bool verbose){
	Mesh mesh;

	// Variables to store polygon soup:
	std::ifstream is (fileName);
	std::vector<Point> points;
	std::vector<std::vector<std::size_t>> polygons;

	if (!CGAL::read_OBJ(is, points, polygons)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
	int nbPtsIni = points.size();

	// Processing polygon soup:
	if (!PMP::is_polygon_soup_a_polygon_mesh(polygons))
	{
		// polygons do not define a valid mesh
		if (verbose) std::cout << "Warning: '" << fileName << "' does not define a valid polygon mesh " << std::endl;
		if ( !PMP::orient_polygon_soup(points, polygons) ){ // try to orient
			// orientation failed
			if (verbose) std::cout << "Warning: orientation of polygon soup from '" << fileName << "' failed\n => "
				<< points.size() - nbPtsIni << " points were duplicated" << std::endl;
			// if (verbose) std::cout << " => Number of points: " << points.size() << std::endl;
			// if (verbose) std::cout << " => Number of faces: " << polygons.size() << std::endl;
		} else {
			// orientation succeeded
			if (verbose) std::cout << "Orientation of the polygon soup [" << fileName << "] was a success." << std::endl;
		}
	} else {
		// polygons define a valid mesh
		if (verbose) std::cout << "Success: [" << fileName << "] defines a valid polygon mesh " << std::endl;
	}

	PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

	return mesh;
}

template <typename Mesh, typename Point>
Mesh read_mesh(const char* fileName, bool verbose){
	if (verbose) std::cout << "\nREADING MESH: " << fileName << std::endl;
	Mesh mesh; std::ifstream is (fileName);

	// Determine extension then extract points and polygons
	std::string extension = getFileExt(fileName);
	if (extension == "obj"){
		mesh = read_OBJ_mesh<Mesh,Point>(fileName, verbose);
	} else if (extension == "ply") {
		if (!read_ply(is, mesh)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
		// std::vector<Point> points;
		// std::vector<std::vector<std::size_t>> polygons;
		// std::vector<CGAL::Color> fcolors;
		// std::vector<CGAL::Color> vcolors;
		// bool success = CGAL::read_PLY (is, points, polygons, fcolors, vcolors);
		// if (success) std::cout << "success" << std::endl;
		// std::cout << "number of points: " << points.size() << std::endl;
		// std::cout << "number of polygons: " << polygons.size() << std::endl;
		// PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
	} else if (extension == "off") {
		if (!read_off(is, mesh)) std::cerr << "\nERROR: impossible to read mesh\n" << std::endl;
	} else {
		std::cout << "\nERROR: extension [." << extension << "] not readable.\n" << std::endl;
	}
	if (verbose) std::cout << "    Number of vertices: " << mesh.number_of_vertices() << std::endl;
	if (verbose) std::cout << "    Number of faces: " << mesh.number_of_faces() << std::endl;
	return mesh;
}

template <typename Point_set>
Point_set read_point_set(const char* fileName, bool verbose){
	if (verbose) std::cout << "READING POINT SET: " << fileName << std::endl;
	Point_set pcd; std::ifstream is (fileName);

	std::string extension = getFileExt(fileName);
	if (extension == "ply") {
		if (!CGAL::read_ply_point_set(is, pcd)) std::cerr
			<< "\nERROR: impossible to read point set\n" << std::endl;
	} else if (extension == "off") {
		if (!CGAL::read_off_point_set(is, pcd)) std::cerr
			<< "\nERROR: impossible to read mesh\n" << std::endl;
	} else {
		std::cerr << "\nERROR: extension [." << extension << "] not readable.\n" << std::endl;
	}
	if (verbose) std::cout << "    Number of points: " << pcd.number_of_points() << std::endl;
	return pcd;
}

template <typename Mesh>
void write_mesh(const char* fileName, Mesh &mesh, bool verbose){
	std::ofstream of(fileName);
	of.precision(15);
	bool success = false;

	// Determine extension then use appropriate function
	std::string extension = getFileExt(fileName);
	if (extension == "ply"){
		success = write_ply(of, mesh);
	} else if (extension == "off") {
		success = write_off(of, mesh);
	} else {
		std::cerr << "ERROR: extension [." << extension << "] not writable." << std::endl;
	}
	// Outcome:
	if (success){
		if (verbose) std::cout << "\nwrote: '" << fileName
			<< "' (" << mesh.number_of_vertices() << " vertices / "
			<< mesh.number_of_faces() << " faces)" << std::endl;
	} else {
		std::cerr << "ERROR: impossible to write '" << fileName << "'" << std::endl;
	}
}

template <typename Point_set>
void write_point_set(const char* fileName, Point_set &pcd, bool verbose){
	std::ofstream of(fileName);
	of.precision(15);
	bool success = false;

	// Determine extension then use appropriate function
	std::string extension = getFileExt(fileName);
	if (extension == "ply"){
		success = CGAL::write_ply_point_set(of, pcd);
	} else if (extension == "off") {
		success = CGAL::write_off_point_set(of, pcd);
	} else {
		std::cerr << "ERROR: extension [." << extension << "] not writable." << std::endl;
	}
	// Outcome:
	if (success){
		if (verbose) std::cout << "\nwrote: '" << fileName
			<< "' (" << pcd.number_of_points() << " points)" << std::endl;
	} else {
		std::cerr << "ERROR: impossible to write '" << fileName << "'" << std::endl;
	}
}

