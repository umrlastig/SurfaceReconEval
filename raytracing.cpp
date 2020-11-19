#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

int main(int argc, char* argv[])
{
	//const int x(1); const int y(-4); const int z(6);
	//std::cout << x << std::endl;
	Point A(1,-4,6);
	std::cout << A << std::endl;

	// CGAL::Point_3<int> P();
}
