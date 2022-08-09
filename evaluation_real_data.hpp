#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

// Geometrical objects
#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>

// Accelerating data structure
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

// Kernel
#include <CGAL/Simple_cartesian.h>

// kernel
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;

// geometrical objects
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Point_set_3<Point> Point_set;

// iterators and descriptors
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

// accelerating data structure
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

// optical centers
typedef Point_set::Property_map<double> X_Origin_Map;
typedef Point_set::Property_map<double> Y_Origin_Map;
typedef Point_set::Property_map<double> Z_Origin_Map;

typedef std::tuple<Point*, Point*, double> Adjacency;
typedef std::vector < Adjacency > Adjacency_Graph;

typedef std::tuple<Point*, double, double, bool> Pt_d2Ori_d2Gt_beforeGt;


// neighbor searching
typedef CGAL::Search_traits_3<K> TreeTraitsNS;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraitsNS> Neighbor_search;
typedef Neighbor_search::Tree TreeNS;

// accelerating data structure
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Tree_traits;
typedef CGAL::AABB_tree<Tree_traits> Tree;

// quality
typedef Point_set::Property_map<double> vertexQuality;

// color
typedef Point_set::Property_map<unsigned char> color_map;


struct Skip {
  face_descriptor fd;
  Skip(const face_descriptor fd)
    : fd(fd)
  {}
  bool operator()(const face_descriptor& t) const
  { if(t == fd){
      std::cerr << "ignore " << t  <<std::endl;
    };
    return(t == fd);
  }
};


bool compare_adjacencies(const Adjacency &a, const Adjacency &b){
	/*
		usage: define comparison of adjacencies based on distance
	*/
	return ( std::get<2>(a) ) < ( std::get<2>(b) );
}

bool compare_intersections(const std::pair<Point*, double> &a,
						   const std::pair<Point*, double> &b){
	/*
		usage: define comparison of adjacencies based on distance
	*/
	return ( std::get<1>(a) ) < ( std::get<1>(b) );
}

bool compare_intersections_wrt_Ori(const Pt_d2Ori_d2Gt_beforeGt &a,
								   const Pt_d2Ori_d2Gt_beforeGt &b){
	/*
		usage: define comparison based on distance to Ori
	*/
	return ( std::get<1>(a) ) < ( std::get<1>(b) );
}

bool compare_intersections_wrt_Gt(const Pt_d2Ori_d2Gt_beforeGt &a,
								  const Pt_d2Ori_d2Gt_beforeGt &b){
	/*
		usage: define comparison based on distance to Gt
	*/
	return ( std::get<2>(a) ) < ( std::get<2>(b) );
}

unsigned int increasing_cmap(double d, double dmax){
	if (d > dmax) d = dmax;
	return (unsigned char) round( 255 * ( d/dmax ) );
}

unsigned char decreasing_cmap(double d, double dmax){
	if (d > dmax) d = dmax;
	return (unsigned char) round( 255 * ( 1 - d/dmax ) );
}

double avg(std::vector<double> const& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

void hist(std::vector<double> distances, double a, double b, const int n, int nbOfRays){
	std::cout << "\n*-----------*" << std::endl;
	std::cout << "* HISTOGRAM *" << std::endl;
	std::cout << "*-----------*" << std::endl;

	// Initialize histogram:
	int histo[n];
	for (int i=0; i<n; i++){
		histo[i]=0;
	}
	std::vector<double> mid_distances;
	for (std::vector<double>::const_iterator vi = distances.begin(); vi != distances.end(); ++vi){
		for (int k = 0; k <= n-1; ++k){
			double left = a + k*(b-a)/n;
			double right = left + (b-a)/n;
			if ( *vi >= left && *vi < right ){
				histo[k] += 1;
			}
		}
		if (*vi == b) histo[n-1] += 1;
	}

	// Display full histogram:
	for (int k = 0; k <= n-1; ++k){
		double left = a + k*(b-a)/n;
		double right = left + (b-a)/n;
		std::cout << "I_"<<k<<": ["<<left<<", "<<right<<"] | pop: "<<histo[k] << std::endl;
	}
	
	// Display vector of mid-distances (y-axis)
	std::cout << "vector of mid-distances: [ ";
	for (int k = 0; k < n-1; ++k){
		double left = a + k*(b-a)/n;
		double right = left + (b-a)/n;
		std::cout << (left+right)/2 << ", ";
	}
	std::cout << (2*a + (b-a)*(2*n-1)/n) / 2 << " ]" << std::endl;

	// Display vector of population (y-axis)
	std::cout << "vector of populations: [ ";
	for (int k = 0; k < n-1; ++k){
		std::cout << histo[k] << ", ";
	}
	std::cout << histo[n-1] << " ]" << std::endl;
	std::cout << "\nTOTAL NUMBER OF RAYS: " << nbOfRays
		<<"\n(Divide histogram populations by this number to get the normed histogram)"
		<< std::endl << "*-----------*" << std::endl;
}

void read_optical_centers(Point_set &pcd, X_Origin_Map &x_origin, Y_Origin_Map &y_origin, Z_Origin_Map &z_origin){
	bool success_x = false; bool success_y = false; bool success_z = false;
	boost::tie(x_origin, success_x) = pcd.property_map<double>("x_origin");
	boost::tie(y_origin, success_y) = pcd.property_map<double>("y_origin");
	boost::tie(z_origin, success_z) = pcd.property_map<double>("z_origin");
	if (!success_x || !success_y || !success_z){
		std::cerr << "ERROR: optical centers absent from in:pcd";
	}
}

void add_color_maps(Point_set &pcd, color_map &red, color_map &green, color_map &blue){
	bool successR = false, successG = false, successB = false;
	boost::tie (red, successR) = pcd.add_property_map<unsigned char>("red", 0);
	boost::tie (green, successG) = pcd.add_property_map<unsigned char>("green", 0);
	boost::tie (blue, successB) = pcd.add_property_map<unsigned char>("blue", 0);
	if (!successR || !successG || !successB) {
		std::cerr << "ERROR: impossible to add color property to point set" << std::endl;
	} else {
		std::cout << "    Successfully added color property to point sets" << std::endl;
	}
}

void compute_all_intersections(Tree &tree, Ray &ray, std::vector<Ray_intersection> &intersections,
				std::chrono::duration<double> &duration
){

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	tree.all_intersections( ray, std::back_inserter(intersections) );

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	duration += (toc - tic);
}

void browse_and_insert(
				std::vector<Ray_intersection> &intersections,
				std::vector< Pt_d2Ori_d2Gt_beforeGt > &inters_and_distances,
				Point &Ori,
				Point &Gt,
				std::chrono::duration<double> &duration
){

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	// Browse the n potential intersections and insert them into [inters_and_dist_to_Ori]:
	for (std::vector<Ray_intersection>::iterator interIt = intersections.begin(); interIt != intersections.end(); ++interIt){
		Ray_intersection &intersection = *interIt;
		if(intersection){
			if(boost::get<Point>(&(intersection->first))){
				Point& I = boost::get<Point>(intersection->first);
				double OriI = CGAL::sqrt(CGAL::squared_distance(Ori, I));
				double GtI = CGAL::sqrt(CGAL::squared_distance(I, Gt));
				bool IisBeforeGt = (OriI <= CGAL::sqrt(CGAL::squared_distance(Ori, Gt)));
				inters_and_distances.push_back( Pt_d2Ori_d2Gt_beforeGt(&I, OriI, GtI, IisBeforeGt) );
			}
		}
	}
	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	duration += (toc - tic);
}

void sort_intersections(
					std::vector< Pt_d2Ori_d2Gt_beforeGt > &inters_and_distances,
					std::chrono::duration<double> &duration
){
	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	// Sort the vector of Point* and distances by distance to Origin:
	std::sort(inters_and_distances.begin(), inters_and_distances.end(), compare_intersections_wrt_Ori);

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	duration += (toc - tic);
}

void print_intersections(std::vector< Pt_d2Ori_d2Gt_beforeGt > &inters_and_distances){
	for (std::vector<Pt_d2Ori_d2Gt_beforeGt>::iterator it = inters_and_distances.begin(); it != inters_and_distances.end(); ++it){
		Point a = *std::get<0>(*it);
		double b = std::get<1>(*it);
		double c = std::get<2>(*it);
		double d = std::get<3>(*it);
		std::cout<<"Sorted\n"<< " - "<<"Point: ["<< a <<"] | distance to Ori: "<<b<<"] | distance to Gt: "<<c<<"] | I before Gt: "<<d<<std::endl;
	}		
}

void update_metrics_FP_pcd(
				std::vector< Pt_d2Ori_d2Gt_beforeGt > &inters_and_distances,
				std::vector<double> &distances, // distances of all nearest intersection points between GT ray and reconstructed surface
				int &FalPosFront, // total number of points from recon in front of the closest point from the GT point
				int &FalPosClosest, // number of points of GT for which ray champfer is above a given threshold
				int &oddIntersections, // number of rays for which the closest intersection with the reconS is odd
				int &intersectedRays, // number of rays that have at least one intersection
				double &maxDistance,
				bool &debug,
				std::chrono::duration<double> &duration,
				Point_set &FP_pcd, // False Positives point cloud
				bool export_FP_pts
){

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	if ( inters_and_distances.size() == 0 ){ // No intersection for this ray
		if (debug) std::cout << "No intersection found!" << std::endl;
	} else {
		intersectedRays++;
		// Compute nearest intersection from GT point:
		std::vector< Pt_d2Ori_d2Gt_beforeGt >::iterator itMin = std::min_element(inters_and_distances.begin(),
																		inters_and_distances.end(),
																		compare_intersections_wrt_Gt);
		double rayDist = std::get<2>(*itMin); // distance to Gt
		bool beforeGt = std::get<3>(*itMin); // I is before Gt
		int min_pos = std::distance(inters_and_distances.begin(), itMin); // number of points in front of nearest intersection wrt Gt
		if (export_FP_pts){
			for (std::vector< Pt_d2Ori_d2Gt_beforeGt >::iterator it_FP = inters_and_distances.begin(); it_FP != itMin; ++it_FP){
				FP_pcd.insert(*std::get<0>(*it_FP)); // points in front of the closest intersection
			}
		}
		FalPosFront += min_pos;
		if ( (min_pos + 1) % 2 != 0 ) oddIntersections++; // +1 since min_pos=0 for the first intersection (which is odd)
		if (rayDist >= maxDistance){
			if (beforeGt){
				FalPosClosest++;
				if (export_FP_pts) FP_pcd.insert(*std::get<0>(*itMin)); // closest intersections before Gt and for which distance >= threshold
			}
		} else {
			distances.push_back(rayDist);
		}
		if (debug) std::cout << "Min element is number "<<min_pos<<": "<< rayDist << std::endl;
	}

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	duration += (toc - tic);
}

void visual_evaluation(
				std::vector< Pt_d2Ori_d2Gt_beforeGt > &inters_and_distances,
				double &OriGt,
				double &maxDistance,
				Point_set &outPcd,
				color_map &red,
				color_map &green,
				color_map &blue,
				std::chrono::duration<double> &duration
){

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	// Browse all the intersections, in sorted-order:
	int k2 = 1;
	for (std::vector< Pt_d2Ori_d2Gt_beforeGt >::iterator it = inters_and_distances.begin(); it != inters_and_distances.end(); ++it){
		Point& I = *std::get<0>(*it);
		double& OriI = std::get<1>(*it);
		double GtI = std::get<2>(*it);
		Point_set::iterator itOut = outPcd.insert( I );
		if (OriI < OriGt){ // intersection is BEFORE GT point
			// Point_set::iterator itOut = outPcd.insert( I );
			if (k2 % 2 == 0){ // even intersection
				// YELLOW
				red[*itOut] = increasing_cmap(GtI, maxDistance);
				green[*itOut] = 255;
				// blue[*itOut] = 0; // true by default
			} else { // odd intersection
				// RED
				red[*itOut] = increasing_cmap(GtI, maxDistance);
				green[*itOut] = decreasing_cmap(GtI, maxDistance);
				// blue[*itOut] = 0; // true by default
			}
		} else {  // intersection is AFTER GT point
			if (k2 % 2 == 0){ // even intersection
				// BLUE
				// red[*itOut] = 0; // true by default
				green[*itOut] = decreasing_cmap(GtI, maxDistance);
				blue[*itOut] = increasing_cmap(GtI, maxDistance);
			} else { // odd intersection
				// CYAN
				// red[*itOut] = 0; // true by default
				green[*itOut] = 255;
				blue[*itOut] = increasing_cmap(GtI, maxDistance);
			}
		}
		k2++;
	}

	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	duration += (toc - tic);
}

void print_metrics_values(
					std::vector<double> &distances, // distances of all nearest intersection points between GT ray and reconstructed surface
					int &FalPosFront, // total number of points from recon in front of the closest point from the GT point
					int &FalPosClosest, // number of points of GT for which ray champfer is above a given threshold
					int &oddIntersections, // number of rays for which the closest intersection with the reconS is odd
					int &intersectedRays, // number of rays that have at least one intersection
					int gtPcd_nb_of_pts,
					double &maxDistance, // threshold distance
					int &nIntervals){ // for histogram

	double a = 0, b = maxDistance; hist(distances, a, b, nIntervals, gtPcd_nb_of_pts); // histogram
	std::cout << "    => Mean distance (saturated at " << maxDistance << " [m]): " << avg(distances)
		<< " (over " << distances.size() << " points)" << std::endl;
	std::cout << "    => Nb of pts in front of the closest intersection: "<< FalPosFront << std::endl;
	std::cout << "    => Nb of closest intersections before Gt and for which ray champfer is above maxDistance: "<< FalPosClosest << std::endl;
	std::cout << "    => ratio of odd (~correct) intersections: "<< (double)oddIntersections / (double)intersectedRays * 100 <<"%"<< std::endl;
	std::cout << "True positives      ; False Positives (alg)        ; TP_gt + FP_gt         ; mean Distance"<< std::endl;
	std::cout <<  distances.size() <<";"<< FalPosFront + FalPosClosest <<";"<< gtPcd_nb_of_pts <<";"<< avg(distances)<< std::endl;
}

void print_time_profile(
					std::chrono::duration<double> &DURATION_compute_all_intersections,
					std::chrono::duration<double> &DURATION_browse_and_insert,
					std::chrono::duration<double> &DURATION_sort_intersections,
					std::chrono::duration<double> &DURATION_update_metrics,
					std::chrono::duration<double> &DURATION_visual_evaluation,
					std::chrono::duration<double> &DURATION_ALL
){
	std::cout << "\nDURATION for:"
		<< "\n - 'compute_all_intersections' : " << DURATION_compute_all_intersections.count() << "[s] (" << DURATION_compute_all_intersections.count() / DURATION_ALL.count() * 100 << " %)"
		<< "\n - 'browse_and_insert' : "         << DURATION_browse_and_insert.count()         << "[s] (" << DURATION_browse_and_insert.count() / DURATION_ALL.count() * 100 << " %)"
		<< "\n - 'sort_intersections' : "        << DURATION_sort_intersections.count()        << "[s] (" << DURATION_sort_intersections.count() / DURATION_ALL.count() * 100 << " %)"
		<< "\n - 'update_metrics' : "            << DURATION_update_metrics.count()            << "[s] (" << DURATION_update_metrics.count() / DURATION_ALL.count() * 100 << " %)"
		<< "\n - 'visual_evaluation' : "         << DURATION_visual_evaluation.count()         << "[s] (" << DURATION_visual_evaluation.count() / DURATION_ALL.count() * 100 << " %)"
		<< "\n -  => TOTAL : "                   << DURATION_ALL.count()                       << "[s]"
		<< std::endl;
}

void eval_raytracing_sequential(Mesh &reconMesh,
								Point_set &gtPcd,
								double maxDistance,
								Point_set &outPcd, color_map &red, color_map &green, color_map &blue,
								Point_set &FP_pcd, 
								std::vector<double> &distances,
								int &FalPosFront,
								int &FalPosClosest,
								int &oddIntersections,
								int &intersectedRays,
								bool visual_eval,
								bool export_FP_pts,
								bool debug
						  ){
	std::cout << "\nEVALUATING RECONSTRUCTION" << std::endl;

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now(); // TIC

	// Profiling study:
	std::chrono::duration<double> DURATION_compute_all_intersections = std::chrono::seconds::zero();
	std::chrono::duration<double> DURATION_browse_and_insert         = std::chrono::seconds::zero();
	std::chrono::duration<double> DURATION_sort_intersections        = std::chrono::seconds::zero();
	std::chrono::duration<double> DURATION_update_metrics            = std::chrono::seconds::zero();
	std::chrono::duration<double> DURATION_visual_evaluation         = std::chrono::seconds::zero();
	std::chrono::duration<double> DURATION_ALL                       = std::chrono::seconds::zero();
	
	// Read optical centers and check validity:
	X_Origin_Map x_origin; Y_Origin_Map y_origin; Z_Origin_Map z_origin;
	read_optical_centers(gtPcd, x_origin, y_origin, z_origin);

	
	Tree tree(faces(reconMesh).first, faces(reconMesh).second, reconMesh);
	for (Point_set::const_iterator it = gtPcd.begin(); it != gtPcd.end(); ++it){
		Point Ori = Point(x_origin[*it], y_origin[*it], z_origin[*it]); // origin
		Point Gt = gtPcd.point(*it);
		Ray ray(Ori, Gt);
		double OriGt = CGAL::sqrt(CGAL::squared_distance(Ori, Gt));
		
		std::vector<Ray_intersection> intersections;
		std::vector< Pt_d2Ori_d2Gt_beforeGt > inters_and_distances;

		if (debug){
			std::cout<<"\n\nTrying to compute new intersection(s)..."
			<<"\nRay: "<<ray
			<<"\nOri: "<<Ori
			<<"\nGt: "<<Gt
			<<std::endl;
		}

		compute_all_intersections(tree, ray, intersections, DURATION_compute_all_intersections);

		if (debug) std::cout<<"Done with computing new intersection(s)."<<std::endl;
		if (debug){
			if (intersections.size() != 0) std::cout << "\nintersection(s) and distance to [" << Gt << "]"  << std::endl;
		}

		// Browse the n potential intersections and insert them into [inters_and_dist_to_Ori]:
		browse_and_insert(intersections, inters_and_distances, Ori, Gt, DURATION_browse_and_insert);

		// Sort the vector of Point* and distances by distance to Origin:
		sort_intersections(inters_and_distances, DURATION_sort_intersections);

		// Print the sorted vector of n potential intersections:
		if (debug) print_intersections(inters_and_distances);

		update_metrics_FP_pcd(
				inters_and_distances,
				distances,
				FalPosFront,
				FalPosClosest,
				oddIntersections,
				intersectedRays,
				maxDistance,
				debug,
				DURATION_update_metrics,
				FP_pcd,	export_FP_pts);

		// Visual evaluation:
		if (visual_eval) visual_evaluation(inters_and_distances, OriGt, maxDistance, outPcd, red, green, blue, DURATION_visual_evaluation);
		if (debug) std::cout<<"Done with current intersection(s)"<<std::endl;
	}
	
	std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now(); // TOC
	DURATION_ALL = (toc - tic);

	// print_time_profile(
	// 				DURATION_compute_all_intersections,
	// 				DURATION_browse_and_insert,
	// 				DURATION_sort_intersections,
	// 				DURATION_update_metrics,
	// 				DURATION_visual_evaluation,
	// 				DURATION_ALL);
}

