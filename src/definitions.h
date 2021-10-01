#ifndef CITYCFD_DEFINITIONS_H
#define CITYCFD_DEFINITIONS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/IO/WKT.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/interpolation_functions.h>

#include "nlohmann/json.hpp"

//-- CGAL Basics
typedef CGAL::Exact_predicates_inexact_constructions_kernel  EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel    EPECK;
typedef CGAL::Projection_traits_xy_3<EPICK>                  iProjection_traits;
typedef CGAL::Projection_traits_xy_3<EPECK>                  Projection_traits;

//-- Kernel converter
typedef CGAL::Cartesian_converter<EPECK, EPICK> EKtoIK;
typedef CGAL::Cartesian_converter<EPICK, EPECK> IKtoEK;

//-- CGAL Point
typedef EPICK::Point_2            Point_2;
typedef EPICK::Point_3            Point_3;
typedef EPECK::Point_3            ePoint_3;
typedef EPICK::Segment_3          Segment_3;
typedef CGAL::Point_set_3<Point_3> Point_set_3;

//-- CGAL Mesh
typedef CGAL::Surface_mesh<Point_3>                      Mesh;
typedef Mesh::Vertex_index                               vertex_descriptor;
typedef Mesh::Face_index                                 face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::string> Face_property;

//-- CGAL normal
namespace PMP = CGAL::Polygon_mesh_processing;
typedef EPICK::Vector_3 Vector;

//-- CGAL tree search
typedef CGAL::Search_traits_3<EPICK>                 Traits;
//typedef CGAL::Kd_tree<Traits>                         SearchTree;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>    Neighbor_search;
typedef Neighbor_search::Tree                         SearchTree;
typedef CGAL::Fuzzy_iso_box<Traits>                   Fuzzy_iso_box;
typedef CGAL::Fuzzy_sphere<Traits>                    Fuzzy_sphere;

struct FaceInfo2
{
    FaceInfo2(){}
    int nesting_level;
    bool in_domain(){
        return nesting_level%2 == 1;
    }
    int surfaceLayer = -9999; // Face handle to output mesh for specific surface layer
};
//-- CGAL triangulation
typedef CGAL::Triangulation_vertex_base_with_id_2<Projection_traits>               Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Projection_traits>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Projection_traits, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                TDS;
//typedef CGAL::Exact_predicates_tag                                                  Itag;
typedef CGAL::Exact_intersections_tag                                              Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Projection_traits, TDS, Itag>   CDTp;
typedef CGAL::Constrained_triangulation_plus_2<CDTp>                               CDT;
typedef CDT::Point                                                                 Point;
typedef CDT::Face_handle                                                           Face_handle;
typedef CDT::Vertex_handle                                                         Vertex_handle;

//testing
typedef CGAL::Delaunay_triangulation_2<CGAL::Projection_traits_xy_3<EPICK>> DT;

typedef CGAL::Polygon_2<EPICK>                                                     Polygon_2;
typedef CGAL::Polygon_2<Projection_traits>                                         Polygon_3;
//typedef CGAL::Polygon_with_holes_2<EPICK>                                          Polygon_with_holes_2;

//-- TopoClasses
typedef enum {
    TERRAIN          = 0,
    BUILDING         = 1,
    BOUNDARY         = 2,
    WATER            = 3,
    BRIDGE           = 4,
    ROAD             = 5,
    FOREST           = 6,
    SIDES            = 7,
    TOP              = 8,
    SURFACELAYER     = 9,
} TopoClass;

const std::map<int, std::string> topoClassName {
        {0, "Terrain"},
        {1, "Building"},
        {2, "Boundary"},
        {3, "Terrain"},
        {4, "Bridge"},
        {5, "Road"},
        {6, "Forest"},
        {7, "Sides"},
        {8, "Top"},
        {9, "SurfaceLayer"},
};

//-- Output formats
typedef enum {
    OBJ       = 0,
    CityJSON  = 1,
    STL       = 2,
} OutputFormat;

//-- CGAL's Polygon_with_holes container sucks for not having iterator over all rings
struct Polygon_with_holes_2 {
    std::vector<Polygon_2> _rings;

    std::vector<Polygon_2>& rings() {return _rings;}
    const std::vector<Polygon_2>& rings() const {return _rings;}

    bool has_holes() const {
        if (_rings.size() > 1) return true;
        return false;
    }
    Polygon_2& outer_boundary() {
        return _rings.front();
    }
    const Polygon_2& outer_boundary() const {
        return _rings.front();
    }

    const auto holes_begin() const {
        if (has_holes()) return _rings.begin() + 1; else return _rings.end();
    }
    const auto holes_end() const {
        return _rings.end();
    }

    const auto bbox() const {return _rings.front().bbox();}
};

//-- Global constants
const double infty     = 1e6;
const double smallnum  = 1e-6;

#endif //CITYCFD_DEFINITIONS_H