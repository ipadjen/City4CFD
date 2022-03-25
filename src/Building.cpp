#include "Building.h"

#include "geomutils.h"
#include "LoD12.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

Building::Building()
        : PolyFeature(1), _height(-g_largnum) {}

Building::Building(const int internalID)
        : PolyFeature(1, internalID), _height(-g_largnum) {}

Building::Building(const nlohmann::json& poly)
        : PolyFeature(poly, 1), _height(-g_largnum) {}

Building::Building(const nlohmann::json& poly, const int internalID)
        : PolyFeature(poly, 1, internalID), _height(-g_largnum) {}

Building::~Building() = default;

void Building::check_feature_scope(const Polygon_2& otherPoly) {
        for (auto& poly: _poly.rings()) {
            for (auto& vert : poly) {
                if (geomutils::point_in_poly(vert, otherPoly))
                    return;
            }
        }
//    std::cout << "Poly ID " << this->get_id() << " is outside the influ region. Deactivating." << std::endl;
    this->deactivate();
}

void Building::refine() {
    typedef Mesh::Halfedge_index           halfedge_descriptor;
    typedef Mesh::Edge_index               edge_descriptor;

    double target_edge_length = 5;//5;
    unsigned int nb_iter =  30;//30;

    PMP::remove_self_intersections(_mesh);

    //-- Set the property map for constrained edges
    Mesh::Property_map<edge_descriptor,bool> is_constrained =
            _mesh.add_property_map<edge_descriptor,bool>("e:is_constrained",false).first;

    //-- Detect sharp features
    for (auto& e : edges(_mesh)) {
        halfedge_descriptor hd = halfedge(e,_mesh);
        if (!is_border(e,_mesh)) {
            double angle = CGAL::Mesh_3::dihedral_angle(_mesh.point(source(hd,_mesh)),
                                                        _mesh.point(target(hd,_mesh)),
                                                        _mesh.point(target(next(hd,_mesh),_mesh)),
                                                        _mesh.point(target(next(opposite(hd,_mesh),_mesh),_mesh)));
            if (CGAL::abs(angle)<179.5)
                is_constrained[e]=true;
        }
    }

    PMP::isotropic_remeshing(faces(_mesh), target_edge_length, _mesh,
                             PMP::parameters::number_of_iterations(nb_iter)
                                     .edge_is_constrained_map(is_constrained));

//    PMP::remove_self_intersections(_mesh);
}

void Building::translate_footprint(const double h) {
    for (auto& ring : _base_heights) {
        for (auto& pt : ring) {
            pt += h;
        }
    }
}

double Building::max_dim() {
    std::vector<double> dims;
    EPICK::Vector_2 diag(_poly.bbox().xmax() - _poly.bbox().xmin(), _poly.bbox().ymax() - _poly.bbox().ymin());

    dims.emplace_back(diag.squared_length() * pow(cos(M_PI_4), 2));
    dims.emplace_back(diag.squared_length() * pow(sin(M_PI_4), 2));
    dims.emplace_back(_height * _height);

    return sqrt(*(std::max_element(dims.begin(), dims.end())));
}

double Building::get_height() const {
    return _height;
}

void Building::get_cityjson_info(nlohmann::json& b) const {
    b["type"] = "Building";
//  b["attributes"];
//    get_cityjson_attributes(b, _attributes);
//    float hbase = z_to_float(this->get_height_base());
//    float h = z_to_float(this->get_height());
//    b["attributes"]["TerrainHeight"] = _baseHeights.back(); // temp - will calculate avg for every footprint
    b["attributes"]["measuredHeight"] = _height - geomutils::avg(_base_heights[0]);
}

void Building::get_cityjson_semantics(nlohmann::json& g) const { // Temp for checking CGAL mesh properties
    Face_property semantics;
    bool foundProperty;
    boost::tie(semantics, foundProperty) = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    //   auto semantics = _mesh.property_map<face_descriptor, std::string>("f:semantics");
    if (!foundProperty) throw std::runtime_error("Semantic property map not found!");

    std::unordered_map<std::string, int> surfaceId;
    surfaceId["RoofSurface"]   = 0; g["semantics"]["surfaces"][0]["type"] = "RoofSurface";
    surfaceId["GroundSurface"] = 1; g["semantics"]["surfaces"][1]["type"] = "GroundSurface";
    surfaceId["WallSurface"]   = 2; g["semantics"]["surfaces"][2]["type"] = "WallSurface";

    for (auto& faceIdx : _mesh.faces()) {
        auto it = surfaceId.find(semantics[faceIdx]);
        if (it == surfaceId.end()) throw std::runtime_error("Could not find semantic attribute!");

        g["semantics"]["values"][faceIdx.idx()] = it->second;
    }
}

std::string Building::get_cityjson_primitive() const {
    return "MultiSurface";
};

TopoClass Building::get_class() const {
    return BUILDING;
}

std::string Building::get_class_name() const {
    return "Building";
}