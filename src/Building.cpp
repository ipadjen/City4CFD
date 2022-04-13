/*
  Copyright (c) 2021-2022,
  Ivan Pađen <i.paden@tudelft.nl>
  3D Geoinformation,
  Delft University of Technology

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "Building.h"

#include "geomutils.h"
#include "LoD12.h"
#include "Terrain.h"

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/clip.h>

Building::Building()
        : PolyFeature(1), _height(-g_largnum) {}

Building::Building(const int internalID)
        : PolyFeature(1, internalID), _height(-g_largnum) {}

Building::Building(const nlohmann::json& poly)
        : PolyFeature(poly, 1), _height(-g_largnum) {}

Building::Building(const nlohmann::json& poly, const int internalID)
        : PolyFeature(poly, 1, internalID), _height(-g_largnum) {}

Building::~Building() = default;

void Building::clip_bottom(const Terrainptr& terrain) {
    //-- Get terrain subset
    Mesh terrainSubsetMesh = terrain->mesh_subset(_poly);
    PMP::reverse_face_orientations(terrainSubsetMesh);

    //-- Set exact point maps
    Exact_point_map mesh1_exact_points =
            terrainSubsetMesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_point_map mesh2_exact_points =
            _mesh.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;
    Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, terrainSubsetMesh);
    Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, _mesh);

    //-- Mesh processing and clip
    PMP::remove_degenerate_faces(_mesh);
    PMP::remove_degenerate_edges(_mesh);
    geomutils::remove_self_intersections(_mesh); //-- Todo handle for hybrid reconstruction
    PMP::clip(_mesh, terrainSubsetMesh, params::vertex_point_map(mesh2_vpm), params::vertex_point_map(mesh1_vpm));
}

void Building::translate_footprint(const double h) {
    for (auto& ring : _base_heights) {
        for (auto& pt : ring) {
            pt += h;
        }
    }
}

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