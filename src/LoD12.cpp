/*
  City4CFD
 
  Copyright (c) 2021-2023, 3D Geoinformation Research Group, TU Delft

  This file is part of City4CFD.

  City4CFD is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  City4CFD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with City4CFD.  If not, see <http://www.gnu.org/licenses/>.

  For any information or further details about the use of City4CFD, contact
  Ivan Pađen
  <i.paden@tudelft.nl>
  3D Geoinformation Research Group
  Delft University of Technology
*/

#include "LoD12.h"

#include "geomutils.h"

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<std::vector<double>>& base_elevations)
        : _elevation(), _poly(poly), _baseElevations(base_elevations) {}

LoD12::LoD12(const Polygon_with_holes_2& poly,
             const std::vector<std::vector<double>>& base_elevations,
             const double elevation)
        : _elevation(elevation), _poly(poly), _baseElevations(base_elevations) {}

void LoD12::set_elevation(const double& elevation) {
    _elevation = elevation;
}

void LoD12::reconstruct(Mesh& mesh) {
    mesh.clear();
    // Add semantics with face properties to the property map
    auto surfaceType = mesh.add_property_map<face_descriptor , std::string>("f:semantics", "").first;
    face_descriptor fIdx;

    CDT cdt_buildings;

    //-- Map CDT and Mesh vertices
    std::map<CDT::Vertex_handle, Mesh::Vertex_index> cdtToMesh;

    int polyCount = 0;
    for (auto& poly : _poly.rings()) { // Loop over polys
        std::vector<Vertex_handle> cdt_handle;
        std::vector<Mesh::Vertex_index> mesh_vertex;
        int count = 0;
        for (auto vert = poly.vertices_begin(); vert != poly.vertices_end(); ++vert) { // Loop over poly vertices
            cdt_handle.emplace_back(cdt_buildings.insert(ePoint_3(vert->x(),
                                                                  vert->y(),
                                                                  _baseElevations[polyCount][count])));
            mesh_vertex.emplace_back(mesh.add_vertex(Point_3(vert->x(),
                                                             vert->y(),
                                                             _baseElevations[polyCount][count++])));
            cdtToMesh[cdt_handle.back()] = mesh_vertex.back();
            mesh_vertex.emplace_back(mesh.add_vertex(Point_3(vert->x(),
                                                             vert->y(),
                                                             _elevation)));
        }
        cdt_handle.emplace_back(cdt_handle.front());
        mesh_vertex.emplace_back(mesh_vertex.front());
        mesh_vertex.emplace_back(mesh_vertex.front() + 1);

        //- Add constraints and create mesh faces for sides
        for (auto i = 0; i < cdt_handle.size() - 1; ++i) {
            cdt_buildings.insert_constraint(cdt_handle[i], cdt_handle[i + 1]);

            auto it1 = mesh.vertices_begin();
            auto it2 = mesh.vertices_begin();

            auto v1 = cdtToMesh[cdt_handle[i]];
            auto v2 = cdtToMesh[cdt_handle[i + 1]];

            std::advance(it1, v1.idx() + 1);
            std::advance(it2, v2.idx() + 1);

            fIdx = mesh.add_face(v1, v2, *it1);
            if (fIdx != Mesh::null_face()) {
                surfaceType[fIdx] = "WallSurface";
            }
            fIdx = mesh.add_face(v2, *it2, *it1);
            if (fIdx != Mesh::null_face()) {
                surfaceType[fIdx] = "WallSurface";
            }
        }
        ++polyCount;
    }
    //- Handle top and bottom
    geomutils::mark_domains(cdt_buildings);
    for (auto& it : cdt_buildings.finite_face_handles()) {
        if (!it->info().in_domain()) continue;

        auto it1 = mesh.vertices_begin();
        auto it2 = mesh.vertices_begin();
        auto it3 = mesh.vertices_begin();

        std::advance(it1, cdtToMesh[it->vertex(0)]);
        std::advance(it2, cdtToMesh[it->vertex(1)]);
        std::advance(it3, cdtToMesh[it->vertex(2)]);

        fIdx = mesh.add_face(*it1, *it3, *it2); // Bottom face
        if (fIdx != Mesh::null_face()) {
            surfaceType[fIdx] = "GroundSurface";
        }
        fIdx = mesh.add_face(*std::next(it1), *std::next(it2), *std::next(it3));
        if (fIdx != Mesh::null_face()) {
            surfaceType[fIdx] = "RoofSurface";
        }
    }
}