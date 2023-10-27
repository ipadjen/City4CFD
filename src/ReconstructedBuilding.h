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

#ifndef CITY4CFD_RECONSTRUCTEDBUILDING_H
#define CITY4CFD_RECONSTRUCTEDBUILDING_H

#include "Building.h"

class ReconstructedBuilding : public Building {
public:
    ReconstructedBuilding();
    ReconstructedBuilding(const int internalID);
    ReconstructedBuilding(const Mesh& mesh);
//    ReconstructedBuilding(const nlohmann::json& poly);
    ReconstructedBuilding(const nlohmann::json& poly, const int internalID);
    ReconstructedBuilding(const Polygon_with_attr& poly, const int internalID);
    ReconstructedBuilding(const std::shared_ptr<ImportedBuilding>& importedBuilding);
    ~ReconstructedBuilding();

    virtual double get_elevation() override;
    virtual void   reconstruct() override;
    virtual void   reconstruct_flat_terrain() override;
    virtual void   get_cityjson_info(nlohmann::json& b) const override;
    virtual void   get_cityjson_semantics(nlohmann::json& g) const override;

protected:
    double _attributeHeight;
    bool   _attributeHeightAdvantage;

    void reconstruct_from_attribute();
    bool reconstruct_again_from_attribute(const std::string& reason);
};


#endif //CITY4CFD_RECONSTRUCTEDBUILDING_H