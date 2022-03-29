/*
  Copyright (c) 2020-2021,
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

#include "TopoFeature.h"

//-- TopoFeature class
TopoFeature::TopoFeature()
        : _mesh(), _id(), _f_active(true), _f_imported(false), _outputLayerID(-1) {}

TopoFeature::TopoFeature(std::string pid)
        : _mesh(), _id(std::move(pid)), _f_active(true), _f_imported(false), _outputLayerID(-1) {}

TopoFeature::TopoFeature(int outputLayerID)
        : _mesh(), _id(), _f_active(true), _f_imported(false), _outputLayerID(outputLayerID) {
    if (_outputLayerID  >= _numOfOutputLayers) _numOfOutputLayers = _outputLayerID + 1;
}

TopoFeature::~TopoFeature() = default;

int TopoFeature::_numOfOutputLayers = 0;

const int TopoFeature::get_internal_id() const {
    return -1;
}

int TopoFeature::get_num_output_layers() {
    return _numOfOutputLayers;
}

Mesh& TopoFeature::get_mesh() {
    return _mesh;
}

const Mesh& TopoFeature::get_mesh() const {
    return _mesh;
}

void TopoFeature::set_id(unsigned long id) {
    _id = std::to_string(id);
}

std::string TopoFeature::get_id() const {
    return _id;
}

const int TopoFeature::get_output_layer_id() const {
    return _outputLayerID;
}

bool TopoFeature::is_active() const {
    return _f_active;
}

bool TopoFeature::is_imported() const {
    return _f_imported;
}

void TopoFeature::deactivate() {
    _f_active = false;
}

void TopoFeature::get_cityjson_info(nlohmann::json& j) const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
}

void TopoFeature::get_cityjson_semantics(nlohmann::json& g) const {
    // TEMP until I figure what to do with this
}


std::string TopoFeature::get_cityjson_primitive() const {
    //TEMP UNTIL ALL FUNCTIONS ARE IMPLEMENTED
    return "Nope";
}