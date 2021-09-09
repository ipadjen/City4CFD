#include "Boundary.h"

// TODO: put more thought into this constructor once the time comes
Boundary::Boundary(const ConfigData& configData)
    :TopoFeature(), _outerPts() {
    if (configData.dimOfDomain != -infty) {
        _dimOfDomain = configData.dimOfDomain;
    }

    if (configData.topHeight != -infty) {
        _topHeight = configData.topHeight;
    }
}

//-- Deactivate point cloud points that are out of bounds
void Boundary::set_bounds_to_pc(Point_set_3& pointCloud) const { // Will try to template it to include CDT
    //- 80% of the total domain size. The rest is left for the buffer zone
    auto it = pointCloud.points().begin();
    int count = 0;
    while (it != pointCloud.points().end()) {
        if (!geomtools::point_in_circle(*it, pointOfInterest, 0.8 * _dimOfDomain)) {
            pointCloud.remove(pointCloud.begin() + count);
        } else {
            ++it;
            ++count;
        }
    }
    pointCloud.collect_garbage(); // Free removed points from the memory
}

//-- Add buffer between the terrain and domain edge
void Boundary::add_buffer(Point_set_3& pointCloud) {
    const int nPts      = 360; // Hardcoded
    const double angInt = 2 * M_PI / (double)nPts;
    double ang = 0;
    for (auto i = 0; i < nPts; ++i) {
        double xPt = pointOfInterest.x() + _dimOfDomain * cos(ang + angInt);
        double yPt = pointOfInterest.y() + _dimOfDomain * sin(ang + angInt);
        ang = ang + angInt;
        Point_3 pt(xPt, yPt, 0.0); // Height set at 0 for now. Could be average of the edge or average of the whole terrain
        _outerPts.push_back(pt); // for top and sides
        pointCloud.insert(pt);
    }
    _outerPts.push_back(_outerPts[0]); // Put first point at the end to close the loop
}

void Boundary::threeDfy() {
    //-- Create mesh out of this building
    std::vector<Mesh::vertex_index> mesh_vertex_side;
    std::vector<Mesh::vertex_index> mesh_vertex_top;

    //-- Top needs DT, sides are done manually
    CDT cdt_top;
    for (auto &pt : _outerPts) {
        cdt_top.insert(Point_3(pt.x(), pt.y(), _topHeight));
    }

    int count = 0;
    //-- Add mesh vertices and store them in a vector
    for (auto it = _outerPts.begin(); it != _outerPts.end(); ++it) {
        mesh_vertex_side.emplace_back(_mesh.add_vertex(*it));
        mesh_vertex_side.emplace_back(_mesh.add_vertex(Point_3(it->x(), it->y(), _topHeight)));
    }

    //-- Add middle top point to mesh
//    Mesh::vertex_index topMiddlePoint = _meshTop.add_vertex(Point_3(bndInfo.xcent, bndInfo.ycent, bndInfo.height));

    //-- Add mesh faces for side
    for (auto i = 0; i < mesh_vertex_side.size() - 3; i= i + 2) {
        // -- i + 1 is i lifted up
        int v1 = i;
        int v2 = i + 2;

        _mesh.add_face(mesh_vertex_side[v2], mesh_vertex_side[v1], mesh_vertex_side[v1 + 1]);
        _mesh.add_face(mesh_vertex_side[v2 + 1], mesh_vertex_side[v2], mesh_vertex_side[v1 + 1]);
    }

    //-- Add mesh faces for top
    geomtools::cdt_to_mesh(cdt_top, _meshTop);
}

const Mesh& Boundary::get_top_mesh() const {
    return _meshTop;
}

TopoClass Boundary::get_class() const {
    return BOUNDARY;
}
