/*
  City4CFD
 
  Copyright (c) 2021-2022, 3D Geoinformation Research Group, TU Delft  

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

#include "PointCloud.h"

#include "io.h"
#include "geomutils.h"
#include "Config.h"
#include "PolyFeature.h"

#include <boost/locale.hpp>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>

PointCloud::PointCloud()  = default;
PointCloud::~PointCloud() = default;

void PointCloud::random_thin_pts() {
    if (Config::get().terrainThinning > 0 + global::smallnum) {
        std::cout <<"\nRandomly thinning terrain points" << std::endl;
        _pointCloudTerrain.remove(CGAL::random_simplify_point_set(_pointCloudTerrain,
                                                                  Config::get().terrainThinning),
                                  _pointCloudTerrain.end());
        _pointCloudTerrain.collect_garbage();
        std::cout << "    Terrain points after thinning: " << _pointCloudTerrain.size() << std::endl;
    }
}

void PointCloud::smooth_terrain() {
    typedef CGAL::Parallel_if_available_tag Concurrency_tag;

    std::cout << "Smoothing terrain" << std::endl;

    //-- WLOP simplification and regularization
    double retain_percentage = 100;
    int& maxTerrainPts = Config::get().maxSmoothPts;
    if (maxTerrainPts > 0 && _pointCloudTerrain.size() > maxTerrainPts) {
        retain_percentage = (double)maxTerrainPts / (double)_pointCloudTerrain.size() * 100.;
        std::cout << "    Performing additional (optimized) terrain thinning to " << maxTerrainPts << " points" << std::endl;
    }

    std::cout << "    Smoothing terrain 1/3..." << std::flush;
    const double neighbor_radius = 0.5;   // neighbors size.
    Point_set_3 simplPts;
    CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>
            (_pointCloudTerrain, simplPts.point_back_inserter(),
             CGAL::parameters::select_percentage(retain_percentage).
                     neighbor_radius (neighbor_radius));
    _pointCloudTerrain.clear();

    std::cout << "\r    Smoothing terrain 2/3..." << std::flush;

    //-- Create CDT of current terrain pts
    DT dt;
    dt.insert(simplPts.points().begin(), simplPts.points().end());
    simplPts.clear(); // end of scope for simplePts

    //-- Make a mesh out of DT
    Mesh mesh;
    geomutils::dt_to_mesh(dt, mesh);
    dt.clear(); // end of scope for the dt

    //-- Isotropic remeshing
    const double target_edge_length = 0;
    const unsigned int nb_iter =  10;
    PMP::remove_degenerate_faces(mesh);
    PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
                             PMP::parameters::number_of_iterations(nb_iter)
                                     );

    //-- Smoothing
    std::cout << "\r    Smoothing terrain 3/3..." << std::flush;
    const double time = 1;
    PMP::smooth_shape(mesh, time, CGAL::parameters::number_of_iterations(Config::get().nSmoothIterations));

    std::cout << "\r    Smoothing terrain...done" << std::endl;

    //-- Mesh back to points
//    _pointCloudTerrain.clear();
    for (auto& pt : mesh.points()) {
        _pointCloudTerrain.insert(pt);
    }
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::create_flat_terrain(const PolyFeaturesPtr& lsFeatures) {
    std::cout << "\nCreating flat terrain" << std::endl;
    for (auto& f : lsFeatures) {
        if (f->get_poly().rings().empty()) {
//            std::cout << "Empty polygon?" << std::endl; //todo investigate this
            continue;
        }
        for (auto& pt : f->get_poly().outer_boundary()) {
            _pointCloudTerrain.insert(Point_3(pt.x(), pt.y(), 0.0));
        }
    }
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::set_flat_terrain() {
    Point_set_3 flatPC;
    for (auto& pt : _pointCloudTerrain.points()) {
        flatPC.insert(Point_3(pt.x(), pt.y(), 0.));
    }
    _pointCloudTerrain = flatPC;
    _pointCloudTerrain.add_property_map<bool> ("is_building_point", false);
}

void PointCloud::flatten_polygon_pts(const PolyFeaturesPtr& lsFeatures,
                                     std::vector<EPECK::Segment_3>& constrainedEdges) {
    std::cout << "\n    Flattening surfaces" << std::endl;
    std::map<int, Point_3> flattenedPts;

    //-- Construct a connectivity map and remove duplicates along the way
    auto is_building_pt = _pointCloudTerrain.property_map<bool>("is_building_point").first;
    std::unordered_map<Point_3, int> pointCloudConnectivity;
    auto it = _pointCloudTerrain.points().begin();
    int count = 0;
    while (it != _pointCloudTerrain.points().end()) {
        auto itPC = pointCloudConnectivity.find(*it);
        if (itPC != pointCloudConnectivity.end()) {
            _pointCloudTerrain.remove(_pointCloudTerrain.begin() + count);
        } else {
            pointCloudConnectivity[*it] = count;
            ++it;
            ++count;
        }
    }
    _pointCloudTerrain.collect_garbage();

    //-- Construct search tree from ground points
    SearchTree searchTree(_pointCloudTerrain.points().begin(),
                          _pointCloudTerrain.points().end(),
                          Config::get().searchtree_bucket_size);

    //-- Perform averaging
    PolyFeaturesPtr avgFeatures;
    std::list<CGAL::Polygon_with_holes_2<EPECK>> avgPolys;
    int outputlayerid;
    for (auto& f : lsFeatures) {
        auto ita = Config::get().flattenSurfaces.find(f->get_output_layer_id());
        if (ita != Config::get().flattenSurfaces.end()) {
            //todo temp
            f->mark_flat();
            avgPolys.push_back(f->get_poly().get_exact());
            outputlayerid = f->get_output_layer_id();
//            if (this->flatten_polygon_inner_points(_pointCloudTerrain, flattenedPts, searchTree,
//                                            pointCloudConnectivity, f->get_poly(), f->get_output_layer_id()))
//            avgFeatures.push_back(f);
//            f->flatten_polygon_inner_points(_pointCloudTerrain, flattenedPts, searchTree, pointCloudConnectivity);
//            if (f->is_flat()) avgFeatures.push_back(f);
        }
    }
    //-- Handle adjacencies
    // traverse over all polygons and try to merge
    auto itOuter = avgPolys.begin();
    while (itOuter != avgPolys.end()) {
        auto itInner = avgPolys.begin();
        if (itOuter == itInner) {
            ++itOuter;
            continue;
        }
        CGAL::Polygon_with_holes_2<EPECK> unionR;
        if (CGAL::join(*itOuter, *itInner, unionR)) {
            (*itOuter) = unionR;
            avgPolys.erase(itInner);
        } else ++itInner;
        ++itOuter;
    }
    // grab edges of new polygons
    Converter<EPECK, EPICK> to_inexact;
    for (auto& poly : avgPolys) {
        Polygon_2 ring;
        for (auto& pt : poly.outer_boundary()) ring.push_back(Point_2(to_inexact(pt.x()), to_inexact(pt.y())));
        std::vector<double> elevations;
        geomutils::interpolate_poly_from_pc(ring, elevations, _pointCloudTerrain);
        Polygon_3 constrainRing;
        for (auto i = 0; i < ring.size() - 1; ++i) {
            auto pt1 = ePoint_3(ring[i].x(), ring[i].y(), elevations[i]);
            auto pt2 = ePoint_3(ring[i + 1].x(), ring[i + 1].y(), elevations[i + 1]);
            constrainedEdges.emplace_back(pt1, pt2);
        }
    }
    // flatten points
    for (auto& poly : avgPolys) {
        this->flatten_polygon_inner_points(_pointCloudTerrain, flattenedPts, searchTree,
                                           pointCloudConnectivity, poly, outputlayerid);
    }

    //-- Handle border for flattened polys
    this->buffer_flat_edges(avgFeatures, constrainedEdges);

    //-- Change points with flattened values
    int pcOrigSize = _pointCloudTerrain.points().size();
    for (auto& it : flattenedPts) {
        _pointCloudTerrain.insert(it.second);
    }
    for (int i = 0; i < pcOrigSize; ++i) {
        auto it = flattenedPts.find(i);
        if (it != flattenedPts.end()) {
            _pointCloudTerrain.remove(i);
            flattenedPts.erase(i);
        }
    }
    _pointCloudTerrain.collect_garbage();
}

bool PointCloud::flatten_polygon_inner_points(const Point_set_3& pointCloud,
                                               std::map<int, Point_3>& flattenedPts,
                                               const SearchTree& searchTree,
                                               const std::unordered_map<Point_3, int>& pointCloudConnectivity,
                                               const CGAL::Polygon_with_holes_2<EPECK>& _poly,
                                               const int outputLayerID) {
    std::vector<int>    indices;
    std::vector<double> originalHeights;
    auto is_building_pt = pointCloud.property_map<bool>("is_building_point").first;
    //-- Take tree subset bounded by the polygon
    std::vector<Point_3> subsetPts;
//    Polygon_2 bbox = geomutils::calc_bbox_poly(_poly.rings().front());
    auto bbox = _poly.bbox();
    Point_2 bbox1(bbox.xmin(), bbox.ymin());
    Point_2 bbox2(bbox.xmax(), bbox.ymax());
    Fuzzy_iso_box pts_range(bbox1, bbox2);
    searchTree.search(std::back_inserter(subsetPts), pts_range);

    //-- Collect points that have not been already flattened
    for (auto& pt3 : subsetPts) {
        ePoint_2 pt(pt3.x(), pt3.y());
        if (CGAL::bounded_side_2(_poly.outer_boundary().begin(),
                                 _poly.outer_boundary().end(),
                                 pt) != CGAL::ON_UNBOUNDED_SIDE) {
            auto itIdx = pointCloudConnectivity.find(pt3);

            auto pointSetIt = pointCloud.begin();
            std::advance(pointSetIt, itIdx->second);
            if (is_building_pt[*pointSetIt])
                return false; //todo temp solution when having adjacent buildings

            auto it = flattenedPts.find(itIdx->second);
            if (it == flattenedPts.end()) {
                indices.push_back(itIdx->second);
                originalHeights.push_back(pointCloud.point(itIdx->second).z());
            }
        }
    }
    //-- Flatten surface points
    if (indices.empty()) {
        return false;
    }
    double avgHeight = geomutils::percentile(originalHeights,
                                             Config::get().flattenSurfaces[outputLayerID] / 100);

    //-- Add new points to the temp map
    for (auto& i : indices) {
        flattenedPts[i] = Point_3(pointCloud.point(i).x(), pointCloud.point(i).y(), avgHeight);
    }
    return true;
}

void PointCloud::buffer_flat_edges(const PolyFeaturesPtr& avgFeatures,
                                   std::vector<EPECK::Segment_3>& constrainedEdges) {
//    std::cout << "\n    Buffering flattening surfaces" << std::endl;
    //-- Add buffer around flattened polygons
    typedef CGAL::Straight_skeleton_builder_traits_2<EPICK> SsBuilderTraits;
    typedef CGAL::Straight_skeleton_2<EPICK> Ss;
    typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss> SsBuilder;
    typedef CGAL::Polygon_offset_builder_traits_2<EPICK> OffsetBuilderTraits;
    typedef CGAL::Polygon_offset_builder_2<Ss, OffsetBuilderTraits, Polygon_2> OffsetBuilder;

    typedef boost::shared_ptr<Polygon_2> ContourPtr;
    typedef std::vector<ContourPtr> ContourSequence;
    // get info using the original point cloud
//    std::vector<double> offsets{0.001, 0.2};
    std::vector<double> offsets{0.001};
    std::vector<Polygon_3> polyList;
    for (auto& f: avgFeatures) {
        auto& poly = f->get_poly().outer_boundary();
        // set the frame
        boost::optional<double> margin = CGAL::compute_outer_frame_margin(poly.begin(), poly.end(), offsets.back());
        CGAL::Bbox_2 bbox = CGAL::bbox_2(poly.begin(), poly.end());
        // Compute the boundaries of the frame
        double fxmin = bbox.xmin() - *margin;
        double fxmax = bbox.xmax() + *margin;
        double fymin = bbox.ymin() - *margin;
        double fymax = bbox.ymax() + *margin;
        // Create the rectangular frame
        Point_2 frame[4] = {Point_2(fxmin, fymin), Point_2(fxmax, fymin), Point_2(fxmax, fymax), Point_2(fxmin, fymax)
        };

        // Instantiate the skeleton builder
        SsBuilder ssb;
        // Enter the frame
        ssb.enter_contour(frame, frame + 4);
        // Enter the polygon as a hole of the frame (NOTE: as it is a hole we insert it in the opposite orientation)
        poly.reverse_orientation();
        ssb.enter_contour(poly.begin(), poly.end());
        // Construct the skeleton
        boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
        // Proceed only if the skeleton was correctly constructed.
        if (ss) {
            for (auto& offset: offsets) {
                // Instantiate the container of offsetSmall contours
                ContourSequence offset_contours;
                // Instantiate the offsetSmall builder with the skeleton
                OffsetBuilder ob(*ss);
                // Obtain the offsetSmall contours
                ob.construct_offset_contours(offset, std::back_inserter(offset_contours));
                // Locate the offsetSmall contour that corresponds to the frame
                // That must be the outmost offsetSmall contour, which in turn must be the one
                // with the largetst unsigned area.
                auto offsend = offset_contours.end();
                double lLargestArea = 0.0;
                for (auto i = offset_contours.begin(); i != offset_contours.end(); ++i) {
                    double lArea = CGAL_NTS abs((*i)->area()); //Take abs() as  Polygon_2::area() is signed.
                    if (lArea > lLargestArea) {
                        offsend = i;
                        lLargestArea = lArea;
                    }
                }
                // Remove the offsetSmall contour that corresponds to the frame.
                offset_contours.erase(offsend);

#ifndef NDEBUG
                if (offset_contours.size() != 1) {
                    std::cout << "DEBUG: SOMETHING WRONG WITH SKELETON, NUMBER OF POLYS: " << offset_contours.size()
                              << std::endl;
                }
#endif
                std::vector<double> height;
                Polygon_2& offsetPoly2 = *(offset_contours.back());
                geomutils::interpolate_poly_from_pc(offsetPoly2, height, _pointCloudTerrain);
                Polygon_3 offsetPoly3;
                for (auto j = 0; j < offsetPoly2.size(); ++j) {
                    offsetPoly3.push_back(ePoint_3(offsetPoly2[j].x(), offsetPoly2[j].y(), height[j]));
                }
                polyList.push_back(offsetPoly3);
            }
        }
    }
    //-- Add only non-adjacent edges to constrain
    CDT cdt;
    for (const auto& ring : polyList) {
        cdt.insert_constraint(ring.begin(), ring.end());
    }
    geomutils::mark_domains(cdt);

    for (CDT::Edge e : cdt.constrained_edges()) {
        CDT::Face_handle f1 = e.first;
        CDT::Face_handle f2 = f1->neighbor(e.second);
        if (f1->info().nesting_level == 0 || f2->info().nesting_level == 0) {
            ePoint_3 p1 = e.first->vertex((e.second + 1) % 3)->point();
            ePoint_3 p2 = e.first->vertex((e.second + 2) % 3)->point();
            constrainedEdges.emplace_back(p1, p2);
        }
    }
}

void PointCloud::read_point_clouds() {
    //- Read ground points
    if (!Config::get().ground_xyz.empty()) {
        std::cout << "Reading ground points" << std::endl;
        IO::read_point_cloud(Config::get().ground_xyz, _pointCloudTerrain);
        _pointCloudTerrain.add_property_map<bool>("is_building_point", false);

        std::cout << "    Points read: " << _pointCloudTerrain.size() << std::endl;
    } else {
        std::cout << "INFO: Did not find any ground points! Will calculate ground as a flat surface." << std::endl;
        std::cout << "WARNING: Ground elevation of buildings can only be approximated. "
                  << "If you are using point cloud to reconstruct buildings, building height estimation can be wrong.\n"
                  << std::endl;
    }

    //- Read building points
    if (!Config::get().buildings_xyz.empty()) {
        std::cout << "Reading building points" << std::endl;
        IO::read_point_cloud(Config::get().buildings_xyz, _pointCloudBuildings);
        if (_pointCloudBuildings.empty()) throw std::invalid_argument("Didn't find any building points!");

        std::cout << "    Points read: " << _pointCloudBuildings.size() << std::endl;
    }
}

Point_set_3& PointCloud::get_terrain() {
    return _pointCloudTerrain;
}
Point_set_3& PointCloud::get_buildings() {
    return _pointCloudBuildings;
}

const Point_set_3& PointCloud::get_terrain() const {
    return _pointCloudTerrain;
}

const Point_set_3& PointCloud::get_buildings() const {
    return _pointCloudBuildings;
}