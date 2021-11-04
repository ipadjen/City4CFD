#ifndef CITYCFD_GEOMUTILS_H
#define CITYCFD_GEOMUTILS_H

#include "config.h"

namespace geomutils {
    double  avg(const std::vector<double>& values);
    double  percentile(std::vector<double> values, const double percentile);
    bool    point_in_circle(const Point_3& pt, const Point_2& center, const double& radius);
    void    cdt_to_mesh(CDT& cdt, Mesh& mesh, const int surfaceLayerID = -9999);
    void    mark_domains(CDT& cdt, PolyFeatures features = {});
    void    mark_domains(CDT& ct, const Face_handle& start, int index,
                         std::list<CDT::Edge>& border, PolyFeatures& features);
    void    shorten_long_poly_edges(Polygon_2& poly, double maxLen = config::edgeMaxLen);
    Point_2 rotate_pt(Point_2& pt, const double angle, Point_2 centerPt = Point_2(0, 0));
    void    interpolate_poly_from_pc(const Polygon_2& poly, std::vector<double>& heights, const Point_set_3& pointCloud);

    //-- Templated functions
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_with_holes_2& polygon);
    template <typename T> bool point_in_poly(const T& pt2, const Polygon_2& polygon);
    template <typename T> void make_round_poly(Point_2& centre, double radius, T& poly);
    template <typename T> void make_round_poly(Point_2& centre, double radius1, double radius2,
                                               int nPts, double angInt, double ang, T& poly);
    template <typename T, typename U> void smooth_dt (const Point_set_3& pointCloud, T& dt);
    template <typename T> Polygon_2 calc_bbox_poly(const T& inputPts);
}

#endif //CITYCFD_GEOMUTILS_H