#define BOOST_TEST_MODULE PolygonTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "Polygon.h"

BOOST_GEOMETRY_REGISTER_POINT_2D(Point2D, Scalar, cs::cartesian, x, y)

BOOST_AUTO_TEST_SUITE (PolygonTest)

BOOST_AUTO_TEST_CASE(IntersectionTest)
{
    typedef boost::geometry::model::polygon<Point2D, false> Polygon2D;

    Polygon2D pgn1, pgn2;

    boost::geometry::append(pgn1, Point2D(0, 0));
    boost::geometry::append(pgn1, Point2D(1.2, 0.3));
    boost::geometry::append(pgn1, Point2D(1.4, 0.91));
    boost::geometry::append(pgn1, Point2D(0, 0.71));
    boost::geometry::append(pgn1, Point2D(0, 0));

    boost::geometry::append(pgn2, Point2D(0, 0));
    boost::geometry::append(pgn2, Point2D(0.42, 0.3));
    boost::geometry::append(pgn2, Point2D(.34, 0.91));
    boost::geometry::append(pgn2, Point2D(0, 0.71));
    boost::geometry::append(pgn2, Point2D(0, 0));

    Point2D c1, c2;

    boost::geometry::centroid(pgn1, c1);
    boost::geometry::centroid(pgn2, c2);

    std::cout << boost::geometry::area(pgn1) << std::endl;
    std::cout << boost::geometry::area(pgn2) << std::endl;
    std::cout << c1 << std::endl;
    std::cout << c2 << std::endl;

    std::vector<Polygon2D> output;

    boost::geometry::intersection(pgn1, pgn2, output);

    std::cout << output.size() << std::endl;

    for(const Polygon2D& pgn: output)
        std::cout << boost::geometry::area(pgn) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
