#define BOOST_TEST_MODULE PolygonTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "Polygon.h"

BOOST_AUTO_TEST_SUITE (PolygonTest)

BOOST_AUTO_TEST_CASE(IntersectionTest)
{
    std::vector<Point2D> vertices = {
      Point2D(0, 0),
        Point2D(1, 0),
        Point2D(1, 1),
        Point2D(0, 1),
    };

    Polygon pgn(vertices);

    std::cout << pgn.area() << std::endl
              << pgn.centroid() << std::endl;

    Point2D pt1(-0.1, 0.5), pt2(0.1, 0.5);

    BOOST_ASSERT(!pgn.isInside(pt1));
    BOOST_ASSERT(pgn.isInside(pt2));

    pgn += Point2D(1, 1);

    std::cout << pgn.area() << std::endl
              << pgn.centroid() << std::endl;

    for(const Point2D &vtx: pgn)
    {
        std::cout << vtx << std::endl;
    }

    pgn -= Point2D(1, 1);

    std::cout << std::endl
              << "Nearerest intersect for " << pt2 << " = " << pgn.nearestIntersect(pt2) << std::endl;
}

BOOST_AUTO_TEST_CASE(AirfoilTest)
{
    // NACA-0009 9.0% smoothed
    std::vector<Point2D> vertices =
    {
         Point2D(1.00000,  0.0),
         Point2D(0.99572,  0.00057),
         Point2D(0.98296,  0.00218),
         Point2D(0.96194,  0.00463),
         Point2D(0.93301,  0.00770),
         Point2D(0.89668,  0.01127),
         Point2D(0.85355,  0.01522),
         Point2D(0.80438,  0.01945),
         Point2D(0.75000,  0.02384),
         Point2D(0.69134,  0.02823),
         Point2D(0.62941,  0.03247),
         Point2D(0.56526,  0.03638),
         Point2D(0.50000,  0.03978),
         Point2D(0.43474,  0.04248),
         Point2D(0.37059,  0.04431),
         Point2D(0.33928,  0.04484),
         Point2D(0.30866,  0.04509),
         Point2D(0.27886,  0.04504),
         Point2D(0.25000,  0.04466),
         Point2D(0.22221,  0.04397),
         Point2D(0.19562,  0.04295),
         Point2D(0.17033,  0.04161),
         Point2D(0.14645,  0.03994),
         Point2D(0.12408,  0.03795),
         Point2D(0.10332,  0.03564),
         Point2D(0.08427,  0.03305),
         Point2D(0.06699,  0.03023),
         Point2D(0.05156,  0.02720),
         Point2D(0.03806,  0.02395),
         Point2D(0.02653,  0.02039),
         Point2D(0.01704,  0.01646),
         Point2D(0.00961,  0.01214),
         Point2D(0.00428,  0.00767),
         Point2D(0.00107,  0.00349),
         Point2D(0.0    ,  0.0),
         Point2D(0.00107, -0.00349),
         Point2D(0.00428, -0.00767),
         Point2D(0.00961, -0.01214),
         Point2D(0.01704 ,-0.01646),
         Point2D(0.02653, -0.02039),
         Point2D(0.03806, -0.02395),
         Point2D(0.05156, -0.02720),
         Point2D(0.06699, -0.03023),
         Point2D(0.08427, -0.03305),
         Point2D(0.10332, -0.03564),
         Point2D(0.12408, -0.03795),
         Point2D(0.14645, -0.03994),
         Point2D(0.17033, -0.04161),
         Point2D(0.19562, -0.04295),
         Point2D(0.22221, -0.04397),
         Point2D(0.25000, -0.04466),
         Point2D(0.27886, -0.04504),
         Point2D(0.30866, -0.04509),
         Point2D(0.33928, -0.04484),
         Point2D(0.37059, -0.04431),
         Point2D(0.43474, -0.04248),
         Point2D(0.50000, -0.03978),
         Point2D(0.56526, -0.03638),
         Point2D(0.62941, -0.03247),
         Point2D(0.69134, -0.02823),
         Point2D(0.75000, -0.02384),
         Point2D(0.80438 ,-0.01945),
         Point2D(0.85355, -0.01522),
         Point2D(0.89668, -0.01127),
         Point2D(0.93301, -0.00770),
         Point2D(0.96194 ,-0.00463),
         Point2D(0.98296, -0.00218),
         Point2D(0.99572, -0.00057),
    };

    Polygon airfoil(vertices);
    airfoil.scale(10);
    airfoil.rotate(-6*M_PI/180.); // six degree AoA

    std::cout << "Airfoil area     : " << airfoil.area() << std::endl
              << "Airfoil centroid : " << airfoil.centroid() << std::endl;

    Point2D pt(0.3, 0.02);

    std::cout << "Test point        : " << pt << std::endl
              << "Inside airfoil?   : " << airfoil.isInside(pt) << std::endl
              << "Nearest intersect : " << airfoil.nearestIntersect(pt) << std::endl
              << "Image point       : " << 2*airfoil.nearestIntersect(pt) - pt << "\n\n";

    for(const Point2D &vtx: airfoil)
        std::cout << vtx.x << " " << vtx.y << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
