#define BOOST_TEST_MODULE PolygonTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "Polygon.h"
#include "Line2D.h"

BOOST_AUTO_TEST_SUITE (PolygonTest)

BOOST_AUTO_TEST_CASE(NearestIntersectionTest)
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

BOOST_AUTO_TEST_CASE(PolygonIntersectionTest)
{
    std::vector<Point2D> verts1 = {
      Point2D(0, 0),
        Point2D(0.452, -0.1),
        Point2D(0.9, 0.1),
        Point2D(0.75, 0.9),
        Point2D(0.2, 0.6),
        Point2D(-0.1, 0.3),
    };

    std::vector<Point2D> verts2 = {
        Point2D(0, 0),
          Point2D(0.2, -0.3),
          Point2D(0.5, 0.4),
          Point2D(1.2, 0.5),
          Point2D(0.5, 1.2),
          Point2D(0.1, 1),
    };

    Polygon pgnA(verts1), pgnB(verts2);

    std::cout << "Area and centroid of polygon A: " << pgnA.centroid() << " " << pgnA.area() << "\n"
              << "Area and centroid of polygon B: " << pgnB.centroid() << " " << pgnB.area() << "\n";

    Polygon pgnC = intersectionPolygon(pgnA, pgnB);

    std::cout << "Area and centroid of polygon C: " << pgnC.centroid() << " " << pgnC.area() << "\n\n";

    for(const auto &vtx: pgnC)
        std::cout << vtx << "\n";
}

BOOST_AUTO_TEST_CASE(PolygonSplitTest)
{
    std::vector<Point2D> verts1 = {
      Point2D(0, 0),
        Point2D(0.452, -0.1),
        Point2D(0.9, 0.1),
        Point2D(0.75, 0.9),
        Point2D(0.2, 0.6),
        Point2D(-0.1, 0.3),
    };

    Polygon pgn(verts1);

    Line2D lA(Point2D(1, 1), Point2D(-1.342, 2)),
            lB(Point2D(0.8, 0.4), Point2D(1.2, 0.7));

    std::cout << "Line A: " << lA << "\n"
              << "Line B: " << lB << "\n"
              << "Intersection: " << Line2D::intersection(lA, lB).first << "\n\n";

    Point2D norm(-0.32, -2), o = pgn.centroid() + norm.unitVec()*0.2;;
    Line2D split(o, norm);

    Polygon clippedPgn = clipPolygon(pgn, split);

    for(const auto &vtx: pgn)
        std::cout << vtx.x << " " << vtx.y << "\n";
    std::cout << "\n";
    for(const auto &vtx: clippedPgn)
        std::cout << vtx.x << " " << vtx.y << "\n";

    std::cout << "Clipped polygon vol frac: " << clippedPgn.area()/pgn.area() << "\n";

    BOOST_REQUIRE_CLOSE(clipPolygon(pgn, Line2D(o, norm.normalVec())).area()/pgn.area(),
                        1. - clipPolygon(pgn, Line2D(o, (-norm).normalVec())).area()/pgn.area(),
                        1e-8);
}


BOOST_AUTO_TEST_SUITE_END()
