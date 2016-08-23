#define BOOST_TEST_MODULE PolygonTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "Polygon.h"
#include "Line2D.h"
#include "Circle.h"

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

    for(const Point2D &vtx: pgn.vertices())
    {
        std::cout << vtx << std::endl;
    }

    pgn -= Point2D(1, 1);

    std::cout << std::endl
              << "Nearerest intersect for " << pt2 << " = " << pgn.nearestIntersect(pt2) << std::endl;
}

BOOST_AUTO_TEST_CASE(FirstIntersectionTest)
{
    Point2D ptA(3, 2), ptB(1, -4);
    Circle circ(Point2D(1, 1), 2.);

    auto intersections = circ.intersections(Line2D(ptA, (ptB - ptA).tangentVec()));

    if(intersections.size() != 0)
    {
        std::cout << "At least one intersection was found. Pt = " << nearestPoint(ptA, intersections) << "\n";
    }
    else
    {
        std::cout << "No intersections were found.\n";
    }
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

    for(const auto &vtx: pgnC.vertices())
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

    for(const auto &vtx: pgn.vertices())
        std::cout << vtx.x << " " << vtx.y << "\n";
    std::cout << "\n";
    for(const auto &vtx: clippedPgn.vertices())
        std::cout << vtx.x << " " << vtx.y << "\n";

    std::cout << "Clipped polygon vol frac: " << clippedPgn.area()/pgn.area() << "\n";

    BOOST_REQUIRE_CLOSE(clipPolygon(pgn, Line2D(o, norm.normalVec())).area()/pgn.area(),
                        1. - clipPolygon(pgn, Line2D(o, (-norm).normalVec())).area()/pgn.area(),
                        1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
