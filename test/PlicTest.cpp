#define BOOST_TEST_MODULE PlicTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "Multiphase.h"
#include "Plic.h"
#include "RectilinearGrid2D.h"
#include "Polygon.h"

BOOST_AUTO_TEST_SUITE (PlicTest)

BOOST_AUTO_TEST_CASE(PlicConstruction)
{
    RectilinearGrid2D grid(10, 10, 0.02, 0.01);

    Scalar gamma;
    Vector2D normal;
    Polygon phasePgn;
    for(const Cell &cell: grid.cells())
    {
        gamma = (Scalar)rand()/RAND_MAX;
        normal = Vector2D(-RAND_MAX/3 + rand(), -RAND_MAX/3 + rand());

        phasePgn = plic::computeInterfacePolygon(cell, gamma, normal);

        BOOST_REQUIRE_CLOSE(gamma, phasePgn.area()/cell.volume(), 1e-10); //- Seems to be the upper limit for accuracy
    }
}

BOOST_AUTO_TEST_CASE(PlicAdvection)
{

}

BOOST_AUTO_TEST_SUITE_END()

