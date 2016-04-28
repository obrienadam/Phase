#define BOOST_TEST_MODULE ImmersedBoundaryTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "ImmersedBoundary.h"
#include "RectilinearGrid2D.h"

BOOST_AUTO_TEST_SUITE (ImmersedBoundaryTest)

BOOST_AUTO_TEST_CASE(CellGroupCreationTest)
{
    RectilinearGrid2D grid(20, 20, 0.05, 0.05);
    ImmersedBoundaryObject ibObj(grid, Point2D(0.5, 0.5), 0.15);

    ibObj.constructStencils();

    for(const auto &ibCell: ibObj.stencils())
    {
        std::cout << "Stencil for cell " << ibCell.cell().id() << ": " << ibCell.boundaryPoint() << " " << ibCell.imagePoint() << "\n";
    }

    std::cout << "Active cell group size   : " << grid.activeCells().size() << "\n"
              << "Inactive cell group size : " << grid.inactiveCells().size() << "\n"
              << "Fluid cell group size    : " << grid.fluidCells().size() << "\n"
              << "IB cell group size       : " << grid.cellGroup("ibCells").size() << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
