#define BOOST_TEST_MODULE RangeSearchTest

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "RectilinearGrid2D.h"
#include "CellSearch.h"

BOOST_AUTO_TEST_SUITE (RangeSearchTest)

BOOST_AUTO_TEST_CASE(RangeSearchTest)
{
    RectilinearGrid2D grid(10, 10, 0.1, 0.1);

    auto result = rangeSearch(grid, 0.3);

    size_t i = 0;
    for(const auto &nbList: result)
    {
        std::cout << "\nFor cell " << i++ << " we have: ";

        for(const Cell &cell: nbList)
        {
            std::cout << cell.id() << ", ";
        }
    }
}

BOOST_AUTO_TEST_CASE(KNearestNeighbourSearchTest)
{
    RectilinearGrid2D grid(5, 5, 0.1, 0.1);

    auto result = kNearestNeighbourSearch(grid, 5);

    size_t i = 0;
    for(const auto &nbList: result)
    {
        std::cout << "\nFor cell " << i++ << " we have: ";

        for(const Cell &cell: nbList)
        {
            std::cout << cell.id() << ", ";
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
