#include <CGAL/Point_set_2.h>

#include "RangeSearch.h"

RangeSearch::RangeSearch(const FiniteVolumeGrid2D &grid, Scalar radius)
{
    search(grid, radius);
}

void RangeSearch::search(const FiniteVolumeGrid2D &grid, Scalar radius)
{
    typedef CGAL::Point_set_2<Kernel>::Vertex_handle VertexHandle;

    std::vector<Kernel::Point_2> centroids(grid.cells().size());
    std::map<Kernel::Point_2, const Cell* > directory;

    CGAL::Point_set_2<Kernel> pointSet;

    //- Construct a point set consisting of the cell centroids
    for(const Cell &cell: grid.activeCells())
    {
        auto centroid = cell.centroid().cgalPoint();

        centroids[cell.id()] = centroid;
        directory.insert(std::pair<Kernel::Point_2, const Cell* >(centroid, &cell));
    }

    pointSet.insert(centroids.begin(), centroids.end());
    result_.resize(centroids.size());

    //- For each cell, find the neighbouring centroids, and then use the directory to construct the result
    for(const Cell &cell: grid.activeCells())
    {
        std::vector<VertexHandle> lv;
        Kernel::Circle_2 circle(cell.centroid().cgalPoint(), radius*radius);
        pointSet.range_search(circle, std::back_inserter(lv));

        for(const VertexHandle &vtxh: lv)
            result_[cell.id()].push_back(Ref<const Cell>(*directory[vtxh->point()]));
    }
}
