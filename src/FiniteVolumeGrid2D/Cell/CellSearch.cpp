#include <CGAL/Point_set_2.h>

#include "CellSearch.h"

std::vector< Ref<const Cell> > rangeSearch(const FiniteVolumeGrid2D &grid, const Circle& circle)
{
    typedef CGAL::Point_set_2<Kernel>::Vertex_handle VertexHandle;

    std::vector<Kernel::Point_2> centroids(grid.cells().size());
    std::map<Kernel::Point_2, const Cell* > directory;
    std::vector< Ref<const Cell> > result;

    CGAL::Point_set_2<Kernel> pointSet;

    //- Construct a point set consisting of the cell centroids
    for(const Cell &cell: grid.activeCells())
    {
        auto centroid = cell.centroid().cgalPoint();

        centroids[cell.id()] = centroid;
        directory.insert(std::pair<Kernel::Point_2, const Cell* >(centroid, &cell));
    }

    pointSet.insert(centroids.begin(), centroids.end());

    //- Find all of the centroids within the circle

    std::vector<VertexHandle> lv;
    Kernel::Circle_2 cgalCircle = circle.cgalCircle();
    pointSet.range_search(cgalCircle, std::back_inserter(lv));

    for(const VertexHandle &vtxh: lv)
        result.push_back(Ref<const Cell>(*directory[vtxh->point()]));

    return result;
}

std::vector< std::vector< Ref<const Cell> > > rangeSearch(const FiniteVolumeGrid2D& grid, Scalar radius)
{
    typedef CGAL::Point_set_2<Kernel>::Vertex_handle VertexHandle;

    std::vector<Kernel::Point_2> centroids(grid.cells().size());
    std::map<Kernel::Point_2, const Cell* > directory;
    std::vector< std::vector< Ref<const Cell> > > result;

    CGAL::Point_set_2<Kernel> pointSet;

    //- Construct a point set consisting of the cell centroids
    for(const Cell &cell: grid.activeCells())
    {
        auto centroid = cell.centroid().cgalPoint();

        centroids[cell.id()] = centroid;
        directory.insert(std::pair<Kernel::Point_2, const Cell* >(centroid, &cell));
    }

    pointSet.insert(centroids.begin(), centroids.end());
    result.resize(centroids.size());

    //- For each cell, find the neighbouring centroids, and then use the directory to construct the result

    std::vector<VertexHandle> lv;
    for(const Cell &cell: grid.activeCells())
    {
        Kernel::Circle_2 circle(cell.centroid().cgalPoint(), radius*radius);
        pointSet.range_search(circle, std::back_inserter(lv));

        for(const VertexHandle &vtxh: lv)
            result[cell.id()].push_back(Ref<const Cell>(*directory[vtxh->point()]));

        lv.clear();
    }

    return result;
}

std::vector< Ref<const Cell> > kNearestNeighbourSearch(const FiniteVolumeGrid2D &grid, const Cell& cell, size_t k)
{
    typedef CGAL::Point_set_2<Kernel>::Vertex_handle VertexHandle;

    std::vector<Kernel::Point_2> centroids(grid.cells().size());
    std::map<Kernel::Point_2, const Cell* > directory;
    std::vector< Ref<const Cell> > result;

    CGAL::Point_set_2<Kernel> pointSet;

    //- Construct a point set consisting of the cell centroids
    for(const Cell &cell: grid.activeCells())
    {
        auto centroid = cell.centroid().cgalPoint();
        centroids[cell.id()] = centroid;
        directory.insert(std::pair<Kernel::Point_2, const Cell* >(centroid, &cell));
    }

    pointSet.insert(centroids.begin(), centroids.end());
    result.reserve(k);

    //- For each cell, find the neighbouring centroids, and then use the directory to construct the result

    std::vector<VertexHandle> lv;

    lv.reserve(k);
    pointSet.nearest_neighbors(cell.centroid().cgalPoint(), k, std::back_inserter(lv));

    for(const VertexHandle &vtxh: lv)
        result.push_back(Ref<const Cell>(*directory[vtxh->point()]));

    return result;
}

std::vector< std::vector< Ref<const Cell> > > kNearestNeighbourSearch(const FiniteVolumeGrid2D& grid, size_t k)
{
    typedef CGAL::Point_set_2<Kernel>::Vertex_handle VertexHandle;

    std::vector<Kernel::Point_2> centroids(grid.cells().size());
    std::map<Kernel::Point_2, const Cell* > directory;
    std::vector< std::vector< Ref<const Cell> > > result;

    CGAL::Point_set_2<Kernel> pointSet;

    //- Construct a point set consisting of the cell centroids
    for(const Cell &cell: grid.activeCells())
    {
        auto centroid = cell.centroid().cgalPoint();

        centroids[cell.id()] = centroid;
        directory.insert(std::pair<Kernel::Point_2, const Cell* >(centroid, &cell));
    }

    pointSet.insert(centroids.begin(), centroids.end());
    result.resize(centroids.size());

    //- For each cell, find the neighbouring centroids, and then use the directory to construct the result

    std::vector<VertexHandle> lv;
    for(const Cell &cell: grid.activeCells())
    {
        pointSet.nearest_neighbors(cell.centroid().cgalPoint(), k, std::back_inserter(lv));

        for(const VertexHandle &vtxh: lv)
            result[cell.id()].push_back(Ref<const Cell>(*directory[vtxh->point()]));

        lv.clear();
    }

    return result;
}
