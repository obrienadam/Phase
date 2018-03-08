#include "FiniteVolumeGrid2DFactory.h"
#include "CgnsUnstructuredGrid.h"
#include "StructuredRectilinearGrid.h"

std::shared_ptr<FiniteVolumeGrid2D> FiniteVolumeGrid2DFactory::create(GridType type, const Input &input)
{
    std::shared_ptr<FiniteVolumeGrid2D> grid;

    switch (type)
    {
        case CGNS:
            grid = std::make_shared<CgnsUnstructuredGrid>(input);
            break;
        case RECTILINEAR:
            grid = std::make_shared<StructuredRectilinearGrid>(input);
            break;
        case RELOAD:
            grid = std::make_shared<CgnsUnstructuredGrid>();
            std::static_pointer_cast<CgnsUnstructuredGrid>(grid)->loadPartitionedGrid(std::make_shared<Communicator>());
            return grid;
        default:
            return nullptr;
    }

    grid->partition(input, std::make_shared<Communicator>());
    return grid;
}

std::shared_ptr<FiniteVolumeGrid2D> FiniteVolumeGrid2DFactory::create(std::string type, const Input &input)
{
    std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c)
    {
        return std::tolower(c);
    });

    if (type == "cgns")
        return create(CGNS, input);
    else if (type == "rectilinear")
        return create(RECTILINEAR, input);
    else if (type == "reload")
        return create(RELOAD, input);

    return nullptr;
}

std::shared_ptr<FiniteVolumeGrid2D> FiniteVolumeGrid2DFactory::create(const Input &input)
{
    return create(input.caseInput().get<std::string>("Grid.type"), input);
}