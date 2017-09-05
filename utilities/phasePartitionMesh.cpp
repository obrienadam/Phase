#include <iostream>

#include <boost/filesystem.hpp>
#include <cgnslib.h>

#include "Input.h"
#include "ConstructGrid.h"
#include "CgnsViewer.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    auto comm = std::make_shared<Communicator>();

    if(comm->isMainProc())
    {
        int fn;
        cg_open(argv[1], CG_MODE_READ, &fn);

        char name[256];
        int cell_dim, phys_dim;
        cg_base_read(fn, 1, name, &cell_dim, &phys_dim);
        int nzones;
        cg_nzones(fn, 1, &nzones);

        for(int zone = 1; zone <= nzones; ++zone)
        {
            cgsize_t size[2];
            cg_zone_read(fn, 1, zone, name, size);
            cgsize_t rmin = 1, rmax = size[0];
            std::vector<double> coordsX(size[0]), coordsY(size[0]);
            cg_coord_read(fn, 1, zone, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, coordsX.data());
            cg_coord_read(fn, 1, zone, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, coordsY.data());
        }

        cg_close(fn);
    }

    Communicator::finalize();

    return 0;
}