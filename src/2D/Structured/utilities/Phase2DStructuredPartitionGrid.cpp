#include <iostream>
#include <numeric>

#include <boost/program_options.hpp>

#include "System/CommandLine.h"
#include "System/CgnsFile.h"

int main(int argc, char **argv)
{
    using namespace std;
    namespace po = boost::program_options;

    CommandLine cl;

    cl.addOptions()
              ("num-partitions,n", po::value<int>()->required(), "number of partitions to generate");

    cl.addSwitch("strip,s", "perform a strip decomposition");

    cl.parseArguments(argc, argv);

    int nparts = cl.get<int>("num-partitions");

    //- Figure out the partitioning
    std::vector<int> factors;

    if (cl.get<bool>("strip"))
        factors.assign({nparts, 1});
    else
    {
        int i = 2, num = nparts;

        while (i * i <= num)
            if (num % i)
                ++i;
            else
            {
                factors.push_back(i);
                num /= i;
            }

        factors.push_back(num);

        while (factors.size() > 2)
        {
            factors.push_back(factors[0] * factors[1]);
            factors.erase(factors.begin(), factors.begin() + 2);
            std::sort(factors.begin(), factors.end());
        }

        if (factors.size() < 2)
            factors.push_back(1);
    }


    int nnodesI = 2001, nnodesJ = 1001;
    int ncellsI = nnodesI - 1, ncellsJ = nnodesJ - 1;


    for (int blockJ = 0; blockJ < factors[1]; ++blockJ)
        for (int blockI = 0; blockI < factors[0]; ++blockI)
        {
            int partNo = blockJ * factors[0] + blockI;

            int ibegin = blockI * nnodesI / factors[0];
            int iend = min((blockI + 1) * nnodesI / factors[0] + 1, nnodesI);

            int jbegin = blockJ * nnodesJ / factors[1];
            int jend = min((blockJ + 1) * nnodesJ / factors[1] + 1, nnodesJ);

            int nPartNodesI = iend - ibegin;
            int nPartNodesJ = jend - jbegin;
            int nPartCellsI = nPartNodesI - 1;
            int nPartCellsJ = nPartNodesJ - 1;

            ostringstream sout;
            sout << "part_" << blockI << "_" << blockJ;

            CgnsFile file("test.cgns", partNo == 0 ? CgnsFile::WRITE : CgnsFile::MODIFY);

            int bid = partNo == 0 ? file.createBase("Grid", 2, 2) : 1;

            int zid = file.createStructuredZone(bid, sout.str(), nPartNodesI, nPartNodesJ, nPartCellsI, nPartCellsJ);

            std::vector<Point2D> nodes;

            for (int j = jbegin; j < jend; ++j)
                for (int i = ibegin; i < iend; ++i)
                    nodes.push_back(Point2D(i, j));

            file.writeCoordinates(bid, zid, nodes);

            int sid = file.writeSolution(bid, zid, "info");

            int fid = file.writeField(bid, zid, sid, "proc", std::vector<int>(nPartCellsI * nPartCellsJ, partNo));

            file.close();
        }

    return 0;
}

