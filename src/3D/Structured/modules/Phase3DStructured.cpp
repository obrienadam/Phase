#include "Structured/Solvers/SolverFactory.h"
#include "Structured/StructuredGrid3D/FaceStencil.h"
#include "Structured/PostProcessing/CgnsViewer.h"

#include "Geometry/Intersections.h"

int main(int argc, char *argv[])
{
    using namespace std;
    Communicator::init(argc, argv);

//

//    CommandLine cl;

//    cl.addSwitch("restart,r", "restart the solution from the latest time point");
//    cl.addSwitch("use-partitioned-grid,g", "use the pre-partitioned grid");

//    cl.parseArguments(argc, argv);

      Input input;
      input.parseInputFile();

      Sphere sp(Point3D(1., 2., 3.), 1.);
      Ray3D r(sp.centroid(), Point3D(254, 54, 135));

      r = Ray3D(r(0.00983), -r.r());

      for(auto pt: intersections(sp, r))
      {
          std::cout << pt << " " << (pt - r.x0()).mag() << std::endl;
      }


      Plane pn(Point3D(3., 0.431, 0.978), Point3D(1., 3., -1.), Point3D(1., 1., 0.2345));

      std::cout << pn.n() << "\n";
      std::cout << pn.nearestPoint(Point3D(4, 2, 4)) << "\n";

      std::cout << intersection(pn, LineSegment3D(Point3D(-1, 0, -2), Point3D(2, 2, -4))).first << "\n";

//    auto grid = std::make_shared<StructuredGrid3D>(input);

//    std::shared_ptr<Solver> solver = SolverFactory::create(input, grid);

//    //PostProcessing postProcessing(input, *solver, std::weak_ptr<const ImmersedBoundary>());

//    //RunControl runControl;

//    //runControl.run(cl, input, *solver, postProcessing);

//

//    auto grid = std::make_shared<StructuredGrid3D>(10, 10, 10, 10, 10, 10);
//    auto solver = SolverFactory().create(SolverFactory::POISSON, input, grid);

//    CgnsViewer v(input, solver);

//    std::cout << "--- Grid initialized ---\n";

//    FaceStencil st((*grid)(4, 4, 4), Face::J_POS, 3, 3);

//    std::cout << (*grid)(6, 9, 5, Face::J_POS).centroid() << std::endl;

//    for(auto c: st.coeffs())
//        std::cout << c << "\n";

//    solver->solve(0);


//    v.write(0);

//    auto st = FaceStencil(grid(0, 3, 4), Cell::I_NEG, 2, 2);

//    for(const Cell &c: st.cells())
//        std::cout << c.i() << "," << c.j() << "," << c.k() << "\n";

//    for(auto c: st.coeffs())
//        std::cout << c << "\n";

//    for(int k = 2; k < grid.nCellsK() - 2; ++k)
//        for(int j = 2; j < grid.nCellsJ() - 2; ++j)
//            for(int i = 2; i < grid.nCellsI() - 2; ++i)
//            {
//                    FaceStencil st1(grid(23, 43, 19), Cell::I_POS, 1, 1);
//                    FaceStencil st2(grid(23, 43, 19), Cell::I_NEG, 1, 1);
//                    FaceStencil st3(grid(23, 43, 19), Cell::J_POS, 1, 1);
//                    FaceStencil st4(grid(23, 43, 19), Cell::J_NEG, 1, 1);
//                    FaceStencil st5(grid(23, 43, 19), Cell::K_POS, 1, 1);
//                    FaceStencil st6(grid(23, 43, 19), Cell::K_NEG, 1, 1);
//            }

    Communicator::finalize();
}
