
<!-- Page Content -->
<div class="container">

    <!-- Page Heading -->
    <div class="row">
        <div class="col-lg-12">
            <h1 class="page-header">Phase
                <small>Unstructured 2D Multiphase Solver</small>
            </h1>
        </div>
    </div>
    <!-- /.row -->

    <div class="row">
        <div class="col-md-12">
            <h3>About</h3>
            <p>
                Phase is a lightweight unstructured finite-volume solver for multiphase fluid equations. The intention of this
                project is to provide a platform allowing for rapid development and testing of new finite-volume methods. The
                framework allows the programmer to implement finite-volume equations almost as they would be written on paper.
                For example, if one wishes to implement the equations for the SIMPLE method of Patankar and Spalding, the following
                code can be used for the momentum equation:
            </p>

            <pre class="prettyprint lang-cpp linenums">
uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) == fv::laplacian(mu, u) - fv::grad(p));
uEqn_.relax(omegaMomentum_);
double error = uEqn_.solve();</pre>

            <p>
                And, for the pressure correction equation:
            </p>

            <pre class="prettyprint lang-cpp linenums">
pCorrEqn_ = (fv::laplacian(rho*d, pCorr) == m);
double error = pCorrEqn_.solve();</pre>

            <p>
                The <code class="prettyprint lang-cpp">fv</code> namespace contains all of the standard finite-volume
                operators. Users can easily define
                their own discretization schemes for either advection or diffusion and combine them with the operators that
                already exist within the framework.
            </p>

            <p>
                Each of the finite-volume functions constructs a sparse matrix representing the discretized form of the dependent variable
                and combines them. The euqation itself is solved using the efficient multi-threaded Bi-Conjugate Gradient Stabilized (BiCGStab)
                iterative method in conjuction with an Incomplete Lower-Upper Triangular (ILUT) factorization preconditioner available in the
                freely available <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a> C++ library.
            </p>

            <h3>Installation</h3>

            <p>
                Phase depends on the following third-party development libaries:
            <ul>
                <li>
                    CMake 2.8.11
                </li>
                <li>
                    Boost 1.53
                </li>
                <li>
                    Eigen 3.2
                </li>
                <li>
                    CGAL 4.7
                </li>
            </ul>
            </p>

            <p>
                The most recent version can be found <a href="https://www.github.com/obrienadam/Phase">here</a>. To clone
                this repository, open a terminal and type the commmand:
            </p>

            <p>
                <samp>> git clone https://www.github.com/obrienadam/Phase</samp>
            </p>

            <p>
                Then, assuming that the repository was cloned into the home directory, navigate to the project directory:
            </p>

            <p>
                <samp>> cd ~/Phase</samp>
            </p>

            <p>
                and run the command
            </p>

            <p>
                <samp>> cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release</samp>
            </p>

            <p>
                to generate the make files. Then, build the project with
            </p>

            <p>
                <samp>> make all</samp>
            </p>

            <p>
                The solver binaries will be placed into <samp>$HOME/Phase/modules</samp>. In order to run a case, such
                as the 2D lid driven cavity case, navigate to the case directory:
            </p>

            <p>
                <samp>> cd ~/Phase/Examples/LidDrivenCavity</samp>
            </p>

            <p>
                and execute the desired solver from the case directory; for example:
            </p>

            <p>
                <samp>> ~/Phase/modules/phasePiso</samp>
            </p>
        </div>
    <hr>

    <!-- Footer -->
    <footer>
        <div class="row">
            <div class="col-lg-12">
                <p>Copyright &copy; Adam O'Brien 2015</p>
            </div>
        </div>
        <!-- /.row -->
    </footer>

</div>
<!-- /.container -->
