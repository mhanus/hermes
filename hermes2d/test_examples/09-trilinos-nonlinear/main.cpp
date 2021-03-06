#include "definitions.h" 

using namespace Teuchos;

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It
//  compares performance of the Newton's method in Hermes (assembling via the DiscreteProblem
//  class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX
//  solver (using the Hermes DiscreteProblem class to assemble discrete problems).
//
//  PDE:  - \nabla (k \nabla u) - f = 0
//  k = (1 + sqr(u_x) + sqr(u_y))^{-0.5}
//
//  Domain: Unit square.
//
//  BC: zero Dirichlet.
//
//  Exact solution: (x - x*x) * (y - y*y).
//
//  Initial guess for the Newton's method: zero function.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.
const int P_INIT = 4;                             // Initial polynomial degree of all mesh elements.

const bool JFNK = true;                          // true = jacobian-free method,
// false = Newton.
const int PRECOND = 1;                            // Preconditioning by jacobian (1) or approximation of jacobian (2)
// in case of JFNK,
// Default ML proconditioner in case of Newton.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
// by the other solvers).
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
// the other solvers).
// Possibilities: none, jacobi, neumann, least-squares, or a
//  preconditioner from IFPACK (see solver/aztecoo.h)
// NOX parameters.
unsigned message_type = 0;//NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
// NOX error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-8;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-8;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space1(new H1Space<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> space2(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = Space<double>::get_num_dofs(space1);

  // Initialize weak formulation,
  CustomWeakForm wf1;

  // Initialize the discrete problem.
  DiscreteProblem<double> dp1(&wf1, space1);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  MeshFunctionSharedPtr<double> sln1(new Hermes::Hermes2D::Solution<double>());
  MeshFunctionSharedPtr<double> sln2(new Hermes::Hermes2D::Solution<double>());

  Hermes::Hermes2D::NewtonSolver<double> newton(&dp1);
  newton.set_verbose_output(true);
  try
  {
    newton.solve();
  }
  catch (Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }

  Solution<double>::vector_to_solution(newton.get_sln_vector(), space1, sln1);

  // Show UMFPACK solution.
  Views::ScalarView view1("Solution 1", new Views::WinGeom(0, 0, 500, 400));
  view1.show(sln1);
  Views::View::wait();

  MeshFunctionSharedPtr<double> ex(new CustomExactSolution(mesh));

  // TRILINOS PART:

  // Initialize the weak formulation for Trilinos.
  CustomWeakForm wf2(JFNK, PRECOND == 1, PRECOND == 2);

  // Initialize DiscreteProblem.
  DiscreteProblemNOX<double> dp2(&wf2, space2);

  // Initialize the NOX solver with the vector "coeff_vec".
  NewtonSolverNOX<double> nox_solver(&dp2);
  nox_solver.set_output_flags(message_type);
  nox_solver.set_ls_tolerance(ls_tolerance);
  nox_solver.set_conv_rel_resid(rel_resid);
  nox_solver.set_conv_iters(max_iters);

  // Choose preconditioning.
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the nonlinear problem using NOX.
  try
  {
    nox_solver.solve(NULL);
  }
  catch (Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }
  Solution<double>::vector_to_solution(nox_solver.get_sln_vector(), space2, sln2);

  // Calculate error.
  errorCalculator.calculate_errors(sln1, ex);
  double rel_err_1 = errorCalculator.get_error_squared(0) * 100;
  errorCalculator.calculate_errors(sln2, ex);
  double rel_err_2 = errorCalculator.get_error_squared(0) * 100;

  // Show NOX solution.
  Views::ScalarView view2("Solution 2", new Views::WinGeom(510, 0, 500, 400));
  view2.show(sln2);

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}