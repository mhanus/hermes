#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//  It is intended to show how evalutation of surface matrix forms that take basis functions defined
//  on different elements work. It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:    Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on[0,0.5] x {0}, and g = 0 anywhere else.
//
//  The following parameters can be changed:


// Set the following flag to "false" to suppress Hermes OpenGL visualization.
const bool HERMES_VISUALIZATION = true;
// Set the following flag to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Number of initial uniform mesh refinements.
const int INIT_REF = 2;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;
// Only for matrix_solver_type == SOLVER_AZTECOO :
const char* iterative_method = "bicgstab";
const char* preconditioner = "jacobi";

int main(int argc, char* args[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();

  // Create an L2 space.
  L2Space<double> space(&mesh, P_INIT);
  int ndof = space.get_num_dofs();

  // Display the mesh.
  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  oview.show(&space);
  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
  bview.show(&space);

  Solution<double> sln;

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_bottom_left", &mesh);

  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);
 
  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
  Vector<double>* rhs = create_vector<double>(matrix_solver_type);
  Solvers::LinearSolver<double>* solver = Solvers::create_linear_solver<double>(matrix_solver_type, matrix, rhs);

  info("Assembling (ndof: %d).", ndof);
  cpu_time.tick();
    dp.assemble(matrix, rhs);
  cpu_time.tick();
  info("Time taken: %lf s", cpu_time.last());
  
  // Solve the linear system. If successful, obtain the solution.
  if(solver->solve())
    Solution<double>::vector_to_solution(solver->get_sln_vector(), &space, &sln);
  else
    throw Hermes::Exceptions::Exception("Matrix solver failed.\n");

  // View the coarse mesh solution.
  view1.show(&sln);

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  // Wait for keyboard or mouse input.
  Views::View::wait();
  return 0;
}