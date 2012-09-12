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
const int INIT_REF_NUM = 6;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;

const int N = 24;

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;
// Only for matrix_solver_type == SOLVER_AZTECOO :
const char* iterative_method = "bicgstab";
const char* preconditioner = "jacobi";

int main(int argc, char* args[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.setParamValue(Hermes::exceptionsPrintCallstack, 1);
  Hermes::Hermes2D::Hermes2DApi.setParamValue(Hermes::Hermes2D::numThreads, 2);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  Hermes::vector<Mesh *> meshes;
  for (int i = 0; i < N; i++)
    meshes.push_back(new Mesh());
  
  MeshReaderH2D mloader;
  mloader.load("square2.mesh", meshes[0]);
  
  for (int i = 1; i < N; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM; j++) 
      meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM; j++) 
    meshes[0]->refine_all_elements();

  Hermes::vector<const Space<double> *> spaces;
  
  for (int n = 0; n < N; n++)
    spaces.push_back(new L2Space<double>(meshes[n], P_INIT));
  
  int ndof =  Space<double>::get_num_dofs(spaces);

  // Display the mesh.
  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  oview.show(spaces[0]);
  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
  bview.show(spaces[0]);

  Hermes::vector<Solution<double>* > slns;
  for (int i = 0; i < N; i++)
    slns.push_back(new Solution<double>());

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_in", meshes[0], N);
 
  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, spaces);
  dp.set_fvm();
  LinearSolver<double> solver(&dp);
  
  Loggable::Static::info("Solving. NDOF = %d", ndof);
  cpu_time.tick();

    // Solve the linear system. If successful, obtain the solution.
    try
    {
      solver.solve();
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, slns);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
  
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());
  
  // Wait for keyboard or mouse input.
  // View the coarse mesh solution.
  for (int n = 0; n < N; n++)
  {
    if(HERMES_VISUALIZATION)
    {
      ScalarView view("Solution", new WinGeom(900, 0, 450, 350));
      view.fix_scale_width(60);
      view.show(slns[n]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(slns[n], (std::string("sln_") + itos(n) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  return 0;
}
