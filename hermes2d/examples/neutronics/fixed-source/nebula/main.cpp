#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"
#include "weakforms_neutronics.h"

#include <iterator>
#include <iomanip>
const bool SAVE_SLN_VECTOR = true;

const bool HERMES_ANG_VISUALIZATION = false;
const bool VTK_ANG_VISUALIZATION = true;
const bool HERMES_SCAL_VISUALIZATION = true;
const bool VTK_SCAL_VISUALIZATION = true;
const bool HERMES_MESH_VISUALIZATION = false;

const bool MULTIMESH = false;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 5;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 1;

const int N = 6;          
const int M = N*(N+2)/2;

//
// Parameters for the source iteration via the Picard's method.
//
// Number of last iterations used.
// 1... standard fixed point.
// >1... Anderson acceleration.
const int PICARD_NUM_LAST_ITER_USED = 3;
// 0 <= beta <= 1... parameter for the Anderson acceleration.
const double PICARD_ANDERSON_BETA = 0;
// Stopping criterion for the Picard's method.
const double PICARD_TOL = 5e-3;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 1000; 
// Value for constant initial condition.
const double INIT_COND_CONST = 1.0; 

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;

int main(int argc, char* args[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 1);
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::matrixSolverType, matrix_solver_type);
  //Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  Hermes::vector<MeshSharedPtr > meshes;
  for (int i = 0; i < M; i++)
    meshes.push_back(new Mesh());
  
  MeshReaderH2D mloader;
  mloader.load(mesh_file.c_str(), meshes[0]);
  
  if (MULTIMESH)
  {  
    for (int i = 1; i < M; i++) 
    {
      // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
      meshes[i]->copy(meshes[0]);
          
      // Initial uniform refinements.
      for (int j = 0; j < INIT_REF_NUM; j++) 
        meshes[i]->refine_all_elements();
    }
  }
  for (int j = 0; j < INIT_REF_NUM; j++) 
    meshes[0]->refine_all_elements();

  if (HERMES_MESH_VISUALIZATION)
  {
    MeshView mview("Coarse mesh", new WinGeom(0, 0, 440, 350));
    mview.show(meshes[0]);
    Views::View::wait();
  }
  
  Loggable::Static::info("%d elements, %d vertices", meshes[0]->get_num_active_elements(), meshes[0]->get_num_vertex_nodes() );
  
  Hermes::vector<Solution<double>* > slns;
  for (int i = 0; i < M; i++)
    slns.push_back(new Solution<double>(MULTIMESH ? meshes[i] : meshes[0]));
 
  // Initialize the weak formulation.
  Hermes::vector<std::string> reflective_boundaries;
  Hermes::vector<std::string> vacuum_boundaries(1);
  vacuum_boundaries.push_back("vacuum");
  
  if (argc > 1)
  {
    bool assemble_matrix = strcmp(args[1], "Q");
    bool assemble_Q = !strcmp(args[1], "Q") || !strcmp(args[1], "LQ");
    
    WeakForm<double> *wf;
    Hermes::vector<const Space<double> *> spaces;
    
    if (!strcmp(args[1], "S") || !strcmp(args[1], "F"))
    {
      wf = new IsotropicScatteringMatrixForm(extinction, thermalization);
      spaces.push_back(new L2Space<double>(meshes[0], P_INIT));
    }
    else
    {
      wf = new SNWeakForm(N, extinction, thermalization, reflective_boundaries, vacuum_boundaries, args[1]);  
      for (int n = 0; n < M; n++)
        spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
    }
    
    Loggable::Static::info("Saving %s. NDOF = %d", args[1], Space<double>::get_num_dofs(spaces));
    
    DiscreteProblem<double> dp(wf, spaces);
    if (P_INIT == 0) dp.set_fvm();
    SourceIteration solver(&dp);
    
    // Perform the source iteration (by Picard's method with Anderson acceleration).
    solver.set_picard_max_iter(1);
    
    if (assemble_Q)  
    {
      solver.output_rhs(1);
      solver.set_rhs_E_matrix_dump_format(DF_HERMES_BIN);;
      solver.set_rhs_filename("Q");
      solver.set_rhs_number_format("%1.15f");
      solver.set_rhs_varname("Q");
    }
    if (assemble_matrix)
    {
      solver.output_matrix(1);
      solver.set_matrix_E_matrix_dump_format(DF_HERMES_BIN);
      solver.set_matrix_filename(!strcmp(args[1], "LQ") ? "L" : args[1]);
      solver.set_matrix_number_format("%1.15f");
      solver.set_matrix_varname(!strcmp(args[1], "LQ") ? "L" : args[1]);
    }
    
    try 
    { 
      solver.solve(); 
    } 
    catch(std::exception& e) { }
    
    delete wf;
    
    return 0;
  }
  
  SNWeakForm wf(N, extinction, thermalization, reflective_boundaries, vacuum_boundaries);
 
  // Initialize the FE problem.
  
  // Approximation spaces.
  Hermes::vector<const Space<double> *> spaces;
  
  for (int n = 0; n < M; n++)
    spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
  
  int ndof =  Space<double>::get_num_dofs(spaces);
  
  // Display the mesh.
//  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
//  oview.show(spaces[0]);
//  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
//  bview.show(spaces[0]);

  // Discrete formulation.
  DiscreteProblem<double> dp(&wf, spaces);
  if (P_INIT == 0) dp.set_fvm();
  
  // Algebraic solver.
  SourceIteration solver(&dp);
  
  solver.use_Anderson_acceleration(false);
  solver.set_picard_tol(PICARD_TOL);
  solver.set_picard_max_iter(PICARD_MAX_ITER);
  solver.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
  solver.set_anderson_beta(PICARD_ANDERSON_BETA);
  solver.set_verbose_output(true);
  
  Loggable::Static::info("Solving. NDOF = %d", ndof);
  cpu_time.tick();
  
  try
  {
    solver.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
  
  double *sln_vector = solver.get_sln_vector();
  Solution<double>::vector_to_solutions(sln_vector, spaces, slns);
  
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());

  // View the coarse mesh solution and/or save it to .vtk files.
  for (int n = 0; n < M; n++)
  {
    if(HERMES_ANG_VISUALIZATION)
    {
      ScalarView view("Solution", new WinGeom(0+450*n, 400, 450, 350));
      view.fix_scale_width(60);
      view.show(slns[n]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_ANG_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(slns[n], (std::string("sln_") + itos(n) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
  // View the coarse mesh scalar flux and/or save it to .vtk files.
  
  Hermes::vector<MeshFunction<double>*> scalar_fluxes;
  SupportClasses::MomentFilter::get_scalar_fluxes(slns, &scalar_fluxes, 1, wf.get_ordinates_data());
  
  if (HERMES_SCAL_VISUALIZATION)
  {
    ScalarView view("Scalar flux", new WinGeom(0, 0, 450, 350));
    view.fix_scale_width(60);
    view.show(scalar_fluxes[0]);
    Views::View::wait();
  }
  
  // VTK output.
  if(VTK_SCAL_VISUALIZATION)
  {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = false;
    lin.save_solution_vtk(scalar_fluxes[0], "scalar_flux.vtk", "Solution", mode_3D);
  }

  SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    
  if (SAVE_SLN_VECTOR)
  {
    std::string file = "x-R"+itos(INIT_REF_NUM)+"P"+itos(P_INIT)+"-S"+itos(N)+".dat";

    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the solution vector to %s", file.c_str());

    fs << setprecision(16);
    std::copy(sln_vector, sln_vector+ndof, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
  }
  
  return 0;
}
