#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"
#include "weakforms_neutronics.h"

#include <iterator>
const bool SAVE_FLUX_PROFILE = true;
const bool SAVE_SLN_VECTOR = true;

const bool HERMES_ANG_VISUALIZATION = false;
const bool VTK_ANG_VISUALIZATION = true;
const bool HERMES_SCAL_VISUALIZATION = true;
const bool VTK_SCAL_VISUALIZATION = true;
const bool HERMES_MESH_VISUALIZATION = false;

const bool MULTIMESH = false;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 4;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
const int N = 16;                    
const int M = N_GROUPS * N*(N+2)/2;

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
const double PICARD_TOL = 1e-9;
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
  
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for (int i = 0; i < M; i++)
    slns.push_back(new ConstantSolution<double>(MULTIMESH ? meshes[i] : meshes[0], INIT_COND_CONST));
  
    // Load material data.
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, std::set<std::string>(materials.begin(), materials.end()));
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_sn(Ssn);
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;

  // Initialize the weak formulation.
  Hermes::vector<std::string> reflective_boundaries(1);
  reflective_boundaries.push_back("reflective");
  Hermes::vector<std::string> vacuum_boundaries(1);
  vacuum_boundaries.push_back("vacuum");
  
  if (argc > 1)
  {
    bool assemble_matrix = strcmp(args[1], "Q");
    bool assemble_Q = !strcmp(args[1], "Q") || !strcmp(args[1], "LQ");
    
    WeakForm<double> *wf;
    Hermes::vector<SpaceSharedPtr<double> > spaces;
    
    if (!strcmp(args[1], "S") || !strcmp(args[1], "F"))
    {
      wf = new IsotropicScatteringAndFissionMatrixForms(matprop, args[1]);
      SupportClasses::AngleGroupFlattener ag(N_GROUPS);
      for (int g = 0; g < N_GROUPS; g++)
        spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[ag(0,g)] : meshes[0], P_INIT));
    }
    else
    {
      wf = new SNWeakForm(N, matprop, reflective_boundaries, vacuum_boundaries, args[1]);  
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
  
  SNWeakForm wf(N, matprop, reflective_boundaries, vacuum_boundaries);
 
  // Initialize the FE problem.
  
  // Approximation spaces.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  
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
    solver.solve(slns);
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
  
  Hermes::vector<MeshFunctionSharedPtr<double> > scalar_fluxes;
  SupportClasses::MomentFilter::get_scalar_fluxes(slns, &scalar_fluxes, N_GROUPS, wf.get_ordinates_data());
  
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    if (HERMES_SCAL_VISUALIZATION)
    {
      ScalarView view("Scalar flux", new WinGeom(0+450*g, 0, 450, 350));
      view.fix_scale_width(60);
      view.show(scalar_fluxes[g]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_SCAL_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin;
      bool mode_3D = true;
      lin.save_solution_vtk(scalar_fluxes[g], (std::string("scalar_flux_g_") + itos(g) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
    
  if (SAVE_SLN_VECTOR)
  {
    std::string file = "x-R"+itos(INIT_REF_NUM)+"P"+itos(P_INIT)+"-S"+itos(N)+".dat";

    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the solution vector to %s", file.c_str());

    fs << setprecision(16);
    std::copy(sln_vector, sln_vector+ndof, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
  }
  
  if (SAVE_FLUX_PROFILE)
  {
    std::string file = "flux_5.625_y-R"+itos(INIT_REF_NUM)+"P"+itos(P_INIT)+"-S"+itos(N)+".dat";

    double x = 5.625;
    
    int nintervals = 100;
    double a = 10.;
    double dy = a / nintervals; // advance by 1 mm
    int npts = nintervals + 1;
    
    double *res = new double [(npts)*N_GROUPS];
        
    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the scalar flux profile at x=5.625cm to %s", file.c_str());
    
    fs << std::setprecision(16);
    
    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::cout << std::endl << "GROUP " << g << std::endl;
      
      double y = 0;
      
      for (int i = 0; i < npts; i++, y+=dy)
      {
        (std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = *scalar_fluxes[g]->get_pt_value(x, y)->val;    
      }
      
      std::copy(res+g*(npts), res+(g+1)*(npts), std::ostream_iterator<double>(fs, "\n"));
    }
    
    fs.close();
    delete [] res;
  }
    
  return 0;
}
