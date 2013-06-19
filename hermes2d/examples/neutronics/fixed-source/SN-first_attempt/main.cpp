#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"
#include "weakforms_neutronics.h"


const bool HERMES_ANG_VISUALIZATION = false;
const bool VTK_ANG_VISUALIZATION = true;
const bool HERMES_SCAL_VISUALIZATION = true;
const bool VTK_SCAL_VISUALIZATION = true;


const bool MULTIMESH = false;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
const int N = 6;                    
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
const double PICARD_TOL = 5e-3;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 150; 
// Value for constant initial condition.
const double INIT_COND_CONST = 1.0; 

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;
// Only for matrix_solver_type == SOLVER_AZTECOO :
const char* iterative_method = "bicgstab";
const char* preconditioner = "jacobi";

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

  MeshView mview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  mview.show(meshes[0]);
  Views::View::wait();
  
  Loggable::Static::info("%d elements, %d vertices", meshes[0]->get_num_active_elements(), meshes[0]->get_num_vertex_nodes() );
  
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  
  for (int n = 0; n < M; n++)
    spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
  
  int ndof =  Space<double>::get_num_dofs(spaces);

  // Display the mesh.
//  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
//  oview.show(spaces[0]);
//  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
//  bview.show(spaces[0]);

  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for (int i = 0; i < M; i++)
    slns.push_back(new ConstantSolution<double>(MULTIMESH ? meshes[i] : meshes[0], INIT_COND_CONST));
  
    // Load material data.
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, rm_map);
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_sn(Ssn);
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;

  // Initialize the weak formulation.
  Hermes::vector<std::string> reflective_boundaries;
  reflective_boundaries.push_back("reflective");
  SNWeakForm wf(N, matprop, reflective_boundaries);
 
  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);
  SourceIteration solver(&dp);
  
  solver.use_Anderson_acceleration(false);
  
  Loggable::Static::info("Solving. NDOF = %d", ndof);
  cpu_time.tick();

  // Perform the source iteration (by Picard's method with Anderson acceleration).
  solver.set_tolerance(PICARD_TOL);
  solver.set_max_allowed_iterations(PICARD_MAX_ITER);
  solver.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
  solver.set_anderson_beta(PICARD_ANDERSON_BETA);
  solver.set_verbose_output(true);
  try
  {
    solver.solve(slns);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, slns);
 
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());
  
  // Wait for keyboard or mouse input.
  // View the coarse mesh solution.
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
      lin.save_solution_vtk(slns[n], (std::string("sln_") + tostr(n) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
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
      bool mode_3D = false;
      lin.save_solution_vtk(scalar_fluxes[g], (std::string("scalar_flux_g_") + tostr(g) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
  
  
  
  
  // Integrate absorption rates and scalar fluxes over specified edit regions and compare with results
  // from the collision probabilities code DRAGON (see problem_data.h).

  PostProcessor pp(NEUTRONICS_SN, HERMES_PLANAR, &wf.get_ordinates_data());
    
  Hermes::vector<double> absorption_rates, integrated_fluxes, areas;
  Hermes::vector<double> absorption_rates_cmp, integrated_fluxes_cmp, areas_cmp;
  pp.get_integrated_reaction_rates(ABSORPTION, slns, &absorption_rates, matprop, edit_regions);
  pp.get_integrated_scalar_fluxes(slns, &integrated_fluxes, N_GROUPS, edit_regions);
  pp.get_areas(meshes[0], edit_regions, &areas); // Areas of the edit regions.
  
  
  
  
  // Multiply integral results by the number of times each region appears in the assembly.
  absorption_rates_cmp.push_back(absorption_rates[0] * 4);
  integrated_fluxes_cmp.push_back(integrated_fluxes[0] * 4);
  areas_cmp.push_back(areas[0] * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[1] + absorption_rates[4]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[1] + integrated_fluxes[4]) * 4);
  areas_cmp.push_back((areas[1] + areas[4]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[2] + absorption_rates[5]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[2] + integrated_fluxes[5]) * 4);
  areas_cmp.push_back((areas[2] + areas[5]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[3] + absorption_rates[6]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[3] + integrated_fluxes[6]) * 4);
  areas_cmp.push_back((areas[3] + areas[6]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[7]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[7]) * 4);
  areas_cmp.push_back((areas[7]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[8] + absorption_rates[10]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[8] + integrated_fluxes[10]) * 4);
  areas_cmp.push_back((areas[8] + areas[10]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[9] + absorption_rates[13]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[9] + integrated_fluxes[13]) * 4);
  areas_cmp.push_back((areas[9] + areas[13]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[11]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[11]) * 4);
  areas_cmp.push_back((areas[11]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[12] + absorption_rates[14]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[12] + integrated_fluxes[14]) * 4);
  areas_cmp.push_back((areas[12] + areas[14]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[15]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[15]) * 4);
  areas_cmp.push_back((areas[15]) * 4);
  
  
  
  
  
  
  
  
  absorption_rates_cmp.push_back(absorption_rates[16] * 4);
  integrated_fluxes_cmp.push_back(integrated_fluxes[16] * 4);
  areas_cmp.push_back(areas[16] * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[17] + absorption_rates[20]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[17] + integrated_fluxes[20]) * 4);
  areas_cmp.push_back((areas[17] + areas[20]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[18] + absorption_rates[21]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[18] + integrated_fluxes[21]) * 4);
  areas_cmp.push_back((areas[18] + areas[21]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[19] + absorption_rates[22]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[19] + integrated_fluxes[22]) * 4);
  areas_cmp.push_back((areas[19] + areas[22]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[23]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[23]) * 4);
  areas_cmp.push_back((areas[23]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[24] + absorption_rates[26]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[24] + integrated_fluxes[26]) * 4);
  areas_cmp.push_back((areas[24] + areas[26]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[25] + absorption_rates[29]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[25] + integrated_fluxes[29]) * 4);
  areas_cmp.push_back((areas[25] + areas[29]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[27]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[27]) * 4);
  areas_cmp.push_back((areas[27]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[28] + absorption_rates[30]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[28] + integrated_fluxes[30]) * 4);
  areas_cmp.push_back((areas[28] + areas[30]) * 4);
  
  absorption_rates_cmp.push_back((absorption_rates[31]) * 4);
  integrated_fluxes_cmp.push_back((integrated_fluxes[31]) * 4);
  areas_cmp.push_back((areas[31]) * 4);
  
  
  
  
  for (int i = 0; i < 2*n_pins; i++)
    Loggable::Static::info("Absorption rate integrated over (DRAGON reg. #%d) = %f (error %g%%)", i+1, absorption_rates_cmp[i], 
        fabs(absorption_rates_cmp[i] - ref_integrated_absorption_rates[i])/ref_integrated_absorption_rates[i] * 100);
  for (int i = 0; i < 2*n_pins; i++)
    Loggable::Static::info("Scalar flux integrated over (DRAGON reg. #%d) = %f (error %g%%)", i+1, integrated_fluxes_cmp[i],
        fabs(integrated_fluxes_cmp[i] - ref_integrated_fluxes[i])/ref_integrated_fluxes[i] * 100);
  for (int i = 0; i < 2*n_pins; i++)
    Loggable::Static::info("Area of (DRAGON reg. #%d) = %f (error %g%%)", i+1, areas_cmp[i],
        fabs(areas_cmp[i] - ref_regions_areas[i])/ref_regions_areas[i] * 100);

  
  return 0;
}
