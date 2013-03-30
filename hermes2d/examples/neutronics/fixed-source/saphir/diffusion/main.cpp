//
//  IAEA EIR-2 benchmark problem, diffusion approximation.
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//       where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//       and Q_{ext}(x,y) external sources.
//
//  Domain: square (0, 96)x(0, 96) (see mesh file ../domain.mesh).
//
//  BC:  Zero Dirichlet on all edges (approximate version of the void boundary condition).
//
#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;
using namespace WeakFormsH1;

// Weak forms.
#include "weakform_library/weakforms_neutronics.h"

using namespace RefinementSelectors;

const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.6;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters
double SIGMA_T_1 = 0.60;                          // Total cross-sections.
double SIGMA_T_2 = 0.48;
double SIGMA_T_3 = 0.70;
double SIGMA_T_4 = 0.65;
double SIGMA_T_5 = 0.90;
double SIGMA_S_1 = 0.53;                          // Scattering cross sections.
double SIGMA_S_2 = 0.20;
double SIGMA_S_3 = 0.66;
double SIGMA_S_4 = 0.50;
double SIGMA_S_5 = 0.89;
double Q_EXT_1 = 1;                               // Nonzero sources in regions 1 and 3 only,
double Q_EXT_3 = 1;                               // sources in other regions are zero.

// Additional constants
double D_1 = 1/(3.*SIGMA_T_1);                    // Diffusion coefficients.
double D_2 = 1/(3.*SIGMA_T_2);
double D_3 = 1/(3.*SIGMA_T_3);
double D_4 = 1/(3.*SIGMA_T_4);
double D_5 = 1/(3.*SIGMA_T_5);
double SIGMA_A_1 = SIGMA_T_1 - SIGMA_S_1;         // Absorption coefficients.
double SIGMA_A_2 = SIGMA_T_2 - SIGMA_S_2;
double SIGMA_A_3 = SIGMA_T_3 - SIGMA_S_3;
double SIGMA_A_4 = SIGMA_T_4 - SIGMA_S_4;
double SIGMA_A_5 = SIGMA_T_5 - SIGMA_S_5;

int main(int argc, char* argv[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mesh_reader;
  mesh_reader.load("../domain.mesh", &mesh);
  
  // Perform initial uniform mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("bottom", "right", "top", "left"), 0.0);
  EssentialBCs<double> bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Associate element markers (corresponding to physical regions) 
  // with material properties (diffusion coefficient, absorption 
  // cross-section, external sources).
  Hermes::vector<std::string> regions("1", "2", "3", "4", "5");
  Hermes::vector<double> D_map(D_1, D_2, D_3, D_4, D_5);
  Hermes::vector<double> Sigma_a_map(SIGMA_A_1, SIGMA_A_2, SIGMA_A_3, SIGMA_A_4, SIGMA_A_5);
  Hermes::vector<double> Sources_map(Q_EXT_1, 0.0, Q_EXT_3, 0.0, 0.0);
  
  // Initialize the weak formulation.
  Neutronics::SimpleMonoenergeticDiffusionWeakForms::FixedSourceProblem
    wf(regions, D_map, Sigma_a_map, Sources_map);

  // Initialize coarse and reference mesh solution.
  Solution<double> sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
#ifndef NOGLUT  
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(450, 0, 400, 350));
#endif  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Create Newton's solver on reference mesh.
  NewtonSolver<double> newton(&dp);
  newton.set_verbose_output(false);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(&mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(&space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

    int ndof_ref = ref_space->get_num_dofs();
    newton.set_space(ref_space);
    
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      Neutronics::ErrorHandling::error_function("Newton's iteration failed.");
    }
      
    // Translate the resulting coefficient vector into instances of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);
      
    // Project the fine mesh solution onto the coarse mesh.
    Solution<double> sln;
    Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double> ogProjection;
    ogProjection.project_global(&space, &ref_sln, &sln); 

    // Time measurement.
    cpu_time.tick();
   
    // View the coarse mesh solution and polynomial orders.
#ifndef NOGLUT
    sview.show(&sln);
    oview.show(&space);
#endif
    // Skip visualization time.
    cpu_time.tick(TimeMeasurable::HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Loggable::Static::info("Calculating error estimate."); 
    Adapt<double> adaptivity(&space);
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln, true, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

    // Report results.
    Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
         space.get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    
    // Skip graphing time.
    cpu_time.tick(TimeMeasurable::HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (space.get_num_dofs() >= NDOF_STOP) done = true;
    
    if(done == false) 
      delete ref_mesh;
    delete ref_space;
  }
  while (done == false);
  
  cpu_time.tick();
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

#ifndef NOGLUT
  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();
#else
  // Output solution in VTK format.
  Hermes::Hermes2D::Views::Linearizer lin;
  bool mode_3D = true;
  lin.save_solution_vtk(&ref_sln, "sln.vtk", "Neutron flux", mode_3D);

  // Output mesh and element orders in VTK format.
  Hermes::Hermes2D::Views::Orderizer ord;
  ord.save_orders_vtk(&space, "ord.vtk");
#endif
  
  return 0;
}
