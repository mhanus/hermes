#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;
using namespace WeakFormsH1;
using namespace RefinementSelectors;

// Weak forms.
#include "weakform_library/weakforms_neutronics.h"

//  This example is a standard nuclear engineering benchmark describing an external-force-driven
//  configuration without fissile materials present, using one-group neutron diffusion approximation.
//  It is very similar to example "saphir", the main difference being that the mesh is loaded in
//  the ExodusII format (created for example by Cubit).
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources.
//
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh).
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary").
//       Zero Neumann for the left and bottom edges ("reflection boundary").
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
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
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 0.01;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
double L = 30;                           // Edge of square.
double L0 = 0.75*0.5*L;                  // End of first water region.
double L1 = 0.5*L;                       // End of second water region.
double L2 = 0.75*L;                      // End of iron region.
double Q_EXT = 1.0;                      // Neutron source (nonzero in region 1 only).
double SIGMA_T_WATER = 3.33;             // Total cross-section.
double SIGMA_T_IRON = 1.33;
double C_WATER = 0.994;                  // Scattering ratio.
double C_IRON = 0.831;
double D_WATER = 1./(3.*SIGMA_T_WATER);  // Diffusion coefficient.
double D_IRON = 1./(3.*SIGMA_T_IRON);
double SIGMA_A_WATER = SIGMA_T_WATER - C_WATER*SIGMA_T_WATER;  // Absorbing cross-section.
double SIGMA_A_IRON = SIGMA_T_IRON - C_IRON*SIGMA_T_IRON;

// Materials. 
// Mesh file generated by Cubit - do not replace the markers.
const std::string WATER_1 = "1";
const std::string WATER_2 = "2";
const std::string IRON = "3";

// Boundaries.
const std::string ZERO_FLUX_BOUNDARY = "2";

int main(int argc, char* argv[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::numThreads, 2);
  
  // Load the mesh.
  Mesh mesh;
  MeshReaderExodusII mloader;
  if (!mloader.load("iron-water.e", &mesh)) Neutronics::ErrorHandling::error_function("ExodusII mesh load failed.");
   
  // Perform initial uniform mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  
  Views::MeshView mv("Mesh", new Views::WinGeom(0, 0, 440, 350));
  mv.show(&mesh);
  Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
  
  // Set essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(ZERO_FLUX_BOUNDARY, 0.0);
  EssentialBCs<double> bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  
  // Associate element markers (corresponding to physical regions) 
  // with material properties (diffusion coefficient, absorption 
  // cross-section, external sources).
  Hermes::vector<std::string> regions(WATER_1, WATER_2, IRON);
  Hermes::vector<double> D_map(D_WATER, D_WATER, D_IRON);
  Hermes::vector<double> Sigma_a_map(SIGMA_A_WATER, SIGMA_A_WATER, SIGMA_A_IRON);
  Hermes::vector<double> Sources_map(Q_EXT, 0.0, 0.0);
  
  // Initialize the weak formulation.
  Neutronics::SimpleMonoenergeticDiffusionWeakForms::FixedSourceProblem
    wf(regions, D_map, Sigma_a_map, Sources_map);
  
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(450, 0, 400, 350));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution<double> ref_sln;
  do
  {
    Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    SpaceSharedPtr<double> ref_space = Space<double>::construct_refined_space(&space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize the FE problem.
    DiscreteProblem<double> dp(&wf, ref_space);

    // Perform Newton's iteration on reference emesh.
    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(false);
    
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      Neutronics::ErrorHandling::error_function("Newton's iteration failed.");
    }
    
    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    Solution<double> sln;
    Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double> ogProjection;
    ogProjection.project_global(&space, &ref_sln, &sln); 

    // Time measurement.
    cpu_time.tick();
   
    // View the coarse mesh solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);

    // Skip visualization time.
    cpu_time.tick(TimeMeasurable::HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Loggable::Static::info("Calculating error estimate."); 
    Adapt<double> adaptivity(&space);
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln) * 100;

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
      delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);
  
  cpu_time.tick();
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
