#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "../problem_data.h"
#include "weakforms_neutronics.h"
#include <iomanip>
#include <iterator>


const bool SAVE_FLUX_PROFILE = true;
const bool SAVE_SLN_VECTOR = false;


const bool VTK_VISUALIZATION = true;
    const bool VTK_ANG_VISUALIZATION = true;
    const bool VTK_SCAL_VISUALIZATION = true;
    const bool VTK_ORD_VISUALIZATION = true;

const bool HERMES_VISUALIZATION = false;
    const bool HERMES_ANG_VISUALIZATION = false;
    const bool HERMES_SCAL_VISUALIZATION = false;
    const bool HERMES_ORD_VISUALIZATION = false;
    const bool HERMES_MESH_VISUALIZATION = false;
    const bool SHOW_INTERMEDIATE_ORDERS = false;     // Set to "true" to display coarse mesh solutions during adaptivity.
    const bool SHOW_INTERMEDIATE_SOLUTIONS = false;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

const bool MULTIMESH = true;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;
// Use vertex based limiter for P_INIT < 3?
const bool USE_LIMITER = false;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
const int N = 8;
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
const double PICARD_TOL = 1e-6;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 1000;
// Value for constant initial condition.
const double INIT_COND_CONST = 1.0;

//
// Error calculation & adaptivity.
//
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.5;
// Error calculator.
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, M);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-2;

// Setup the convergence graphs.
void setup_convergence_graph(GnuplotGraph *graph)
{
  graph->add_row("angular", "r", "-", "o");
//  graph->add_row("scalar", "b", "-", "+");

  graph->set_log_x();
  graph->set_log_y();
  graph->show_legend();
  graph->show_grid();
}

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;

int main(int argc, char* args[])
{
  // Set the number of threads used in Hermes.
  
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::matrixSolverType, matrix_solver_type);
  //Hermes::HermesCommonApi.set_integral_param_value(Hermes::numThreads, 2);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  Hermes::vector<MeshSharedPtr > meshes;
  for (int i = 0; i < M; i++)
    meshes.push_back(MeshSharedPtr(new Mesh()));
  
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
  
  Hermes::vector<MeshFunctionSharedPtr<double> > slns, ref_slns;
  for (int i = 0; i < M; i++)
  {
    ref_slns.push_back(new Solution<double>);
    slns.push_back(new ConstantSolution<double>(MULTIMESH ? meshes[i] : meshes[0], INIT_COND_CONST));
  }

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
  
  SNWeakForm wf(N, matprop, reflective_boundaries, vacuum_boundaries);
 
  SupportClasses::Visualization views(N_GROUPS, N, wf.get_ordinates_data());

  // Approximation spaces.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  
  for (int n = 0; n < M; n++)
    spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
    //spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT, new L2ShapesetTaylor));
  
  int ndof =  Space<double>::get_num_dofs(spaces);
  
  // Display the mesh.
//  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
//  oview.show(spaces[0]);
//  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
//  bview.show(spaces[0]);

  // Discrete problem solver.
  SourceIteration solver;
  solver.set_weak_formulation(&wf);

  solver.use_Anderson_acceleration(false);
  solver.set_tolerance(PICARD_TOL, Hermes::Solvers::ResidualNormRatioToInitial);
  solver.set_max_allowed_iterations(PICARD_MAX_ITER);
  solver.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
  solver.set_anderson_beta(PICARD_ANDERSON_BETA);
  solver.set_verbose_output(true);

  // Initialize refinement selectors.
  L2ProjBasedSelector<double> selector(CAND_LIST);
  selector.set_error_weights(1., 1., 1.);
  
  Hermes::vector<RefinementSelectors::Selector<double> *> selectors(M);
  for (int i=0; i < M; i++)
    selectors.push_back(&selector);

  // DOF and CPU convergence graphs.
  GnuplotGraph graph_dof_est("Relative L2 error of flux approximation", "NDOF", "error [%]");
  GnuplotGraph graph_cpu_est("Relative L2 error of flux approximation", "CPU time [s]", "error [%]");
  setup_convergence_graph(&graph_dof_est);
  setup_convergence_graph(&graph_cpu_est);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Hermes::vector<SpaceSharedPtr<double> > ref_spaces;
  ref_spaces.resize(M);
  double *sln_vector;
  do
  {
    Loggable::Static::info("---- Adaptivity step %d:", as);

    // Initialize the fine mesh problem.

    for (unsigned int i = 0; i < M; i++)
    {
      Mesh::ReferenceMeshCreator ref_mesh_creator(meshes[i]);
      MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator ref_space_creator(spaces[i], ref_mesh);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
      ref_spaces[i] = ref_space;
    }

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces);

    report_num_dof("Solving on reference spaces, #DOF: ", ref_spaces);

    // Time measurement.
    cpu_time.tick();

    // Iterate on reference spaces.
    solver.set_spaces(ref_spaces);

    try
    {
      solver.solve(slns);
    }
    catch (Hermes::Exceptions::Exception& e)
    {
      std::cout << e.info();
      return -1;
    }
    catch (std::exception& e)
    {
      std::cout << e.what();
      return -1;
    }

    sln_vector = solver.get_sln_vector();

    if (P_INIT < 3 && USE_LIMITER)
    {
      PostProcessing::VertexBasedLimiter limiter(ref_spaces, sln_vector, P_INIT);
      limiter.get_solutions(ref_slns);
    }
    else
    {
      Solution<double>::vector_to_solutions(sln_vector, ref_spaces, ref_slns);
    }

    Loggable::Static::info("Projecting reference solutions on coarse mesh.");
    OGProjection<double> ogProjection; ogProjection.project_global(spaces, ref_slns, slns);

    adaptivity.set_spaces(spaces);

    errorCalculator.calculate_errors(slns, ref_slns);
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    // Add the results to convergence graphs.
    cpu_time.tick();

//    Hermes::vector<MeshFunctionSharedPtr<double> > scalar_fluxes, ref_scalar_fluxes;
//    SupportClasses::MomentFilter::get_scalar_fluxes_with_derivatives(slns, &scalar_fluxes, N_GROUPS, wf.get_ordinates_data());
//    SupportClasses::MomentFilter::get_scalar_fluxes_with_derivatives(ref_slns, &ref_scalar_fluxes, N_GROUPS, wf.get_ordinates_data());

//    errorCalculator.calculate_errors(scalar_fluxes, ref_scalar_fluxes, false);
//    double scalar_flux_err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    graph_dof_est.add_values(0, ndof_ref, err_est_rel_total);
    graph_cpu_est.add_values(0, cpu_time.accumulated(), err_est_rel_total);
 //   graph_dof_est.add_values(1, ndof_ref, scalar_flux_err_est_rel_total);
 //   graph_cpu_est.add_values(1, cpu_time.accumulated(), scalar_flux_err_est_rel_total);

    cpu_time.tick(TimeMeasurable::HERMES_SKIP);

    std::cout << "Discretization error (angular fluxes): " << err_est_rel_total << "%." << std::endl;
//    std::cout << "Discretization error (scalar fluxes): " << scalar_flux_err_est_rel_total << "%." << std::endl;

    // If err_est too large, adapt spaces
    if (err_est_rel_total < ERR_STOP)
      done = true;
    else
    {
      Loggable::Static::info("Adapting coarse mesh.");

      try
      {
        done = adaptivity.adapt(selectors);
      }
      catch (Hermes::Exceptions::Exception& e)
      {
        std::cout << e.info();
        return -1;
      }
      catch (std::exception& e)
      {
        std::cout << e.what();
        return -1;
      }

    }

    if (!done) // Increase the adaptivity step counter.
      as++;
    else
    {
      cpu_time.tick();
      Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
      cpu_time.reset();

      // Visualization.

      if (HERMES_VISUALIZATION)
      {
        if (HERMES_ANG_VISUALIZATION)
            views.show_solutions(ref_slns);
        if (HERMES_SCAL_VISUALIZATION)
            views.show_scalar_fluxes(ref_slns);
        if (HERMES_ORD_VISUALIZATION)
            views.show_orders(spaces);

        // Wait for the view to be closed.
        Views::View::wait();
      }
      if (VTK_VISUALIZATION)
      {
        if (VTK_ANG_VISUALIZATION)
            views.save_solutions_vtk("flux", "flux", ref_slns);
        if (VTK_SCAL_VISUALIZATION)
            views.save_scalar_fluxes_vtk("scal_flux", "scalar_flux", ref_slns);
        if (VTK_ORD_VISUALIZATION)
            views.save_orders_vtk("space", spaces);
      }

      // Make the fine-mesh spaces the final spaces for further analyses.
      spaces = ref_spaces;
    }
  }
  while (done == false);

  // Save the convergence graphs.
  graph_dof_est.save(("conv_dof_est-R"+tostr(INIT_REF_NUM)+"P"+tostr(P_INIT)+"-S"+tostr(N)+".gp").c_str());
  graph_cpu_est.save(("conv_cpu_est-R"+tostr(INIT_REF_NUM)+"P"+tostr(P_INIT)+"-S"+tostr(N)+".gp").c_str());

  cpu_time.tick();
  







  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());

  if (SAVE_SLN_VECTOR)
  {
    std::string file = "x-R"+tostr(INIT_REF_NUM)+"P"+tostr(P_INIT)+"-S"+tostr(N)+".dat";

    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the solution vector to %s", file.c_str());

    fs << std::setprecision(16);
    std::copy(sln_vector, sln_vector+ndof, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
  }
  
  if (SAVE_FLUX_PROFILE)
  {
      Hermes::vector<MeshFunctionSharedPtr<double> > scalar_fluxes;
      SupportClasses::MomentFilter::get_scalar_fluxes(slns, &scalar_fluxes, N_GROUPS, wf.get_ordinates_data());

    std::string file = "flux_x_0.1667-R"+tostr(INIT_REF_NUM)+"P"+tostr(P_INIT)+"-S"+tostr(N)+".dat";

    double y = 0.1667;
    
    int nintervals = 100;
    double a = 10.;
    double dx = a / nintervals; // advance by 1 mm
    int npts = nintervals + 1;
    
    double *res = new double [(npts)*N_GROUPS];
        
    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the scalar flux profile at y=0.1667cm to %s", file.c_str());
    
    fs << std::setprecision(16);
    
    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::cout << std::endl << "GROUP " << g << std::endl;
      
      double x = 0;
      
      for (int i = 0; i < npts; i++, x+=dx)
      {
        res[i+g*npts] = *scalar_fluxes[g]->get_pt_value(x, y)->val;    
        (std::cout << "(" << x << "," << y << ") = " << res[i+g*npts] << std::endl).flush();
      }
      
      std::copy(res+g*(npts), res+(g+1)*(npts), std::ostream_iterator<double>(fs, "\n"));
    }
    
    fs.close();
    delete [] res;
  }
  
  return 0;
}
