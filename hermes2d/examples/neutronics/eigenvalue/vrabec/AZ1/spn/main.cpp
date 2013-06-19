#define HERMES_REPORT_ALL
#include "../problem_data.h"
#include "definitions.h"
#include "../../../../utils.h"

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 2;  // Monoenergetic (single group) problem.
const unsigned int SPN_ORDER = 3; // SP3 approximation

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    0,  0,     // SP1
    0,  0,     // SP3
//     1,  1,     // SP5
//     2,  2,     // SP7
//     2,  2      // SP9
  };
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    1,  1,     // SP1
    1,  1,     // SP3
//     1,  1,     // SP5
//     2,  2,     // SP7
//     2,  2      // SP9
  };
        
const double THRESHOLD = 0.2;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = -1;                 // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  
                                                  
const VisualizationOptions visualization = HERMES_SCALAR_VISUALIZATION;
const bool HERMES_VISUALIZATION = visualization & (HERMES_SCALAR_VISUALIZATION | HERMES_ANGULAR_VISUALIZATION);
const bool VTK_VISUALIZATION = visualization & (VTK_SCALAR_VISUALIZATION | VTK_ANGULAR_VISUALIZATION);
const bool DISPLAY_MESHES = false;       // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool INTERMEDIATE_VISUALIZATION = true; // Set to "true" to display coarse mesh solutions during adaptivity.

const int SAVE_MATRICES = 0; // If non-zero, save algebraic representation of individual parts comprising the weak formulation.
                             // If SAVE_MATRICES == 2, the program ends right after saving the matrices.

//
// Eigenvalue iteration control.
//
double TOL_PIT_CM = 1e-5;   // Tolerance for convergence on the coarse mesh.
double TOL_PIT_FM = 1e-9;   // Tolerance for convergence on the fine mesh.
const int MAX_PIT = 1000;   // Maximal number of iterations.
const bool MEASURE_CONVERGENCE_BY_RESIDUAL = false;  // When 'false', eigenvalue difference will be used to monitor convergence. 

const bool USE_RAYLEIGH_QUOTIENT = false; // Use Rayleigh quotient to estimate the eigenvalue in each iteration. 
                                          // When 'false', the reciprocal norm of current eigenvector iterate will be used.

// Shifting strategy for the inverse iteration.
const KeffEigenvalueIteration::ShiftStrategies SHIFT_STRATEGY = KeffEigenvalueIteration::NO_SHIFT;
const double FIXED_SHIFT = 0.98;
const bool MODIFY_SHIFT_DURING_ADAPTIVITY = true;

int main(int argc, char* argv[])
{  
    // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 2);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  ProblemData data("../data.txt");
  MaterialProperties::MaterialPropertyMaps *matprop = new MaterialProperties::MaterialPropertyMaps(N_GROUPS, SPN_ORDER, data.mat_list);
    
  matprop->set_D(data.D);
  matprop->set_Sigma_f(data.Sf);
  matprop->set_nu(data.nu);
  matprop->set_Sigma_a(data.Sa);
  matprop->set_Sigma_s(data.Ss);
  matprop->set_fission_materials(fission_materials);
  
  matprop->validate();
  
  std::cout << *matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<MeshSharedPtr > meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load((std::string("../") + mesh_file).c_str(), meshes[0]);
  
  // Convert the mesh so that it has one type of elements (optional). 
  //meshes[0]->convert_quads_to_triangles();
  //meshes[0]->convert_triangles_to_quads();
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  
  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);
  
  // Initialize the weak formulation.
  // CustomWeakForm wf(*matprop, SPN_ORDER, HermesMultiArray<std::string>(bdy_vacuum));
  CustomWeakForm wf(*matprop, SPN_ORDER, Hermes::vector<std::string>());

  // Create pointers to solutions from the latest power iteration and approximation spaces with default shapeset.
  Hermes::vector<MeshFunctionSharedPtr<double> > power_iterates;  
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    spaces.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
    power_iterates.push_back(new Solution<double>());
  }
  
  if (argc == 3 && !strcmp(argv[1], "-sln"))
  {
    Neutronics::Common::SupportClasses::Visualization* vis;
    load_solution(argv[2], spaces, matprop, visualization);
    return 0;
  }
    
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  report_num_dof("Coarse mesh power iteration, NDOF: ", spaces);
  
  if (SHIFT_STRATEGY)
    wf.add_fission_sparse_structure();
    
    Neutronics::KeffEigenvalueIteration keff_eigenvalue_iteration(&wf, spaces);
  keff_eigenvalue_iteration.set_max_allowed_iterations(MAX_PIT);  
  keff_eigenvalue_iteration.measure_convergence_by_residual(MEASURE_CONVERGENCE_BY_RESIDUAL);
  
  if (MEASURE_CONVERGENCE_BY_RESIDUAL)
    keff_eigenvalue_iteration.set_tolerance(STRATEGY >= 0 ? TOL_PIT_CM : TOL_PIT_FM);
  else
    keff_eigenvalue_iteration.set_keff_tol(STRATEGY >= 0 ? TOL_PIT_CM : TOL_PIT_FM);
  
  WeakForm<double> prod_wf(wf.get_neq());
  if (USE_RAYLEIGH_QUOTIENT || SHIFT_STRATEGY || SAVE_MATRICES)
    wf.get_fission_yield_part(&prod_wf);
  if (USE_RAYLEIGH_QUOTIENT || SHIFT_STRATEGY)
    keff_eigenvalue_iteration.set_production_weakform(&prod_wf);
  
  if (SAVE_MATRICES)
  {
    WeakForm<double> diff_wf(wf.get_neq());
    WeakForm<double> scat_wf(wf.get_neq());
//    wf.get_diffusion_reaction_part(&diff_wf, HermesMultiArray<std::string>(bdy_vacuum));
//    wf.get_scattering_part(&scat_wf);
    
    // NOTE: If called after a call to 'solve', a 'false' argument could be added to prevent reassigning dofs
    
    save_algebraic_representation(&diff_wf, spaces, "A");
    save_algebraic_representation(&prod_wf, spaces, "B");
    save_algebraic_representation(&scat_wf, spaces, "S");

    if (SAVE_MATRICES == 2)
      return 0;
  }
/*  
  keff_eigenvalue_iteration.set_matrix_E_matrix_dump_format(Hermes::Algebra::DF_HERMES_BIN);
  keff_eigenvalue_iteration.set_matrix_filename("A");
  keff_eigenvalue_iteration.set_matrix_number_format("%1.15f");
  keff_eigenvalue_iteration.output_matrix(1);
  keff_eigenvalue_iteration.set_rhs_E_matrix_dump_format(Hermes::Algebra::DF_HERMES_BIN);
  keff_eigenvalue_iteration.set_rhs_filename("b");
  keff_eigenvalue_iteration.set_rhs_number_format("%1.15f");
  keff_eigenvalue_iteration.output_rhs(1);
*/    
  keff_eigenvalue_iteration.use_rayleigh_quotient(USE_RAYLEIGH_QUOTIENT);  
  keff_eigenvalue_iteration.set_spectral_shift_strategy(SHIFT_STRATEGY, FIXED_SHIFT);
    
  keff_eigenvalue_iteration.solve();
  Solution<double>::vector_to_solutions(keff_eigenvalue_iteration.get_sln_vector(), spaces, power_iterates);

  
  if (STRATEGY >= 0)
  {
    // DOF and CPU convergence graphs
    GnuplotGraph graph_dof("Error convergence", "NDOF", "log(error)");
    graph_dof.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_dof.add_row("L2 err. est. [%]", "g", "-", "s");
    graph_dof.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_dof.set_log_y();
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Error convergence", "CPU time [s]", "log(error)");
    graph_cpu.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_cpu.add_row("L2 err. est. [%]", "g", "-", "s");
    graph_cpu.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_cpu.set_log_y();
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    // Create pointers to coarse mesh solutions used for error estimation.
    Hermes::vector<MeshFunctionSharedPtr<double> > coarse_solutions;
    
    // Initialize all the new solution variables.
    for (unsigned int i = 0; i < N_EQUATIONS; i++) 
      coarse_solutions.push_back(new Solution<double>());
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector<double>*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
      
    if (MEASURE_CONVERGENCE_BY_RESIDUAL)
      keff_eigenvalue_iteration.set_tolerance(TOL_PIT_FM);
    else
      keff_eigenvalue_iteration.set_keff_tol(TOL_PIT_FM);
    
    // Adaptivity loop:
    int as = 1; bool done = false; 
    Hermes::vector<SpaceSharedPtr<double> > ref_spaces;    
    ref_spaces.resize(N_EQUATIONS);
    do 
    {
      Loggable::Static::info("---- Adaptivity step %d:", as);
      
      Loggable::Static::info("Solving on fine meshes.");
      
      for (unsigned int i = 0; i < N_EQUATIONS; i++)
      {
        Mesh::ReferenceMeshCreator ref_mesh_creator(meshes[i]);
        MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
        Space<double>::ReferenceSpaceCreator ref_space_creator(spaces[i], ref_mesh);
        SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
        ref_spaces[i] = ref_space;
      }
      
      // Solve the fine mesh problem.
      report_num_dof("Fine mesh power iteration, NDOF: ", ref_spaces);     
      keff_eigenvalue_iteration.set_spaces(ref_spaces);
      keff_eigenvalue_iteration.solve(power_iterates);
      Solution<double>::vector_to_solutions(keff_eigenvalue_iteration.get_sln_vector(), ref_spaces, power_iterates);
              
      report_num_dof("Projecting fine mesh solutions on coarse meshes, NDOF: ", spaces);
      OGProjection<double> ogProjection;
      ogProjection.project_global(spaces, power_iterates, coarse_solutions);
      
      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION && INTERMEDIATE_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        views.show_solutions(coarse_solutions);
        views.show_orders(spaces);
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }
      
      Adapt<double> adaptivity(spaces);
      
      Loggable::Static::info("Calculating errors.");
      Hermes::vector<double> h1_group_errors;
      double h1_err_est = adaptivity.calc_err_est(coarse_solutions, power_iterates, &h1_group_errors) * 100;
      
      // Time measurement.
      cpu_time.tick();
      double cta = cpu_time.accumulated();
      
      // Report results.
      
      // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
      double keff_err = 1e5*fabs(keff_eigenvalue_iteration.get_keff() - REF_K_EFF)/REF_K_EFF;
      
      report_errors("group err_est_coarse (H1): ", h1_group_errors);
      Loggable::Static::info("total err_est_coarse (H1): %g%%", h1_err_est);
      Loggable::Static::info("k_eff err: %g milli-percent", keff_err);
      
      // Add entry to DOF convergence graph.
      int ndof_coarse = Space<double>::get_num_dofs(spaces);
      graph_dof.add_values(0, ndof_coarse, h1_err_est);
      graph_dof.add_values(1, ndof_coarse, keff_err);
      
      // Add entry to CPU convergence graph.
      graph_cpu.add_values(0, cta, h1_err_est);
      graph_cpu.add_values(1, cta, keff_err);
            
      cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (h1_err_est < ERR_STOP || as == MAX_ADAPT_NUM) 
        done = true;
      else 
      {
        Loggable::Static::info("Adapting the coarse meshes.");
        done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);        
        if (Space<double>::get_num_dofs(spaces) >= NDOF_STOP) 
          done = true;
      }
      
      if (!done) as++;
    }
    while (!done);
  }
  else
  {
    if (HERMES_VISUALIZATION)
    {
      views.show_solutions(power_iterates);
      views.show_orders(spaces);
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", power_iterates);
      views.save_orders_vtk("space", spaces);
    }
    
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double keff_err = 1e5*fabs(keff_eigenvalue_iteration.get_keff() - REF_K_EFF)/REF_K_EFF;
    Loggable::Static::info("K_eff error = %g pcm", keff_err);
  }
    
  cpu_time.tick();
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
  
  if (HERMES_VISUALIZATION)
    Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
  
  // Test the unit source normalization.
  Neutronics::PostProcessor pp(NEUTRONICS_SPN);
  pp.normalize_to_unit_fission_source(&power_iterates, *matprop);
  
  SupportClasses::SourceFilter sf(power_iterates, *matprop);
  Loggable::Static::info("Total fission source by normalized flux: %g.", sf.integrate());
  
  delete matprop;
  
  if (HERMES_VISUALIZATION)
  {
    views.show_all_flux_moments(power_iterates, *matprop);
    Views::View::wait();
  }
  if (VTK_VISUALIZATION)
  {
    views.save_solutions_vtk("flux", "flux", power_iterates);
    views.save_orders_vtk("space", spaces);
  }
  
  return 0;
}
