#define HERMES_REPORT_ALL

// Define USE_SPN in problem_data.h in order to use the SPN approximation. 
// Otherwise, the diffusion approximation will be used.
//
#include "problem_data.h"
#include "../../utils.h"
#include "error_estimators.h"

#ifdef USE_SPN
  using namespace Neutronics::SPN;
#else // DIFFUSION
  using namespace Neutronics::Diffusion;
#endif

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.

#ifdef USE_SPN
  const unsigned int SPN_ORDER = 3; // SP3 approximation
#else // DIFFUSION
  const unsigned int SPN_ORDER = 1;
#endif

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

// Initial uniform mesh refinement and initial polynomial orders for the individual solution components. 
#ifdef USE_SPN
  const int INIT_REF_NUM[N_EQUATIONS] = {
    1, 1//, 1//, 1
  };
  const int P_INIT[N_EQUATIONS] = {
    2, 2//, 3//, 1
  };

  const int INIT_REF_NUM_REF[N_EQUATIONS] = {
    3, 3//, 3//, 1
  };
  const int P_INIT_REF[N_EQUATIONS] = {
    6, 6//, 6//, 1
  };
#else // DIFFUSION
  const int INIT_REF_NUM[N_EQUATIONS] = {
    1
  };
  const int P_INIT[N_EQUATIONS] = {
    2
  };

  const int INIT_REF_NUM_REF[N_EQUATIONS] = {
    3
  };
  const int P_INIT_REF[N_EQUATIONS] = {
    5
  };
#endif

const int MESH_REGULARITY = -1;
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Visualisation options.
const bool HERMES_VISUALIZATION = true;        // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = true;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = false; // Set to "true" to display solutions on intermediate meshes during adaptivity.

// There are several run cases of this example (mainly for testing purposes). Check the switches on variable 'run'
// throughout the code to see what each of them exactly does.
//
#ifdef USE_SPN
  const int NUM_RUN_CASES = 1;
  const int RUN_CASES[NUM_RUN_CASES] = { 1 };
#else
  const int NUM_RUN_CASES = 1;
  const int RUN_CASES[NUM_RUN_CASES] = { 1 };
#endif

const bool CALCULATE_REFERENCE_SOLUTION = false;  // Set to true to calculate a reference solution on a globally highly refined mesh
                                                  // before evaluating all the adaptivity run cases.

// Setup the convergence graphs to display appropriate data for selected run cases.
void setup_convergence_graph(GnuplotGraph *graph, const std::set<int>& run_cases)
{
  if (run_cases.find(1) != run_cases.end())
  {
    graph->add_row("H1 err. est. 1", "r", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("H1 err. 1", "r", "-", "o");
  }    
  if (run_cases.find(2) != run_cases.end())
  {
    graph->add_row("L2 err. est. 1", "g", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("L2 err. 1", "g", "-", "o");
  }
  if (run_cases.find(3) != run_cases.end())
  {
    graph->add_row("H1 err. est. 2", "b", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("H1 err. 2", "b", "-", "o");
  }
  if (run_cases.find(4) != run_cases.end())
  {
    graph->add_row("L2 err. est. 2", "k", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("L2 err. 2", "k", "-", "o");
  }
  if (run_cases.find(5) != run_cases.end())
  {
    graph->add_row("H1 err. est. 3", "c", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("H1 err. 3", "c", "-", "o");
  }
  if (run_cases.find(6) != run_cases.end())
  {
    graph->add_row("H1 err. est. 4", "m", "-", ".");
    if (CALCULATE_REFERENCE_SOLUTION)
      graph->add_row("H1 err. 4", "m", "-", "o");
  }

  graph->set_log_x();
  graph->set_log_y();
  graph->show_legend();
  graph->show_grid();
}

int main(int argc, char* argv[])
{ 
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);

  double err_stop;        // Stopping criterion for adaptivity (rel. error tolerance between the
                          // fine mesh and coarse mesh solution in percent).
  int ndof_stop;          // Adaptivity process stops when the number of degrees of freedom grows over
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.
    
  int max_adapt_num;      // Adaptivity process stops when the number of adaptation steps grows over this limit.

  if (CALCULATE_REFERENCE_SOLUTION)
  {
    err_stop = 1e-15;
    ndof_stop = 256e3;
    max_adapt_num = 1e3;
  }
  else
  {
    err_stop = 0.33;
    ndof_stop = 60000;
    max_adapt_num = 30;
  }
  
#ifdef USE_SPN
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, rm_map);
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
#else // DIFFUSION
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, rm_map);  
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_s(Ss);
#endif
  matprop.set_iso_src(src);
  matprop.validate();
  
  std::cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group and pseudo flux.
  Hermes::vector<MeshSharedPtr > meshes, ref_meshes, basic_meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    meshes.push_back(new Mesh());
    basic_meshes.push_back(new Mesh());
    if (CALCULATE_REFERENCE_SOLUTION)
      ref_meshes.push_back(new Mesh());
  }
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load(mesh_file.c_str(), meshes[0]);

  if (CALCULATE_REFERENCE_SOLUTION)
    mesh_reader.load(mesh_file.c_str(), ref_meshes[0]);
  
  meshes[0]->refine_element_id(0,2);
  meshes[0]->refine_element_id(1);
  
  if (CALCULATE_REFERENCE_SOLUTION)
  {
    ref_meshes[0]->refine_element_id(0,2);
    ref_meshes[0]->refine_element_id(1);
  }
  
  // Convert the mesh so that it has one type of elements (optional). 
  //meshes[0]->convert_quads_to_triangles();
  //meshes[0]->convert_triangles_to_quads();
  //ref_meshes[0]->convert_quads_to_triangles();
  //ref_meshes[0]->convert_triangles_to_quads();
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
    if (CALCULATE_REFERENCE_SOLUTION)
      ref_meshes[i]->copy(ref_meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();

    if (CALCULATE_REFERENCE_SOLUTION)
      for (int j = 0; j < INIT_REF_NUM_REF[i]; j++) 
        ref_meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();

  if (CALCULATE_REFERENCE_SOLUTION)
    for (int j = 0; j < INIT_REF_NUM_REF[0]; j++) 
      ref_meshes[0]->refine_all_elements();
  
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    basic_meshes[i]->copy(meshes[i]);
  
#ifdef USE_SPN
  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
#else // DIFFUSION
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
#endif
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);

  // Create pointers to final solutions (the fine mesh solutions in case of strategy >= 0).
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions, ref_solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    solutions.push_back(new Solution<double>());
    if (CALCULATE_REFERENCE_SOLUTION)
      ref_solutions.push_back(new Solution<double>());
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<SpaceSharedPtr<double> > spaces, ref_spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    spaces.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
    if (CALCULATE_REFERENCE_SOLUTION)
      ref_spaces.push_back(new H1Space<double>(ref_meshes[i], P_INIT_REF[i]));
  }
  
  // Initialize the weak formulation.
#ifdef USE_SPN
  WeakForms::FixedSourceProblem wf(matprop, SPN_ORDER);
#else // DIFFUSION
  WeakForms::FixedSourceProblem wf(matprop);
#endif

  // Adaptivity runcases.
  std::set<int> run_cases(RUN_CASES, RUN_CASES + NUM_RUN_CASES);

  // DOF and CPU convergence graphs
  GnuplotGraph graph_dof_est("Error convergence", "NDOF", "error");
  GnuplotGraph graph_cpu_est("Error convergence", "CPU time [s]", "error");
  setup_convergence_graph(&graph_dof_est, run_cases);
  setup_convergence_graph(&graph_cpu_est, run_cases);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();

  if (CALCULATE_REFERENCE_SOLUTION)
  {
    report_num_dof("Solving the reference problem for errors comparison, #DOF: ", ref_spaces);

    // Initialize the discrete algebraic representation of the problem.
    DiscreteProblem<double> dp(&wf, ref_spaces);
    
    // Initial coefficient vector for the Newton's method.  
    
    // Perform Newton's iteration.
    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(true);
    
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      ErrorHandling::error_function("Newton's iteration failed.");
    }
    
    // Translate the resulting coefficient vector into instances of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces, ref_solutions);
    
    cpu_time.tick();

    Loggable::Static::info("\t --- time taken: %g s", cpu_time.accumulated());
    //Loggable::Static::info("\t -------- setup: %g s", newton.get_setup_time());
    //Loggable::Static::info("\t -------- assem: %g s", newton.get_assemble_time());
    //Loggable::Static::info("\t -------- solve: %g s", newton.get_solve_time());

    cpu_time.reset();
  }

  // For each run case defined in RUN_CASES.
  unsigned int run_number = 0;  
  for (std::set<int>::const_iterator run = run_cases.begin(); run != run_cases.end(); ++run, ++run_number)
  {
    Loggable::Static::info("_______________________________________________________________________________________________");
    Loggable::Static::info(("____________________________________________"+tostr(*run)+"__________________________________________________").c_str());

    // Norm in which to calculate errors.
    ProjNormType error_norm;
    
    // Adaptive strategy.
    //  strategy = 0 ... Refine elements until sqrt(threshold) times total error is processed. 
    //                   If more elements have similar errors, refine all to keep the mesh symmetric.
    //  strategy = 1 ... Refine all elements whose error is larger than threshold times maximum element error.
    //  strategy = 2 ... refine all elements whose error is larger than threshold.
    //
    double threshold; 
    int strategy;             
      
    bool use_scalar_flux_factor;       // Multiply the estimate for a given component (pseudo-flux) by the factor with which
                                       // it contributes to the total scalar flux.
    unsigned int error_flags;
    bool use_residual_error_estimator; // Set to "true" to use the volumetric (residual) part of the "Kelly" estimator.
    
    switch (*run)
    {
      case 1:
      {
        strategy = 0;
        threshold = 0.6;
        error_norm = HERMES_H1_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS;
        use_scalar_flux_factor = false;
        use_residual_error_estimator = true;
        break;
      }
      case 2:
      {
        strategy = 0;
        threshold = 0.6;
        error_norm = HERMES_L2_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS;
        use_scalar_flux_factor = false;
        use_residual_error_estimator = true;
        break;
      }
      case 3:
      {
        strategy = 0;
        threshold = 0.6;
        error_norm = HERMES_H1_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS;
        use_scalar_flux_factor = true;
        use_residual_error_estimator = true;
        break;
      }
      case 4:
      {
        strategy = 0;
        threshold = 0.6;
        error_norm = HERMES_L2_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS;
        use_scalar_flux_factor = true;
        use_residual_error_estimator = true;
        break;
      }
      case 5:
      {
        strategy = 1;
        threshold = 0.3;
        error_norm = HERMES_H1_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_REL;
        use_scalar_flux_factor = false;
        use_residual_error_estimator = true;
        break;
      }
      case 6:
      {
        strategy = 1;
        threshold = 0.3;
        error_norm = HERMES_H1_NORM;
        error_flags = HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_REL;
        use_scalar_flux_factor = true;
        use_residual_error_estimator = true;
        break;
      }
    };
    
    cpu_time.tick(TimeMeasurable::HERMES_SKIP);
    
    // Adaptivity loop:
    int as = 1; bool done = false;
    do 
    {
      Loggable::Static::info("---- Adaptivity step %d:", as);
                        
      // Initialize the discrete problem.
      DiscreteProblem<double> dp(&wf, spaces);
      report_num_dof("Solving, #DOF: ", spaces);
      
      // Assemble the Jacobian and perform Newton's iteration. Since the Jacobian may become
      // very large during the h-adaptivity, we destroy the NewtonSolver instance right after
      // it finishes its job to save some memory for error calculations.

      NewtonSolver<double> *newton = new NewtonSolver<double>(&dp);
      newton->set_verbose_output(false);
    
      try
      {
        newton->solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.print_msg();
        ErrorHandling::error_function("Newton's iteration failed.");
      }
      
      // Translate the resulting coefficient vector into instances of Solution.
      Solution<double>::vector_to_solutions(newton->get_sln_vector(), spaces, solutions);
      
      // Clean up.
      delete newton;

      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        if (SHOW_INTERMEDIATE_SOLUTIONS)
          views.show_solutions(solutions);
        if (SHOW_INTERMEDIATE_ORDERS)
          views.show_orders(spaces);
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }
      
      // Calculate element errors.
      Loggable::Static::info("Calculating error estimate.");
      
      Hermes::vector<ProjNormType> norms;
      if (error_norm != HERMES_UNSET_NORM)
        for (unsigned int i = 0; i < N_EQUATIONS; i++)
          norms.push_back(error_norm);
      
      bool ignore_visited_segments = true; 
      KellyTypeAdapt<double> adaptivity(spaces, ignore_visited_segments, 
                                        Hermes::vector<const InterfaceEstimatorScalingFunction*>(), norms);
                                              
  #ifdef USE_SPN
  
      MomentGroupFlattener mg(N_GROUPS); // Creates a single index from the moment-group pair.

      for (unsigned int gto = 0; gto < N_GROUPS; gto++)
        for (unsigned int m = 0; m < N_ODD_MOMENTS; m++)
        {
          adaptivity.add_error_estimator_surf(new InterfaceEstimatorSPN(m, gto, mg, matprop, use_scalar_flux_factor));
          adaptivity.add_error_estimator_surf(new BoundaryEstimatorSPN(m, gto, mg, matprop, use_scalar_flux_factor));
          if(use_residual_error_estimator)
            adaptivity.add_error_estimator_vol(new VolumetricEstimatorSPN(m, gto, mg, matprop, use_scalar_flux_factor));
        }
  #else // DIFFUSION

      for (unsigned int g = 0; g < N_GROUPS; g++)
      {
        adaptivity.add_error_estimator_surf(new InterfaceEstimatorDiffusion(g, matprop));
        adaptivity.add_error_estimator_surf(new BoundaryEstimatorDiffusion(g, matprop));
      }
      
  #endif
        
      double err_scalar_flux = 0.0;
        
      if (CALCULATE_REFERENCE_SOLUTION)
      {
        cpu_time.tick();
        
        Loggable::Static::info("  --- Calculating group-wise errors of scalar flux.");
                    
  #ifdef USE_SPN        
  
        Hermes::vector< MeshFunctionSharedPtr<double> > scalar_fluxes, ref_scalar_fluxes;
        MomentFilter::get_scalar_fluxes_with_derivatives(solutions, &scalar_fluxes, N_GROUPS);
        MomentFilter::get_scalar_fluxes_with_derivatives(ref_solutions, &ref_scalar_fluxes, N_GROUPS);
        
        for (unsigned int g = 0; g < N_GROUPS; g++)
        {
          double group_err = Hermes::sqr(Global<double>::calc_abs_error(scalar_fluxes[g], ref_scalar_fluxes[g], error_norm));
          
          if (error_flags & HERMES_TOTAL_ERROR_REL)
            group_err /= Hermes::sqr(Global<double>::calc_norm(ref_scalar_fluxes[g], error_norm));
          
          err_scalar_flux += group_err;
        }
                
  #else // DIFFUSION
  
        for (unsigned int g = 0; g < N_GROUPS; g++)
        {
          double group_err = Hermes::sqr(Global<double>::calc_abs_error(solutions[g], ref_solutions[g], error_norm));
          if (error_flags & HERMES_TOTAL_ERROR_REL)
            group_err /= Hermes::sqr(Global<double>::calc_norm(ref_solutions[g], error_norm));
          err_scalar_flux += group_err;
        }
        
  #endif
  
        err_scalar_flux = sqrt(err_scalar_flux);
        
        if (error_flags & HERMES_TOTAL_ERROR_REL)
          Loggable::Static::info("  --- total H1 relative error (w.r.t. reference solution): %g%%", err_scalar_flux*100);
        else
          Loggable::Static::info("  --- total H1 relative error (w.r.t. reference solution): %g", err_scalar_flux);
        
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }
      
      Loggable::Static::info("Calculating element error estimates for adaptivity and total error for each pseudo-flux approximation.");
      Hermes::vector<double> err_est_rel;
      
      double err_solutions = adaptivity.calc_err_est(solutions, &err_est_rel, error_flags);

      if (error_flags & HERMES_TOTAL_ERROR_REL)
        Loggable::Static::info("  --- total H1 relative error estimate: %g%%", err_solutions * 100);
      else
        Loggable::Static::info("  --- total H1 absolute error estimate: %g", err_solutions);
      
      cpu_time.tick();

      int ndof = Space<double>::get_num_dofs(spaces);
      if (CALCULATE_REFERENCE_SOLUTION)
      {
        graph_dof_est.add_values(2*run_number, ndof, err_solutions);
        graph_cpu_est.add_values(2*run_number, cpu_time.accumulated(), err_solutions);
        graph_dof_est.add_values(2*run_number+1, ndof, err_scalar_flux);
        graph_cpu_est.add_values(2*run_number+1, cpu_time.accumulated(), err_scalar_flux);
      }
      else
      {
        graph_dof_est.add_values(run_number, ndof, err_solutions);
        graph_cpu_est.add_values(run_number, cpu_time.accumulated(), err_solutions);
      }
      
      cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (err_solutions < err_stop || as == max_adapt_num || ndof >= ndof_stop) 
        done = true;
      else 
      {
        Loggable::Static::info("Adapting the coarse meshes.");
        done = adaptivity.adapt(threshold, strategy, MESH_REGULARITY);        
      }
      
      if (!done) // Increase counter.
        as++;
    }
    while (done == false);
    
    graph_dof_est.save(("conv_dof_est_sp"+tostr(SPN_ORDER)+".gp").c_str());
    graph_cpu_est.save(("conv_cpu_est_sp"+tostr(SPN_ORDER)+".gp").c_str());
          
    cpu_time.tick();
    Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
    cpu_time.reset();
    
    //
    // Analysis of the solution.
    //
    
    // Visualization.
    
    if (HERMES_VISUALIZATION)
    {
      views.show_solutions(solutions);
      views.show_orders(spaces);
      Views::View::wait();
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", solutions);
      views.save_orders_vtk("space", spaces);
    }
    
    // Integrate absorption rates and scalar fluxes over specified edit regions and compare with results
    // from the collision probabilities code DRAGON (see problem_data.h).

  #ifdef USE_SPN
    PostProcessor pp(NEUTRONICS_SPN);
  #else // DIFFUSION
    PostProcessor pp(NEUTRONICS_DIFFUSION);
  #endif
    
    Hermes::vector<double> absorption_rates, integrated_fluxes, areas;
    pp.get_integrated_reaction_rates(ABSORPTION, solutions, &absorption_rates, matprop, edit_regions);
    pp.get_integrated_scalar_fluxes(solutions, &integrated_fluxes, N_GROUPS, edit_regions);
    pp.get_areas(spaces[0]->get_mesh(), edit_regions, &areas); // Areas of the edit regions.
    
    // Multiply integral results by the number of times each region appears in the assembly.
    for (int i = 0; i < 2*n_pins; i++)
    {
      absorption_rates[i] *= 8;
      integrated_fluxes[i] *= 8;
      areas[i] *= 8;
    }
    
    for (int i = 0; i < 2*n_pins; i++)
      Loggable::Static::info("Absorption rate integrated over %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, absorption_rates[i], 
          fabs(absorption_rates[i] - ref_integrated_absorption_rates[i])/ref_integrated_absorption_rates[i] * 100);
    for (int i = 0; i < 2*n_pins; i++)
      Loggable::Static::info("Scalar flux integrated over %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, integrated_fluxes[i],
          fabs(integrated_fluxes[i] - ref_integrated_fluxes[i])/ref_integrated_fluxes[i] * 100);
    for (int i = 0; i < 2*n_pins; i++)
      Loggable::Static::info("Area of %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, areas[i],
          fabs(areas[i] - ref_regions_areas[i])/ref_regions_areas[i] * 100);
          
    for(unsigned int i = 0; i < N_EQUATIONS; i++)
    {
        meshes[i]->copy(basic_meshes[i]);
        
        H1SpaceSharedPtr<double> sp = static_cast< H1SpaceSharedPtr<double> >(spaces[i]->duplicate(meshes[i]));
        delete spaces[i];
        spaces[i] = sp;
    }
  }        
    
  // Final clean up.
  for(unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    delete solutions[i];
    delete meshes[i];
    delete spaces[i];
    delete basic_meshes[i];
    
    if (CALCULATE_REFERENCE_SOLUTION)
    {
      delete ref_meshes[i];
      delete ref_spaces[i];
      delete ref_solutions[i];
    }
  }
  
  return 0;
}
