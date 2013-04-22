#define HERMES_REPORT_ALL

// Define USE_SPN in problem_data.h in order to use the SPN approximation. 
// Otherwise, the diffusion approximation will be used.
//
#include "problem_data.h"
#include "../../utils.h"
#include "error_estimators.h"

// Include namespace with variables for controlling the projection-based adaptivity.
using namespace RefinementSelectors;

#ifdef USE_SPN 
  using namespace Neutronics::SPN;
#else // DIFFUSION
  using namespace Neutronics::Diffusion;
#endif

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.

#ifdef USE_SPN
  const unsigned int SPN_ORDER = 3; // Currently implemented maximum is 9.
#else // DIFFUSION
  const unsigned int SPN_ORDER = 1;
#endif

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

// Initial uniform mesh refinement and initial polynomial orders for the individual solution components. 
#ifdef USE_SPN
  const int INIT_REF_NUM[N_EQUATIONS] = {
    1, 1//, 1//, 1//, 1
  };
  const int P_INIT[N_EQUATIONS] = {
    1, 1//, 1//, 1//, 1
  };
#else // DIFFUSION
  const int INIT_REF_NUM[N_EQUATIONS] = {
    1
  };
  const int P_INIT[N_EQUATIONS] = {
    1
  };
#endif

const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = -1... do not use adaptivity.
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric. Absolute element errors should be
                                         //   used in Adapt::calc_err_est for optimal multi-mesh adaptivity.
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
const double ERR_STOP = 0.33;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent). Setting a very small value
                                         // effectively disables this convergence test in favor of the following ones.
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit.
const int MAX_ADAPT_NUM = 60;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
                                         
// Solvers used for solving the discrete problems. Possibilities (depends on how you configured the hermes library): 
//  SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                         
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

// Visualisation options.
const bool HERMES_VISUALIZATION = true;        // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = true;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = true;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

// There are several run cases of this example (mainly for testing purposes). Check the switches on variable 'run'
// throughout the code to see what each of them exactly does. A short description follows.
//
// runcase    description of error/norm setting
// ----------------------------------------------------------------------------------------------------------------------------
//    0       E/N: all terms in the scalar flux expansion in terms of pseudo-fluxes, including the expansion coeffs.
//            This is the only viable option in case of diffusion calculation (i.e. undefined USE_SPN), where default Hermes
//            setting (standard H1 norm) is used for both error evaluation and normalization.
//    1       E: all terms with coeffs., N: standard H1 norm covering all terms (i.e. without the expansion coeffs.)
//    2       E/N: only unmixed terms (i.e. F1F1, not F1F3), including the appropriate expansion coeffs. (i.e., squared)
//    3       E: unmixed terms with coeffs., N: H1 norm covering only the unmixed terms, without coeffs.
//    4       E/N: unmixed H1, without coeffs. (default in Hermes, i.e. if the set_error/norm_form wasn't called at all)
//    5       E/N: mixed H1, without coeffs. (second choice in Hermes, typically used when solving a system with interactions
//            between solution components)
//
#ifdef USE_SPN
  const int NUM_RUN_CASES = 1;
  const int RUN_CASES[NUM_RUN_CASES] = { 0 };
#else
  const int NUM_RUN_CASES = 1;
  const int RUN_CASES[NUM_RUN_CASES] = { 0 };
#endif

// Setup the convergence graphs to display appropriate data for selected run cases.
void setup_convergence_graph(GnuplotGraph *graph, const std::set<int>& run_cases)
{
  static const char* colors[] = { "r", "g", "b", "c", "m", "k" };
  
  for (std::set<int>::const_iterator run = run_cases.begin(); run != run_cases.end(); ++run)
  {
#ifdef USE_SPN
    graph->add_row(("pseudo-fluxes "+itos(*run)).c_str(), colors[*run], "-", "o");
#endif
    graph->add_row(("scalar fluxes "+itos(*run)).c_str(), colors[*run], "-", "+");
    graph->add_row(("absorb. rates "+itos(*run)).c_str(), colors[*run], "-", "p");
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

  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load material data.
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
  
  // Print material data.
  std::cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group and pseudo-flux.
  Hermes::vector<MeshSharedPtr > meshes, basic_meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    meshes.push_back(new Mesh());
    basic_meshes.push_back(new Mesh()); // Before each of the repeated runs of the main calculation 
                                        // loop, meshes will be reset to 'basic_meshes'.
  }
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load(mesh_file.c_str(), meshes[0]);
  meshes[0]->refine_element_id(0,2);
  meshes[0]->refine_element_id(1);
   
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
  
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    basic_meshes[i]->copy(meshes[i]);
  
#ifdef USE_SPN
  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
#else // DIFFUSION
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
#endif
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);

  // Create pointers to final solutions (the fine mesh solutions in case of STRATEGY >= 0).
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    solutions.push_back(new Solution<double>());
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<SpaceSharedPtr<double> > spaces, basic_spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    spaces.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
    basic_spaces.push_back(new H1Space<double>(basic_meshes[i], P_INIT[i]));
  }
  
  // Initialize the weak formulation.
#ifdef USE_SPN
  WeakForms::FixedSourceProblem wf(matprop, SPN_ORDER);
#else // DIFFUSION
  WeakForms::FixedSourceProblem wf(matprop);
#endif
  
  if (STRATEGY >= 0)  // Use adaptivity.
  {
    std::set<int> run_cases(RUN_CASES, RUN_CASES + NUM_RUN_CASES);
    
    // DOF and CPU convergence graphs.
    GnuplotGraph graph_dof_est("Convergence of H1 estimates of rel. error", "NDOF", "error [%]");    
    GnuplotGraph graph_cpu_est("Convergence of H1 estimates of rel. error", "CPU time [s]", "error [%]");
    setup_convergence_graph(&graph_dof_est, run_cases);
    setup_convergence_graph(&graph_cpu_est, run_cases);
        
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
    
    // For each run case defined in RUN_CASES.
    unsigned int run_number = 0;
    for (std::set<int>::const_iterator run = run_cases.begin(); run != run_cases.end(); ++run, ++run_number)
    {
      Loggable::Static::info("_______________________________________________________________________________");
      Loggable::Static::info(("____________________________________"+itos(*run)+"__________________________________________").c_str());

      ProjNormType error_norm = HERMES_H1_NORM;

      // Adaptivity loop:
      int as = 1; bool done = false;
      do 
      {
        Loggable::Static::info("---- Adaptivity step %d:", as);
        
        // Initialize the fine mesh problem.
        Hermes::vector<SpaceSharedPtr<double> > fine_spaces = Space<double>::construct_refined_spaces(spaces);
        int ndof_fine = Space<double>::get_num_dofs(fine_spaces);
        
        report_num_dof("Solving on fine meshes, #DOF: ", fine_spaces);
      
        DiscreteProblem<double> dp(&wf, fine_spaces);

        // Perform Newton's iteration on reference mesh.
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
        Solution<double>::vector_to_solutions(newton->get_sln_vector(), fine_spaces, solutions);

        // Solution of discrete problem completed, delete the solver object.
        delete newton;
        
        // Project the fine mesh solution onto the coarse mesh.
        report_num_dof("Projecting fine-mesh solutions onto coarse meshes, #DOF: ", spaces);
        OGProjection<double> ogProjection;
        ogProjection.project_global(spaces, solutions, coarse_solutions);

        // View the coarse-mesh solutions and polynomial orders.
        if (HERMES_VISUALIZATION)
        {
          cpu_time.tick();
          Loggable::Static::info("Visualizing.");
          if (SHOW_INTERMEDIATE_SOLUTIONS)
            views.show_solutions(coarse_solutions);
          if (SHOW_INTERMEDIATE_ORDERS)
            views.show_orders(spaces);
          cpu_time.tick(TimeMeasurable::HERMES_SKIP);
        }
        
        // Calculate element errors.
        Loggable::Static::info("Calculating error estimate."); 
        Adapt<double> adaptivity(spaces);
                
  #ifdef USE_SPN            

        cpu_time.tick();
        TimeMeasurable cpu_time2;
        cpu_time2.tick();
        
        Loggable::Static::info("  --- Calculating total relative error of scalar flux approximation.");
        Hermes::vector< MeshFunctionSharedPtr<double> > coarse_scalar_fluxes, fine_scalar_fluxes;
        // If error_norm == HERMES_L2_NORM, MomentFilter::get_scalar_fluxes can be used instead (derivatives will not be needed).
        MomentFilter::get_scalar_fluxes_with_derivatives(coarse_solutions, &coarse_scalar_fluxes, N_GROUPS);
        MomentFilter::get_scalar_fluxes_with_derivatives(solutions, &fine_scalar_fluxes, N_GROUPS);
        
        // Calculate relative error estimate of the scalar flux approximation (linear comb. of the actual solutions).
        MomentGroupFlattener mg(N_GROUPS);  // Creates a single index from the moment-group pair.
        double scalar_flux_err_est_rel = 0.0;
        for (unsigned int g = 0; g < N_GROUPS; g++)
        {
          double group_err_est = Hermes::sqr(Global<double>::calc_abs_error(coarse_scalar_fluxes[g], fine_scalar_fluxes[g], error_norm));
          double group_norm = Hermes::sqr(Global<double>::calc_norm(fine_scalar_fluxes[g], error_norm));
          scalar_flux_err_est_rel += group_err_est/group_norm;
        }  
        
        scalar_flux_err_est_rel = sqrt(scalar_flux_err_est_rel) * 100;
        
        // Calculate relative error estimate of the integrated absorption rate (summed over all groups).
        PostProcessor pp(NEUTRONICS_SPN);
        double abs_rate_coarse = pp.get_integrated_reaction_rates(ABSORPTION, coarse_solutions, matprop, edit_regions);
        double abs_rate_fine = pp.get_integrated_reaction_rates(ABSORPTION, solutions, matprop, edit_regions);
        double absorb_rate_err_est_rel = fabs((abs_rate_coarse - abs_rate_fine)/abs_rate_fine) * 100;
                
        cpu_time2.tick();
        Loggable::Static::info("  ------- Time taken (not included in the total running time): %g s", cpu_time2.last());
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
        
        // Calculate error estimate for each solution component and the total error estimate.
        Loggable::Static::info("  --- Calculating total relative error of pseudo-fluxes approximation.");
        
        // Set the combination of error estimation form / error normalization form corresponding to the current run case.
        for (unsigned int g = 0; g < N_GROUPS; g++)
        {
          for (unsigned int mrow = 0; mrow < N_ODD_MOMENTS; mrow++)
          {
            for (unsigned int mcol = 0; mcol < N_ODD_MOMENTS; mcol++)
            {
              if (*run == 0 || *run == 1)
              {
                adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, error_norm));
                if (*run == 1)
                  adaptivity.set_norm_form(mg.pos(mrow,g), mg.pos(mcol,g), new Adapt<double>::MatrixFormVolError(mg.pos(mrow,g), mg.pos(mcol,g), error_norm));
                else
                  adaptivity.set_norm_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, error_norm));
              }
              else if (*run == 2 || *run == 3)
              {
                if (mcol == mrow)
                {
                  adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, error_norm));
                  if (*run == 3)
                    adaptivity.set_norm_form(mg.pos(mrow,g), mg.pos(mcol,g), new Adapt<double>::MatrixFormVolError(mg.pos(mrow,g), mg.pos(mcol,g), error_norm));
                  else
                    adaptivity.set_norm_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, error_norm));
                }
              }
              else if (*run == 5)
              {
                adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new Adapt<double>::MatrixFormVolError(mg.pos(mrow,g), mg.pos(mcol,g), error_norm));
                adaptivity.set_norm_form(mg.pos(mrow,g), mg.pos(mcol,g), new Adapt<double>::MatrixFormVolError(mg.pos(mrow,g), mg.pos(mcol,g), error_norm));
              }
            }
          }
        }
        
        // Calculate the element and total error estimates and make them available for mesh adaptation.
        double pseudo_flux_err_est_rel = adaptivity.calc_err_est(coarse_solutions, solutions, NULL, true,
                                                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100; 
        
  #else // DIFFUSION

        if (*run != 0)
          ErrorHandling::error_function("Only one run case possible for the diffusion model.");
        
        double scalar_flux_err_est_rel = adaptivity.calc_err_est(coarse_solutions, solutions, NULL, true, 
                                                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;
        
        PostProcessor pp(NEUTRONICS_DIFFUSION);
        double abs_rate_coarse = pp.get_integrated_reaction_rates(ABSORPTION, coarse_solutions, matprop, edit_regions);
        double abs_rate_fine = pp.get_integrated_reaction_rates(ABSORPTION, solutions, matprop, edit_regions);
        double absorb_rate_err_est_rel = fabs((abs_rate_coarse - abs_rate_fine)/abs_rate_fine) * 100;

  #endif
        
        // Report results.
  #ifdef USE_SPN
        Loggable::Static::info("  --- pseudo-fluxes rel. error estimate: %g%%", pseudo_flux_err_est_rel);
  #endif
        Loggable::Static::info("  --- scalar flux rel. error estimate: %g%%", scalar_flux_err_est_rel);
        Loggable::Static::info("  --- integrated absorption rate rel. error estimate: %g%%", absorb_rate_err_est_rel);
        
        // Add the results to convergence graphs.
        cpu_time.tick();
  #ifdef USE_SPN   
        graph_dof_est.add_values(3*run_number, ndof_fine, pseudo_flux_err_est_rel);
        graph_cpu_est.add_values(3*run_number, cpu_time.accumulated(), pseudo_flux_err_est_rel);
        graph_dof_est.add_values(3*run_number+1, ndof_fine, scalar_flux_err_est_rel);
        graph_cpu_est.add_values(3*run_number+1, cpu_time.accumulated(), scalar_flux_err_est_rel);
        graph_dof_est.add_values(3*run_number+2, ndof_fine, absorb_rate_err_est_rel);
        graph_cpu_est.add_values(3*run_number+2, cpu_time.accumulated(), absorb_rate_err_est_rel);
  #else // DIFFUSION
        graph_dof_est.add_values(2*run_number, ndof_fine, scalar_flux_err_est_rel);
        graph_cpu_est.add_values(2*run_number, cpu_time.accumulated(), scalar_flux_err_est_rel);
        graph_dof_est.add_values(2*run_number+1, ndof_fine, absorb_rate_err_est_rel);
        graph_cpu_est.add_values(2*run_number+1, cpu_time.accumulated(), absorb_rate_err_est_rel);
  #endif
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
        
        // If err_est is too large, adapt the mesh.
        if (scalar_flux_err_est_rel < ERR_STOP || as == MAX_ADAPT_NUM || ndof_fine >= NDOF_STOP) 
          done = true;
        else 
        {
          Loggable::Static::info("Adapting the coarse meshes.");
          done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);
        }
        
        if (!done)
        {
          for(unsigned int i = 0; i < N_EQUATIONS; i++)
            delete fine_spaces[i]->get_mesh();
          delete &fine_spaces;
          
          // Increase the adaptivity step counter.
          as++;
        }
        else
        {
          cpu_time.tick();
          Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
          cpu_time.reset();
          
          // Visualization.
          
          if (HERMES_VISUALIZATION)
          {
            views.show_solutions(solutions);
            views.show_orders(spaces);
            
            // Wait for the view to be closed.  
            Views::View::wait();
          }
          if (VTK_VISUALIZATION)
          {
            views.save_solutions_vtk("flux", "flux", solutions);
            views.save_orders_vtk("space", spaces);
          }
      
          // Make the fine-mesh spaces the final spaces for further analyses.
          for (unsigned int i = 0; i < N_EQUATIONS; i++)
            delete spaces[i]->get_mesh(); // Alternatively "delete meshes[i]".   
          spaces = fine_spaces; 
        }
      }
      while (done == false);
      
      // Save the convergence graphs.
      graph_dof_est.save(("conv_dof_est_sp"+itos(SPN_ORDER)+".gp").c_str());
      graph_cpu_est.save(("conv_cpu_est_sp"+itos(SPN_ORDER)+".gp").c_str());
      
      
      //
      // Analysis of the solution.
      //
      
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
        
      // Reset the meshes and spaces for the next run case.
      for(unsigned int i = 0; i < N_EQUATIONS; i++)
      {       
        // Destroy the mesh created by construct_refined_spaces and set it to a new copy of basic_mesh.
        delete spaces[i]->get_mesh(); 
        meshes[i] = new Mesh();
        meshes[i]->copy(basic_meshes[i]); 
        
        delete spaces[i]; // Destroy the space created by construct_refined_spaces.
        spaces[i] = basic_spaces[i]->duplicate(meshes[i]); // Create a new space over the fresh copy of basic_mesh, duplicating the basic_space.
      }
      delete &spaces; // Destroy the Hermes::vector<H1SpaceSharedPtr<double> > created by construct_refined_spaces.
      spaces.set(&spaces); // Set a new ConstantableSpacesVector to point to the newly created spaces.
            
      // Advance to the next run case ...
    }
  }
  else  // Do not use adaptivity.
  { 
    report_num_dof("Solving - #DOF: ", spaces);
    
    cpu_time.tick();
    
    DiscreteProblem<double> dp(&wf, spaces);
  
    // Perform Newton's iteration on reference mesh.
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
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, solutions);
    
    cpu_time.tick();
    Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
    
    //
    // Analysis of the solution (the same as if adaptivity was used).
    //
    
    // Visualization.
    
    if (HERMES_VISUALIZATION)
    {
      views.show_solutions(solutions);
      views.show_orders(spaces);
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

  } 
  
  // Final clean up.
  for(unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    delete spaces[i]->get_mesh(); // In case of STRATEGY == -1, this is the same as "delete meshes[i]"; 
                                        // otherwise, it will delete "fine_spaces[i]->get_mesh()".
    delete spaces[i];
    delete solutions[i];
    
    delete basic_meshes[i];
    delete basic_spaces[i];
  }
  
  return 0;
}
