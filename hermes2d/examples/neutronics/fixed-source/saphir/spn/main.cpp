//
//  IAEA EIR-2 benchmark problem, SPN approximation. Note that when SPN_ORDER == 1, this is approximately
//  equal to the diffusion version (the only difference being the correct handling of the void boundary 
//  conditions in this version as opposed to the zero-Dirichlet approximation used in the diffusion version).
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
#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
const unsigned int SPN_ORDER = 3;   // Currently implemented maximum is 9.

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

// Initial uniform mesh refinement and initial polynomial orders for the individual solution components. 
const int INIT_REF_NUM[N_EQUATIONS] = {
  1, 1//, 1, 1//, 1
};
const int P_INIT[N_EQUATIONS] = {
  1, 1//, 1, 1//, 1
};

const double THRESHOLD = 0.66;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
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
const double ERR_STOP = 0.5;             // Stopping criterion for adaptivity (rel. error tolerance between the
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
const bool HERMES_VISUALIZATION = false;         // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = false;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = false;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

// Problem parameters:

// Isotropic distributed sources:
const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "1", 1.00
)(
  "2", 0.00
)(
  "3", 1.00
)(
  "4", 0.00
)(
  "5", 0.00
);

// Total cross-sections:
const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "1", row(0.60)
)(
  "2", row(0.48)
)(
  "3", row(0.70)  
)(
  "4", row(0.65)
)(
  "5", row(0.90)  
);

// Isotropic scattering cross-sections:
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "1", page(matrix(row(0.53)))
)(
  "2", page(matrix(row(0.20)))
)(
  "3", page(matrix(row(0.66)))
)(
  "4", page(matrix(row(0.50)))
)(
  "5", page(matrix(row(0.89)))
);

// Reference solutions - average fluxes by composition:
const Hermes::vector<std::string> edit_regions = HermesMultiArray<std::string>
  ("1")("2")("3")("4")("5");  
const double ref_average_fluxes_sp3[5] = {
  1.1973E+01, 5.3613E-01, 1.9222E+01, 8.2946E-01, 1.5318E+00
};
const double ref_average_fluxes_sp7[5] = {
  1.1970E+01, 5.3690E-01, 1.9218E+01, 8.3028E-01, 1.5314E+00
};
const double ref_average_fluxes_s8[5] = {
  1.1960E+01, 5.3968E-01, 1.9202E+01, 8.3364E-01, 1.5263E+00
};

int main(int argc, char* argv[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  //Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);

  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load material data (last argument specifies a list of material, which is required
  // for automatic fill-in of missing standard data (fission properties in this case)).
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, 
                                                   std::set<std::string>(edit_regions.begin(), edit_regions.end()));
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group and pseudo-flux.
  Hermes::vector<MeshSharedPtr > meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load("../domain.mesh", meshes[0]);
  
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

  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, 510, 510, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);
  
  // Get region areas for flux averaging.
  PostProcessor pp(NEUTRONICS_SPN);      
  Hermes::vector<double> areas;
  pp.get_areas(meshes[0], edit_regions, &areas);
  
  // Create pointers to the coarse and fine mesh solutions.
  Hermes::vector<MeshFunctionSharedPtr<double> > coarse_solutions, solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    coarse_solutions.push_back(new Solution<double>());
    solutions.push_back(new Solution<double>());
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<SpaceSharedPtr<double> > spaces_;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces_.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
  
  ConstantableSpacesVector spaces(&spaces_);
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, SPN_ORDER);
  
  DiscreteProblem<double> dp(&wf, spaces.get_const());
    
  // Create the solver based on Newton's method.
  NewtonSolver<double> *newton = new NewtonSolver<double>(&dp);
  newton->set_verbose_output(false);
    
  if (STRATEGY < 0)
  {
    report_num_dof("Solving on unadapted meshes, #DOF: ", spaces.get());
  
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
    Solution<double>::vector_to_solutions(newton->get_sln_vector(), spaces.get_const(), solutions);
    
    cpu_time.tick();
    Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
    cpu_time.reset();
    
    // Visualization.
    
    if (HERMES_VISUALIZATION)
    {
      views.show_all_flux_moments(solutions, matprop);    
      // Wait for the view to be closed.  
      Views::View::wait();
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", solutions);
      views.save_orders_vtk("space", spaces.get());
    }
  }
  else
  {
    // DOF and CPU convergence graphs initialization.
    GnuplotGraph graph_dof("Convergence of rel. errors", "NDOF", "error [%]");
    graph_dof.add_row("pseudo-fluxes (H1)", "r", "-", "+");
    graph_dof.add_row("scalar-fluxes (H1)", "g", "-", "+");
    graph_dof.add_row("avg. fluxes w.r.t. S8", "b", "-", "+");
    graph_dof.add_row("avg. fluxes est.", "b", "--", "o");
    
    graph_dof.set_log_x();
    graph_dof.set_log_y();
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Convergence of rel. errors", "CPU time [s]", "error [%]");
    graph_cpu.add_row("pseudo-fluxes (H1)", "r", "-", "+");
    graph_cpu.add_row("scalar-fluxes (H1)", "g", "-", "+");
    graph_cpu.add_row("l2 error w.r.t. S8", "b", "-", "+");
    graph_cpu.add_row("avg. fluxes est.", "b", "--", "o");
    
    graph_cpu.set_log_x();
    graph_cpu.set_log_y();
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector<double>*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
    {
      selectors.push_back(&selector);
      selector.set_error_weights( 0.5*H2DRS_DEFAULT_ERR_WEIGHT_H + 0.5*H2DRS_DEFAULT_ERR_WEIGHT_P );
    }
    
    // Adaptivity loop.
    int as = 1; bool done = false;
    Hermes::vector<SpaceSharedPtr<double> > ref_spaces_;
    ref_spaces_.resize(N_EQUATIONS);
    do 
    {
      Loggable::Static::info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      
      for (unsigned int i = 0; i < N_EQUATIONS; i++)
      {
        Mesh::ReferenceMeshCreator ref_mesh_creator(meshes[i]);
        MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
        Space<double>::ReferenceSpaceCreator ref_space_creator(spaces.get_const()[i], ref_mesh);
        SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
        ref_spaces_[i] = ref_space;
      }
    
      ConstantableSpacesVector ref_spaces(&ref_spaces_);
      int ndof_fine = Space<double>::get_num_dofs(ref_spaces.get());
      
      report_num_dof("Solving on fine meshes, #DOF: ", ref_spaces.get());

      // Perform Newton's iteration on reference mesh.
      newton->set_spaces(ref_spaces.get_const());
    
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
      Solution<double>::vector_to_solutions(newton->get_sln_vector(), ref_spaces.get_const(), solutions);
      
      // Project the fine mesh solution onto the coarse mesh.
      report_num_dof("Projecting fine-mesh solutions into coarse meshes, #DOF: ", spaces.get());
      OGProjection<double> ogProjection;
      ogProjection.project_global(spaces.get_const(), solutions, coarse_solutions);

      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        if (SHOW_INTERMEDIATE_SOLUTIONS)
          views.show_solutions(coarse_solutions);
        if (SHOW_INTERMEDIATE_ORDERS)
          views.show_orders(spaces.get());
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }
      
      // Calculate element errors.
      Loggable::Static::info("Calculating error estimate.");    

      cpu_time.tick();
      
      Loggable::Static::info("  --- Calculating total relative error of scalar flux approximation.");
      Hermes::vector< MeshFunctionSharedPtr<double> > coarse_scalar_fluxes, fine_scalar_fluxes;
      // If error_norm == HERMES_L2_NORM, MomentFilter::get_scalar_fluxes can be used instead (derivatives will not be needed).
      MomentFilter::get_scalar_fluxes_with_derivatives(coarse_solutions, &coarse_scalar_fluxes, N_GROUPS);
      MomentFilter::get_scalar_fluxes_with_derivatives(solutions, &fine_scalar_fluxes, N_GROUPS);
      
      double scalar_flux_err_est_rel = 0.0;
      double avg_flux_err_est_rel = 0.0;
      double avg_flux_err_s8_rel = 0.0;
      double scalar_flux_norm = 0.0;
      for (unsigned int g = 0; g < N_GROUPS; g++)
      {
        // Calculate relative error (squared) of the scalar flux approximation (linear comb. of the actual solutions) in specified norm.
        double group_err_est = Hermes::sqr(Global<double>::calc_abs_error(coarse_scalar_fluxes[g], fine_scalar_fluxes[g], HERMES_H1_NORM));
        double group_norm = Hermes::sqr(Global<double>::calc_norm(fine_scalar_fluxes[g], HERMES_H1_NORM));
        scalar_flux_err_est_rel += group_err_est/group_norm;
        scalar_flux_norm += group_norm;
        
        // Calculate relative error (squared) of the region-averaged fluxes (in l2 norm).
        double avg_group_flux_err_est = 0.0;
        double avg_group_flux_err_s8 = 0.0;      
        double avg_group_fine_flux_norm = 0.0;
        double avg_group_ref_flux_norm = 0.0;
        for (unsigned int i = 0; i < 5; i++)
        {
          double avg_group_coarse_flux = pp.integrate(coarse_scalar_fluxes[g], edit_regions[i]) / areas[i];
          double avg_group_fine_flux = pp.integrate(fine_scalar_fluxes[g], edit_regions[i]) / areas[i];
          
          // Estimate.
          avg_group_flux_err_est += Hermes::sqr(avg_group_coarse_flux - avg_group_fine_flux);
          avg_group_fine_flux_norm += Hermes::sqr(avg_group_fine_flux);
          
          // "Exact".
          avg_group_flux_err_s8 += Hermes::sqr(avg_group_fine_flux - ref_average_fluxes_s8[i]);
          avg_group_ref_flux_norm += Hermes::sqr(ref_average_fluxes_s8[i]);
        }
        avg_flux_err_est_rel += avg_group_flux_err_est / avg_group_fine_flux_norm;
        avg_flux_err_s8_rel += avg_group_flux_err_s8 / avg_group_ref_flux_norm;      
      }  
      
      scalar_flux_err_est_rel = sqrt(scalar_flux_err_est_rel) * 100;
      avg_flux_err_est_rel = sqrt(avg_flux_err_est_rel) * 100;
      avg_flux_err_s8_rel = sqrt(avg_flux_err_s8_rel) * 100;
      scalar_flux_norm = sqrt(scalar_flux_norm);
          
      MomentFilter::clear_scalar_fluxes(&coarse_scalar_fluxes);
      MomentFilter::clear_scalar_fluxes(&fine_scalar_fluxes);
      
      cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      
      // Calculate error estimate for each solution component and the total error estimate.
      Loggable::Static::info("  --- Calculating total relative error of pseudo-fluxes approximation.");
          
      Adapt<double> adaptivity(spaces.get());  
         
      MomentGroupFlattener mg(N_GROUPS);  // Creates a single index from the moment-group pair.
      // Set the error estimation/normalization form.
      for (unsigned int g = 0; g < N_GROUPS; g++)
        for (unsigned int mrow = 0; mrow < N_ODD_MOMENTS; mrow++)
          for (unsigned int mcol = 0; mcol < N_ODD_MOMENTS; mcol++)
            adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, HERMES_H1_NORM));
                
      // Calculate the element and total error estimates and make them available for mesh adaptation.
      double pseudo_flux_err_est_rel = adaptivity.calc_err_est(coarse_solutions, solutions, NULL, true,
                                                              HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS) / scalar_flux_norm * 100; 
              
      // Report results.
      Loggable::Static::info("  --- pseudo-fluxes rel. error estimate: %g%%", pseudo_flux_err_est_rel);
      Loggable::Static::info("  --- scalar flux rel. error estimate: %g%%", scalar_flux_err_est_rel);
      Loggable::Static::info("  --- l2 norm of estimated rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_est_rel);    
      Loggable::Static::info("  --- l2 norm of \"exact\" rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_s8_rel);
      
      // Add the results to convergence graphs.
      cpu_time.tick();
      graph_dof.add_values(0, ndof_fine, pseudo_flux_err_est_rel);
      graph_cpu.add_values(0, cpu_time.accumulated(), pseudo_flux_err_est_rel);
      graph_dof.add_values(1, ndof_fine, scalar_flux_err_est_rel);
      graph_cpu.add_values(1, cpu_time.accumulated(), scalar_flux_err_est_rel);
      graph_dof.add_values(2, ndof_fine, avg_flux_err_s8_rel);
      graph_cpu.add_values(2, cpu_time.accumulated(), avg_flux_err_s8_rel);
      graph_dof.add_values(3, ndof_fine, avg_flux_err_est_rel);
      graph_cpu.add_values(3, cpu_time.accumulated(), avg_flux_err_est_rel);
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
        {
          delete ref_spaces_[i]->get_mesh();
          delete ref_spaces_[i];
        }
        
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
          views.show_all_flux_moments(solutions, matprop);    
          // Wait for the view to be closed.  
          Views::View::wait();
        }
        if (VTK_VISUALIZATION)
        {
          views.save_solutions_vtk("flux", "flux", solutions);
          views.save_orders_vtk("space", spaces.get());
        }
    
        for (unsigned int i = 0; i < N_EQUATIONS; i++)
        {
          // Make the fine-mesh spaces the final spaces for further analyses.
          delete spaces.get()[i]->get_mesh(); // Alternatively "delete meshes[i]".
          delete spaces.get()[i];
          spaces.get()[i] = ref_spaces.get()[i]; 

          // Delete the auxiliary coarse-mesh solutions.
          delete coarse_solutions[i];   
        }
      }
    }
    while (done == false);
    
    // Save the convergence graphs.
    graph_dof.save(("conv_dof_sp"+itos(SPN_ORDER)+".gp").c_str());
    graph_cpu.save(("conv_cpu_sp"+itos(SPN_ORDER)+".gp").c_str());
  }
        
  //
  // Analysis of the solution.
  //
    
  Hermes::vector<double> integrated_fluxes;
  pp.get_integrated_scalar_fluxes(solutions, &integrated_fluxes, N_GROUPS, edit_regions);
        
  for (int i = 0; i < 5; i++)
  {
    double average_flux = integrated_fluxes[i] / areas[i];
    Loggable::Static::info("Scalar flux averaged over region %s = %f (error w.r.t S8 = %g%%)", edit_regions[i].c_str(), average_flux,
         fabs(average_flux - ref_average_fluxes_s8[i])/ref_average_fluxes_s8[i] * 100);
    if (SPN_ORDER == 3)
      Loggable::Static::info("\t (solution by the SP3-equivalent A2 method = %f, error w.r.t S8 = %g%%)", ref_average_fluxes_sp3[i],
           fabs(ref_average_fluxes_sp3[i] - ref_average_fluxes_s8[i])/ref_average_fluxes_s8[i] * 100);
    else if (SPN_ORDER == 7)
      Loggable::Static::info("\t (solution by the SP7-equivalent A4 method = %f, error w.r.t S8 = %g%%)", ref_average_fluxes_sp7[i],
           fabs(ref_average_fluxes_sp7[i] - ref_average_fluxes_s8[i])/ref_average_fluxes_s8[i] * 100);
  }

  // Final clean up.
  for(unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    delete spaces.get()[i]->get_mesh();
    delete spaces.get()[i];
    delete solutions[i];
  }

  return 0;
}