
//
#define HERMES_REPORT_ALL
#include "definitions.h"
#include <iterator>

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 2;    // Monoenergetic (single group) problem.
#ifdef USE_SPN
  const unsigned int SPN_ORDER = 9; // Currently implemented maximum is 9.
#else // DIFFUSION
  const unsigned int SPN_ORDER = 1;
#endif

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = N_MOMENTS/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

// Initial uniform mesh refinement and initial polynomial orders for the individual solution components. 

#ifdef USE_SPN
  const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    6,  6,     // SP1
    6,  6,     // SP3
    6,  6,     // SP5
    6,  6,     // SP7
    6,  6      // SP9
  };
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    2,  2,     // SP1
    2,  2,     // SP3
    2,  2,     // SP5
    2,  2,     // SP7
    2,  2      // SP9
  };
#else // DIFFUSION
  const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    6,  6
  };
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    2,  2
  };
#endif

const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = -1;                  // Adaptive strategy:
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
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
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
const bool HERMES_VISUALIZATION = true;         // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = false;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = true;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

// Problem parameters:

const std::vector<std::string> regions = StdMultiArray<std::string>("Region 1")("Region 2");

// Isotropic distributed sources:
const MaterialPropertyMap1 src = material_property_map<rank1>
(
  "Region 1", row(1.0)(1.0)
)(
  "Region 2", row(0.0)(0.0)
);

// Total cross-sections:
const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "Region 1", row(1.0)(1.0)
)(
  "Region 2", row(1.0)(1.0)
);

// Isotropic scattering cross-sections:
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "Region 1", 
  page(
    matrix(
      row(0.5)(0.0)
    )(
      row(0.5)(0.5)
    )
  )(
    matrix(
      row(0.3)(0.0)
    )(
      row(0.3)(0.3)
    )
  )(
    matrix(
      row(0.2)(0.0)
    )(
      row(0.2)(0.2)
    )
  )(
    matrix(
      row(3./35.)(0.0)
    )(
      row(3./35.)(3./35.)
    )
  )
)(
  "Region 2", 
  page(
    matrix(
      row(0.5)(0.0)
    )(
      row(0.5)(0.5)
    )
  )(
    matrix(
      row(0.3)(0.0)
    )(
      row(0.3)(0.3)
    )
  )(
    matrix(
      row(0.2)(0.0)
    )(
      row(0.2)(0.2)
    )
  )(
    matrix(
      row(3./35.)(0.0)
    )(
      row(3./35.)(3./35.)
    )
  )
);

// Reference solutions - scalar flux values at specified points of interest:
const std::vector<double> points_of_interest = StdMultiArray<double>
  (0.125)(0.375)(1.875)(2.125)(2.375)(2.625)(18.875)(19.125)(19.375)(19.625)(19.875);
const std::vector<rank1> ref_solutions = StdMultiArray<rank1>(
  row(1.75949)(1.74706)(1.16971)(0.78535)(0.55585)(0.41958)(7.77495e-7)(6.35673e-7)(5.18211e-7)(4.19663e-7)(3.33956e-7)
)(
  row(3.20460)(3.17762)(2.16312)(1.64672)(1.29081)(1.04638)(7.21376e-6)(5.95205e-6)(4.89303e-6)(3.99218e-6)(3.19123e-6)
);
  
int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::TimePeriod cpu_time;
  double total_cpu_time;
  cpu_time.tick();

  // Load material data (last argument specifies a list of material, which is required
  // for automatic fill-in of missing standard data (fission properties in this case)).
  
#ifndef USE_SPN
  // Extract the isotropic part from the scattering matrix.
  
  MaterialPropertyMap2 Ss0;
  MaterialPropertyMap3::const_iterator it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss0[it->first] = it->second[0];
#endif

#ifdef USE_SPN  

  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, 
                                                   std::set<std::string>(regions.begin(), regions.end()));
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);

#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)  

  MaterialPropertyMap2 Ss1;
  it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss1[it->first] = it->second[1];
  
  MaterialProperties::TransportCorrectedMaterialPropertyMaps matprop(N_GROUPS, Ss1);
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_s(Ss0);
  
#else // SIMPLE DIFFUSION
  
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, 
                                                   std::set<std::string>(regions.begin(), regions.end()));
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_s(Ss0);

#endif
    
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group and pseudo-flux.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load("domain.mesh", meshes[0]);
  
  int refinement_type = 2;            // Split elements vertically.
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements(refinement_type);
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements(refinement_type);

#ifdef USE_SPN
  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
#else // DIFFUSION
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
#endif

  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);
  
  // Create pointers to the coarse and fine mesh solutions.
  Hermes::vector<Solution<double>*> coarse_solutions, solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    coarse_solutions.push_back(new Solution<double>());
    solutions.push_back(new Solution<double>());
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space<double> *> spaces_;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces_.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
  
  ConstantableSpacesVector spaces(&spaces_);
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, SPN_ORDER);

  if (STRATEGY >= 0)
  {
    // DOF and CPU convergence graphs initialization.
    GnuplotGraph graph_dof("Convergence of rel. errors", "NDOF", "error [%]");
  #ifdef USE_SPN
    graph_dof.add_row("pseudo-fluxes (H1)", "r", "-", "+");
  #else // DIFFUSION
    graph_dof.add_row("fluxes (H1)", "r", "-", "+");
  #endif
    graph_dof.add_row("l2 error w.r.t. ANISN", "b", "-", "+");
    
    graph_dof.set_log_x();
    graph_dof.set_log_y();
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Convergence of rel. errors", "CPU time [s]", "error [%]");
  #ifdef USE_SPN  
    graph_cpu.add_row("pseudo-fluxes (H1)", "r", "-", "+");
  #else // DIFFUSION
    graph_cpu.add_row("fluxes (H1)", "r", "-", "+");
  #endif
    graph_cpu.add_row("l2 error w.r.t. ANISN", "b", "-", "+");
    
    graph_cpu.set_log_x();
    graph_cpu.set_log_y();
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector<double>*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
    
    // Adaptivity loop.
    int as = 1; bool done = false;
    do 
    {
      info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      int order_increase = 1;          // FIXME: This should be increase in the x-direction only.
      int refinement_type = 0;
      ConstantableSpacesVector fine_spaces(Space<double>::construct_refined_spaces(spaces.get(), order_increase, refinement_type));
      int ndof_fine = Space<double>::get_num_dofs(fine_spaces.get());
          
      report_num_dof("Solving on fine meshes, #DOF: ", fine_spaces.get());
    
      DiscreteProblem<double> dp(&wf, fine_spaces.get_const());
    
      // Perform Newton's iteration on reference mesh.
      NewtonSolver<double> *newton = new NewtonSolver<double>(&dp, matrix_solver);
      newton->set_verbose_output(false);
    
      try
      {
        newton->solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.printMsg();
        error("Newton's iteration failed.");
      }
      
      // Translate the resulting coefficient vector into instances of Solution.
      Solution<double>::vector_to_solutions(newton->get_sln_vector(), fine_spaces.get_const(), solutions);
      
      // Clean up.
      delete newton;
      
      // Project the fine mesh solution onto the coarse mesh.
      report_num_dof("Projecting fine-mesh solutions onto coarse meshes, #DOF: ", spaces.get());
      OGProjection<double>::project_global(spaces.get_const(), solutions, coarse_solutions, matrix_solver);

      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        info("Visualizing.");
        if (SHOW_INTERMEDIATE_SOLUTIONS)
          views.show_solutions(coarse_solutions);
        if (SHOW_INTERMEDIATE_ORDERS)
          views.show_orders(spaces.get());
        cpu_time.tick(Hermes::HERMES_SKIP);
      }   

      cpu_time.tick();
      
      // Calculate element errors.
      info("Calculating error estimate.");
      
      info("  --- Calculating total relative error of scalar flux approximation.");
      
  #ifdef USE_SPN    
      Hermes::vector< MeshFunction<double>* >* scalar_fluxes = new Hermes::vector< MeshFunction<double>* >();    
      MomentFilter::get_scalar_fluxes(solutions, scalar_fluxes, N_GROUPS);
  #else // DIFFUSION
      Hermes::vector< Solution<double>* >* scalar_fluxes = &solutions;
  #endif

      double pointwise_flux_err_ref = 0.0;
      for (unsigned int g = 0; g < N_GROUPS; g++)
        for (unsigned int i = 0; i < points_of_interest.size(); i++)
          pointwise_flux_err_ref += Hermes::sqr((scalar_fluxes->at(g)->get_pt_value(points_of_interest[i], 0) - ref_solutions[g][i]) / ref_solutions[g][i]);
      
      double pointwise_flux_rel_err_ref = sqrt(pointwise_flux_err_ref) * 100;

  #ifdef USE_SPN    
      MomentFilter::clear_scalar_fluxes(scalar_fluxes);
      delete scalar_fluxes;
  #endif

      cpu_time.tick(Hermes::HERMES_SKIP);
      
      // Calculate error estimate for each solution component and the total error estimate.
      info("  --- Calculating total relative error of the solution approximation.");
          
      Adapt<double> adaptivity(spaces.get());  
      
  #ifdef USE_SPN    
      MomentGroupFlattener mg(N_GROUPS);  // Creates a single index from the moment-group pair.
      // Set the error estimation/normalization form.
      for (unsigned int g = 0; g < N_GROUPS; g++)
        for (unsigned int mrow = 0; mrow < N_ODD_MOMENTS; mrow++)
          for (unsigned int mcol = 0; mcol < N_ODD_MOMENTS; mcol++)
            adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, HERMES_H1_NORM));
  #endif            
      
      // Calculate the element and total error estimates and make them available for mesh adaptation.
      double solution_err_est_rel = adaptivity.calc_err_est(coarse_solutions, solutions, NULL, true,
                                                            HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100; 
              
      // Report results.
      info("  --- H1 norm of rel. solution error estimates : %g%%", solution_err_est_rel);
      info("  --- l2 norm of rel. errors in pointwise scalar fluxes: %g%%", pointwise_flux_rel_err_ref);
      
      // Add the results to convergence graphs.
      cpu_time.tick();
      graph_dof.add_values(0, ndof_fine, solution_err_est_rel);
      graph_cpu.add_values(0, cpu_time.accumulated(), solution_err_est_rel);
      graph_dof.add_values(1, ndof_fine, pointwise_flux_rel_err_ref);
      graph_cpu.add_values(1, cpu_time.accumulated(), pointwise_flux_rel_err_ref);
      cpu_time.tick(Hermes::HERMES_SKIP);
      
      // If err_est is too large, adapt the mesh.
      if (solution_err_est_rel < ERR_STOP || as == MAX_ADAPT_NUM || ndof_fine >= NDOF_STOP) 
        done = true;
      else 
      {
        info("Adapting the coarse meshes.");
        done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);
      }
      
      if (!done)
      {
        for(unsigned int i = 0; i < N_EQUATIONS; i++)
          delete fine_spaces.get()[i]->get_mesh();
        delete &fine_spaces.get();
        
        // Increase the adaptivity step counter.
        as++;
      }
      else
      {
        cpu_time.tick();
        total_cpu_time = cpu_time.accumulated();
        verbose("Total running time: %g s", total_cpu_time);
        cpu_time.reset();
        
        if (HERMES_VISUALIZATION)
        {
          info("Visualizing final solutions.");
          views.inspect_solutions(solutions);
          views.inspect_orders(fine_spaces.get());
        }
        
        for (unsigned int i = 0; i < N_EQUATIONS; i++)
        {
          // Make the fine-mesh spaces the final spaces for further analyses.
          delete spaces.get()[i]->get_mesh(); // Alternatively "delete meshes[i]".
          delete spaces.get()[i];
          spaces.get()[i] = fine_spaces.get()[i]; 

          // Delete the auxiliary coarse-mesh solutions.
          delete coarse_solutions[i];   
        }
        delete &fine_spaces.get();
      }
    }
    while (done == false);
    
    // Save the convergence graphs.
  #ifdef USE_SPN  
    graph_dof.save(("conv_dof_sp"+itos(SPN_ORDER)+".gp").c_str());
    graph_cpu.save(("conv_cpu_sp"+itos(SPN_ORDER)+".gp").c_str());
  #elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
    graph_dof.save("conv_dof_diffusion_trc.gp");
    graph_cpu.save("conv_cpu_diffusion_trc.gp");
  #else // SIMPLE DIFFUSION
    graph_dof.save("conv_dof_diffusion.gp");
    graph_cpu.save("conv_cpu_diffusion.gp");
  #endif
  }
  else
  {
    report_num_dof("Solving - #DOF: ", spaces.get());
    
    cpu_time.tick();
    
    int ndof = Space<double>::get_num_dofs(spaces.get());
    
    DiscreteProblem<double> dp(&wf, spaces.get_const());

    // Perform Newton's iteration on reference mesh.
    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(true);
  
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
    
    // Translate the resulting coefficient vector into instances of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces.get_const(), solutions);
    
    cpu_time.tick();
    total_cpu_time = cpu_time.accumulated();
    verbose("Total running time: %g s", total_cpu_time);
    cpu_time.reset();
        
    if (HERMES_VISUALIZATION)
    {
      info("Visualizing solutions.");
      views.inspect_solutions(solutions);
      views.inspect_orders(spaces.get());
    }
  }
  
  //
  // Analysis of the solution.
  //
  
  // Visualization.

#ifdef USE_SPN
  if (HERMES_VISUALIZATION)
  {
    views.show_all_flux_moments(solutions, matprop);    
    // Wait for the view to be closed.  
    Views::View::wait();
  }
#endif

  if (VTK_VISUALIZATION)
  {
    views.save_solutions_vtk("flux", "flux", solutions);
    views.save_orders_vtk("space", spaces.get());
  }
  
  // Output file names.
#ifdef USE_SPN  
  std::string file = std::string("pointwise_flux-sp")+itos(SPN_ORDER)+std::string(".dat");
  std::string file_err = std::string("pointwise_flux-sp")+itos(SPN_ORDER)+std::string(".err");
#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
  std::string file = "pointwise_flux-diffusion_trc.dat";
  std::string file_err = "pointwise_flux-diffusion_trc.err";
#else // SIMPLE DIFFUSION
  std::string file = "pointwise_flux-diffusion.dat";
  std::string file_err = "pointwise_flux-diffusion.err";
#endif

  std::cout << std::endl;
  info("Calculating errors of pointwise fluxes w.r.t. ANISN and saving to %s (errors only: %s)", file.c_str(), file_err.c_str());
  std::cout << std::endl;

  FILE *fp = fopen(file.c_str(), "wt");
  FILE *fp_err = fopen(file_err.c_str(), "wt");
  fprintf(fp, "Total running time: %g s\n\n", total_cpu_time);

  Filter<double>* test;
#ifdef USE_SPN
    Hermes::vector< Filter<double>* >* scalar_fluxes =  new Hermes::vector< Filter<double>* >();
    MomentFilter::get_scalar_fluxes(solutions, scalar_fluxes, N_GROUPS);
#else // DIFFUSION
    Hermes::vector< Solution<double>* >* scalar_fluxes;
    scalar_fluxes = &solutions;
#endif

  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    info(" GROUP %d ", g);
    std::cout << std::endl;
    info("Scalar flux at");
    
    fprintf(fp, " GROUP %d \n", g);
    fprintf(fp, "Scalar flux at\n");
    fprintf(fp_err, " GROUP %d \n", g);
    
    for (int i = 0; i < points_of_interest.size(); i++)
    {
      double pt_flux = scalar_fluxes->at(g)->get_pt_value(points_of_interest[i], 0);
      info("\t x = %.3f : %g (error w.r.t ANISN = %g%%)", points_of_interest[i], pt_flux, fabs(pt_flux - ref_solutions[g][i])/ref_solutions[g][i] * 100);
      fprintf(fp, "\t x = %.3f : %g (error w.r.t ANISN = %g%%)\n", points_of_interest[i], pt_flux, fabs(pt_flux - ref_solutions[g][i])/ref_solutions[g][i] * 100);
      fprintf(fp_err, "%g\n", fabs(pt_flux - ref_solutions[g][i])/ref_solutions[g][i] * 100);
    }
    
    std::cout << std::endl;
    fprintf(fp, "\n");
    fprintf(fp_err, "\n");
  }
  
  fclose(fp_err);
  fclose(fp);
  
#ifdef USE_SPN
    MomentFilter::clear_scalar_fluxes(scalar_fluxes);
    delete scalar_fluxes;
#endif    

  // Final clean up.
  for(unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    delete spaces.get()[i]->get_mesh();
    delete spaces.get()[i];
    delete solutions[i];
  }

  return 0;
}
