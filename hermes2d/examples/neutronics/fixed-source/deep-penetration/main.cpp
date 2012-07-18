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
#include <iterator>

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 2;    // Monoenergetic (single group) problem.
#ifdef USE_SPN
  const unsigned int SPN_ORDER = 5; // Currently implemented maximum is 9.
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
    1,  1,     // SP1
    1,  1,     // SP3
    1,  1/*,     // SP5
    5,  5,     // SP7
    5,  5      // SP9
  */};
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    1,  1,     // SP1
    1,  1,     // SP3
    1,  1/*,     // SP5
    2,  2,     // SP7
    2,  2      // SP9
  */};
#else // DIFFUSION
  const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    5,  5
  };
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    2,  2
  };
#endif

const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
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
const bool HERMES_VISUALIZATION = false;         // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = true;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = true;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

const bool SAVE_FLUX_PROFILES = true;

// Problem parameters:

// Isotropic distributed sources:
const MaterialPropertyMap1 src = material_property_map<rank1>
(
  "source", row(0.006546)(0.017701)
)(
  "media", row(0.0)(0.0)
);

// Total cross-sections:
const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "source", row(0.10108)(0.108529)
)(
  "media", row(0.10108)(0.108529)
);

// Isotropic scattering cross-sections:
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "source", 
  page(
    matrix(
      row(0.015923)(0.0)
    )(
      row(0.023434)(0.012502)
    )
  )(
    matrix(
      row(0.008976)(0.0)
    )(
      row(0.009016)(0.003914)
    )
  )
)(
  "media", 
  page(
    matrix(
      row(0.015923)(0.0)
    )(
      row(0.023434)(0.012502)
    )
  )(
    matrix(
      row(0.008976)(0.0)
    )(
      row(0.009016)(0.003914)
    )
  )
);

// Reference solutions - average fluxes by composition:
const Hermes::vector<std::string> edit_regions = HermesMultiArray<std::string>
  ("source")("media");  
const double average_fluxes_dragon[N_GROUPS][2] = {
  { 6.9860904E-02, 1.8559563E-03 },
  { 1.8624046E-01, 4.4639211E-03 }
};
const double areas_dragon[2] = {
  3.90000E+03, 1.47200E+04
};

int main(int argc, char* argv[])
{
  // Time measurement.
  TimeMeasurable cpu_time;
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
                                                   std::set<std::string>(edit_regions.begin(), edit_regions.end()));
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
                                                   std::set<std::string>(edit_regions.begin(), edit_regions.end()));
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
  
  // Get region areas for flux averaging.
#ifdef USE_SPN
  PostProcessor pp(NEUTRONICS_SPN);      
#else // DIFFUSION
  PostProcessor pp(NEUTRONICS_DIFFUSION);
#endif

  Hermes::vector<double> areas;
  pp.get_areas(meshes[0], edit_regions, &areas);
  for (unsigned int i = 0; i < edit_regions.size(); i++)
    Loggable::Static::info("Region \"%s\": A=%gcm2 (A_DRAGON=%gcm2)", edit_regions[i].c_str(), areas[i], areas_dragon[i]);
  
  if (STRATEGY >= 0)
  { 
    // DOF and CPU convergence graphs initialization.
    GnuplotGraph graph_dof("Convergence of rel. errors", "NDOF", "error [%]");
  #ifdef USE_SPN
    graph_dof.add_row("pseudo-fluxes (H1)", "r", "-", "+");
  #else // DIFFUSION
    graph_dof.add_row("fluxes (H1)", "r", "-", "+");
  #endif
    graph_dof.add_row("avg. fluxes w.r.t. DRAGON", "b", "-", "+");
    graph_dof.add_row("avg. fluxes est.", "b", "--", "o");
    
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
    graph_cpu.add_row("l2 error w.r.t. DRAGON", "b", "-", "+");
    graph_cpu.add_row("avg. fluxes est.", "b", "--", "o");
    
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
      Loggable::Static::info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      ConstantableSpacesVector fine_spaces(Space<double>::construct_refined_spaces(spaces.get()));
      int ndof_fine = Space<double>::get_num_dofs(fine_spaces.get());
    
      report_num_dof("Solving on fine meshes, #DOF: ", fine_spaces.get());
    
      DiscreteProblem<double> dp(&wf, fine_spaces.get_const());
    
      // Perform Newton's iteration on reference mesh.
      NewtonSolver<double> *newton = new NewtonSolver<double>(&dp);
      newton->set_verbose_output(false);
    
      try
      {
        newton->solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.printMsg();
        error_function("Newton's iteration failed.");
      }
      
      // Translate the resulting coefficient vector into instances of Solution.
      Solution<double>::vector_to_solutions(newton->get_sln_vector(), fine_spaces.get_const(), solutions);
      
      // Clean up.
      delete newton;
      
      // Project the fine mesh solution onto the coarse mesh.
      report_num_dof("Projecting fine-mesh solutions onto coarse meshes, #DOF: ", spaces.get());
      OGProjection<double>::project_global(spaces.get_const(), solutions, coarse_solutions);

      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        if (SHOW_INTERMEDIATE_SOLUTIONS)
          views.show_solutions(coarse_solutions);
        if (SHOW_INTERMEDIATE_ORDERS)
          views.show_orders(spaces.get());
        cpu_time.tick(Hermes::HERMES_SKIP);
      }   

      cpu_time.tick();
      
      // Calculate element errors.
      Loggable::Static::info("Calculating error estimate.");
      
      Loggable::Static::info("  --- Calculating total relative error of scalar flux approximation.");
      
  #ifdef USE_SPN    
      Hermes::vector< MeshFunction<double>* >* coarse_scalar_fluxes = new Hermes::vector< MeshFunction<double>* >();
      Hermes::vector< MeshFunction<double>* >* fine_scalar_fluxes = new Hermes::vector< MeshFunction<double>* >();
      
      MomentFilter::get_scalar_fluxes(coarse_solutions, coarse_scalar_fluxes, N_GROUPS);
      MomentFilter::get_scalar_fluxes(solutions, fine_scalar_fluxes, N_GROUPS);
  #else // DIFFUSION
      Hermes::vector< Solution<double>* >* coarse_scalar_fluxes = &coarse_solutions;
      Hermes::vector< Solution<double>* >* fine_scalar_fluxes = &solutions;
  #endif

      double avg_flux_err_est_rel = 0.0;
      double avg_flux_err_dragon_rel = 0.0;
        
      for (unsigned int g = 0; g < N_GROUPS; g++)
      {      
        // Calculate relative error (squared) of the region-averaged fluxes (in l2 norm).
        for (unsigned int i = 0; i < 2; i++)
        {
          double avg_group_coarse_flux = pp.integrate(coarse_scalar_fluxes->at(g), edit_regions[i]) / areas[i];
          double avg_group_fine_flux = pp.integrate(fine_scalar_fluxes->at(g), edit_regions[i]) / areas[i];
          
          // Estimate.
          avg_flux_err_est_rel += Hermes::sqr(avg_group_coarse_flux - avg_group_fine_flux) / Hermes::sqr(avg_group_fine_flux);
          
          // "Exact".
          avg_flux_err_dragon_rel += Hermes::sqr(avg_group_fine_flux - average_fluxes_dragon[g][i]) / Hermes::sqr(average_fluxes_dragon[g][i]);
        }
      }  
      
      avg_flux_err_est_rel = sqrt(avg_flux_err_est_rel) * 100;
      avg_flux_err_dragon_rel = sqrt(avg_flux_err_dragon_rel) * 100;    

  #ifdef USE_SPN    
      MomentFilter::clear_scalar_fluxes(coarse_scalar_fluxes);
      MomentFilter::clear_scalar_fluxes(fine_scalar_fluxes);
      delete coarse_scalar_fluxes;
      delete fine_scalar_fluxes;
  #endif

      cpu_time.tick(Hermes::HERMES_SKIP);
      
      // Calculate error estimate for each solution component and the total error estimate.
      Loggable::Static::info("  --- Calculating total relative error of the solution approximation.");
          
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
      Loggable::Static::info("  --- solutions rel. error estimate: %g%%", solution_err_est_rel);
      Loggable::Static::info("  --- l2 norm of estimated rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_est_rel);    
      Loggable::Static::info("  --- l2 norm of \"exact\" rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_dragon_rel);
      
      // Add the results to convergence graphs.
      cpu_time.tick();
      graph_dof.add_values(0, ndof_fine, solution_err_est_rel);
      graph_cpu.add_values(0, cpu_time.accumulated(), solution_err_est_rel);
      graph_dof.add_values(1, ndof_fine, avg_flux_err_dragon_rel);
      graph_cpu.add_values(1, cpu_time.accumulated(), avg_flux_err_dragon_rel);
      graph_dof.add_values(2, ndof_fine, avg_flux_err_est_rel);
      graph_cpu.add_values(2, cpu_time.accumulated(), avg_flux_err_est_rel);
      cpu_time.tick(Hermes::HERMES_SKIP);
      
      // If err_est is too large, adapt the mesh.
      if (solution_err_est_rel < ERR_STOP || as == MAX_ADAPT_NUM || ndof_fine >= NDOF_STOP) 
        done = true;
      else 
      {
        Loggable::Static::info("Adapting the coarse meshes.");
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
        Loggable::Static::info("Total running time: %g s", total_cpu_time);
        cpu_time.reset();
        
        if (HERMES_VISUALIZATION)
        {
          Loggable::Static::info("Visualizing final solutions.");
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
        
    DiscreteProblem<double> dp(&wf, spaces.get_const());
  
    // Perform Newton's iteration on reference mesh.
    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(true);
    
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error_function("Newton's iteration failed.");
    }
  
    // Translate the resulting coefficient vector into instances of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces.get_const(), solutions);
      
    cpu_time.tick();
    total_cpu_time = cpu_time.accumulated();
    Loggable::Static::info("Total running time: %g s", total_cpu_time);
    cpu_time.reset();
        
    if (HERMES_VISUALIZATION)
    {
      Loggable::Static::info("Visualizing solutions.");
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
  std::string file = "integrated_flux-sp"+itos(SPN_ORDER)+".dat";
#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
  std::string file = "integrated_flux-diffusion_trc.dat";
#else // SIMPLE DIFFUSION
  std::string file = "integrated_flux-diffusion.dat";
#endif

  std::cout << std::endl;
  Loggable::Static::info("Calculating errors of integrated fluxes w.r.t. DRAGON and saving to %s", file.c_str());

  FILE *fp = fopen(file.c_str(), "wt");
  fprintf(fp, "Total running time: %g s\n\n", total_cpu_time);
  
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    Hermes::vector<double> group_fluxes;
    pp.get_integrated_group_scalar_fluxes(solutions, &group_fluxes, g, N_GROUPS, edit_regions);
        
    for (int i = 0; i < 2; i++)
    {
      double average_flux = group_fluxes[i] / areas[i];
      Loggable::Static::info("--- Scalar flux in group %d averaged over region \"%s\" = %f (error w.r.t DRAGON = %g%%)", g, edit_regions[i].c_str(), average_flux,
           fabs(average_flux - average_fluxes_dragon[g][i])/average_fluxes_dragon[g][i] * 100);
      fprintf(fp, "Scalar flux in group %d averaged over region \"%s\" = %f (error w.r.t DRAGON = %g%%)\n", g, edit_regions[i].c_str(), average_flux,
              fabs(average_flux - average_fluxes_dragon[g][i])/average_fluxes_dragon[g][i] * 100);
    }
  }
  
  fclose(fp);
  std::cout << std::endl;
  
  if (SAVE_FLUX_PROFILES)
  {
#ifdef USE_SPN
    Hermes::vector< Filter<double>* >* scalar_fluxes =  new Hermes::vector< Filter<double>* >();
    MomentFilter::get_scalar_fluxes(solutions, scalar_fluxes, N_GROUPS);
#else // DIFFUSION
    Hermes::vector< Solution<double>* >* scalar_fluxes;
    scalar_fluxes = &solutions;
#endif

    // Output file names.
#ifdef USE_SPN  
    std::string file1 = "flux_x_1.5625-sp"+itos(SPN_ORDER)+".dat";
    std::string file2 = "flux_65.5_y-sp"+itos(SPN_ORDER)+".dat";
#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
    std::string file1 = "flux_x_1.5625-diffusion_trc.dat";
    std::string file2 = "flux_65.5_y-diffusion_trc.dat";
#else // SIMPLE DIFFUSION
    std::string file1 = "flux_x_1.5625-diffusion.dat";
    std::string file2 = "flux_65.5_y-diffusion.dat";
#endif

    double y = 1.5625;
    
    int npts1 = 16;
    int npts2 = 15;
    int npts3 = 15;
    int npts4 = 40;
    double a1 = 50.;
    double a2 = 65.;
    double a3 = 80.;
    double a4 = 133.;
    double d1 = a1 / npts1;
    double d2 = (a2 - a1) / npts2;
    double d3 = (a3 - a2) / npts3;
    double d4 = (a4 - a3) / npts4;
    
    int npts = npts1 + npts2 + npts3 + npts4;
    double *res = new double [npts*N_GROUPS];
        
    std::ofstream fs1(file1.c_str());
    Loggable::Static::info("Saving the scalar flux profile at y=1.5625cm to %s", file1.c_str());
    
    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::cout << std::endl << "GROUP " << g << std::endl;
      
      double x = d1/2.;
      
      for (int i = 0; i < npts1; i++, x+=d1)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }

      x = x - d1/2. + d2/2.;

      for (int i = npts1; i < npts1 + npts2; i++, x+=d2)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }
      
      x = x - d2/2. + d3/2.;
      
      for (int i = npts1 + npts2; i < npts1 + npts2 + npts3; i++, x+=d3)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }

      x = x - d3/2. + d4/2.;

      for (int i = npts1 + npts2 + npts3; i < npts; i++, x+=d4)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }
      
      std::copy(res+g*npts, res+(g+1)*npts, std::ostream_iterator<double>(fs1, "\n"));
      fs1 << std::endl;
    }
    
    fs1.close();
    delete [] res;
    std::cout << std::endl << std::endl;
    
    double x = 65.5;
    
    npts1 = 16;
    npts2 = 10;
    npts3 = 10;
    npts4 = 50;
    a1 = 50.;
    a2 = 60.;
    a3 = 70.;
    a4 = 140.;
    d1 = a1 / npts1;
    d2 = (a2 - a1) / npts2;
    d3 = (a3 - a2) / npts3;
    d4 = (a4 - a3) / npts4;
    
    npts = npts1 + npts2 + npts3 + npts4;
    res = new double [npts*N_GROUPS];
    
    std::ofstream fs2(file2.c_str());
    Loggable::Static::info("Saving the scalar flux profile at x=65.5cm to %s", file2.c_str());

    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::cout << std::endl << "GROUP " << g << std::endl;
      
      double y = d1/2.;
      
      for (int i = 0; i < npts1; i++, y+=d1)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }

      y = y - d1/2. + d2/2.;

      for (int i = npts1; i < npts1 + npts2; i++, y+=d2)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }
      
      y = y - d2/2. + d3/2.;
      
      for (int i = npts1 + npts2; i < npts1 + npts2 + npts3; i++, y+=d3)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }

      y = y - d3/2. + d4/2.;

      for (int i = npts1 + npts2 + npts3; i < npts; i++, y+=d4)
      {
        //(std::cout << "(" << x << "," << y << ")" << std::endl).flush();
        res[i+g*npts] = scalar_fluxes->at(g)->get_pt_value(x, y);    
      }
      
      std::copy(res+g*npts, res+(g+1)*npts, std::ostream_iterator<double>(fs2, "\n"));
      fs2 << std::endl;
    }
    
    fs2.close();
    delete [] res;

#ifdef USE_SPN
    MomentFilter::clear_scalar_fluxes(scalar_fluxes);
    delete scalar_fluxes;
#endif    
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
