#define HERMES_REPORT_ALL
#include "problem_data.h"
#include "../../utils.h"
#include "definitions.h"

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 2; 

const unsigned int N_EQUATIONS = N_GROUPS;

const int INIT_REF_NUM[N_EQUATIONS] = {  // Initial uniform mesh refinement for the individual solution components.
  2, 2
};
const int P_INIT[N_EQUATIONS] = {        // Initial polynomial orders for the individual solution components. 
  1, 1
}; 
const int SELECTIVE_REF_NUM[N_EQUATIONS] = {
  1, 2
};
        
const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = -1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;  // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 0.1;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  
                                                  
const bool HERMES_VISUALIZATION = true;  // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = true;       // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool INTERMEDIATE_VISUALIZATION = true; // Set to "true" to display coarse mesh solutions during adaptivity.

// Power iteration control.
double k_eff = 1.0;         // Initial eigenvalue approximation.
double TOL_PIT_CM = 1e-5;   // Tolerance for eigenvalue convergence on the coarse mesh.
double TOL_PIT_FM = 1e-5;   // Tolerance for eigenvalue convergence on the fine mesh.

int ref_fn(Element* e)
{
  if ( e->marker == 1 || 
       e->marker == 7 || 
      e->marker == 11 ||
      e->marker == 20 ||
      e->marker == 21 ||
      e->marker == 28 ||
      e->marker == 29 ||
      e->marker == 34 || 
      e->marker == 35 ||
      e->marker == 39 ||
      e->marker == 40 ||
      e->marker == 42 ||
      e->marker == 43 ||
      e->marker == 44 )
    return 0;
  return -1;
}

int main(int argc, char* argv[])
{  
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_param_value(Hermes::exceptionsPrintCallstack, 1);
  Hermes::Hermes2D::Hermes2DApi.set_param_value(Hermes::Hermes2D::numThreads, 2);

  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
    
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load(mesh_file.c_str(), meshes[0]);
  
  // Convert the mesh so that it has one type of elements (optional). 
  // meshes[0]->convert_quads_to_triangles();
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
    
    meshes[i]->refine_by_criterion(ref_fn, SELECTIVE_REF_NUM[i]);
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  meshes[0]->refine_by_criterion(ref_fn, SELECTIVE_REF_NUM[0]);
    
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);

  // Create pointers to solutions on coarse and fine meshes and from the latest power iteration, respectively.
  Hermes::vector<Solution<double>*> power_iterates;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    power_iterates.push_back(new ConstantSolution<double>(meshes[i], 1.0));   
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space<double> *> spaces_;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces_.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
  
  ConstantableSpacesVector spaces(&spaces_);
  
  MaterialProperties::MaterialPropertyMaps *matprop = new MaterialProperties::MaterialPropertyMaps(N_GROUPS, rm_map);
   
  matprop->set_nuSigma_f(nSf);
  matprop->set_Sigma_f(Sf);
  matprop->set_Sigma_t(St);
  matprop->set_Sigma_s(Ss);
  matprop->set_D(D);
  
  matprop->validate();
  
  std::cout << *matprop;
  
  // Initialize the weak formulation.
  CustomWeakForm wf(*matprop, power_iterates, fission_materials, k_eff, bdy_vacuum);
    
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  report_num_dof("Coarse mesh power iteration, NDOF: ", spaces.get());
  Neutronics::keff_eigenvalue_iteration(power_iterates, &wf, spaces.get_const(), matrix_solver, TOL_PIT_CM, 0, true);
  
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
    Hermes::vector<Solution<double>*> coarse_solutions;
    
    // Initialize all the new solution variables.
    for (unsigned int i = 0; i < N_EQUATIONS; i++) 
      coarse_solutions.push_back(new Solution<double>());
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector<double>*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
    
    // Adaptivity loop:
    int as = 1; bool done = false; std::vector<Mesh*> old_meshes(power_iterates.size());
    do 
    {
      Loggable::Static::info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      Loggable::Static::info("Solving on fine meshes.");
      
      ConstantableSpacesVector fine_spaces(Space<double>::construct_refined_spaces(spaces.get(),1));
            
      // Solve the fine mesh problem.
      report_num_dof("Fine mesh power iteration, NDOF: ", fine_spaces.get());
      Neutronics::keff_eigenvalue_iteration(power_iterates, &wf, fine_spaces.get_const(), matrix_solver, TOL_PIT_FM);
      
      // Delete meshes dynamically created in 'construct_refined_spaces' in previous adaptativity iteration
      // (they are still needed in current 'keff_eigenvalue_iteration', but pointers to them get replaced 
      // in this function by pointers to 'fine_spaces', so we have to keep track of them via 'old_meshes').
      if (as > 1)
        for(unsigned int i = 0; i < power_iterates.size(); i++)
          delete old_meshes[i];
        
      report_num_dof("Projecting fine mesh solutions on coarse meshes, NDOF: ", spaces.get());
      OGProjection<double> ogProjection;
      ogProjection.project_global(spaces.get_const(), power_iterates, coarse_solutions);
      
      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION && INTERMEDIATE_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        views.show_solutions(coarse_solutions);
        views.show_orders(spaces.get());
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }
      
      Adapt<double> adaptivity(spaces.get());
      
      Loggable::Static::info("Calculating errors.");
      Hermes::vector<double> h1_moment_errors;
      double h1_err_est = adaptivity.calc_err_est(coarse_solutions, power_iterates, &h1_moment_errors) * 100;
      
      // Time measurement.
      cpu_time.tick();
      double cta = cpu_time.accumulated();
      
      // Report results.
      
      // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
      double keff_err = 1e5*fabs(wf.get_keff() - REF_K_EFF)/REF_K_EFF;
      
      report_errors("odd moment err_est_coarse (H1): ", h1_moment_errors);
      Loggable::Static::info("total err_est_coarse (H1): %g%%", h1_err_est);
      Loggable::Static::info("k_eff err: %g milli-percent", keff_err);
      
      // Add entry to DOF convergence graph.
      int ndof_coarse = Space<double>::get_num_dofs(spaces.get());
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
        if (Space<double>::get_num_dofs(spaces.get()) >= NDOF_STOP) 
          done = true;
      }
      
      if (!done)
      {
        for(unsigned int i = 0; i < power_iterates.size(); i++)
          old_meshes[i] = const_cast<Mesh*>(power_iterates[i]->get_mesh());
        
        // Increase counter.
        as++;
      }
    }
    while (!done);
  }
  else
  {
    if (HERMES_VISUALIZATION)
    {
      views.show_solutions(power_iterates);
      views.show_orders(spaces.get());
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", power_iterates);
      views.save_orders_vtk("space", spaces.get());
    }
    
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double keff_err = 1e5*fabs(wf.get_keff() - REF_K_EFF)/REF_K_EFF;
    Loggable::Static::info("K_eff error = %g pcm", keff_err);
  }
    
  cpu_time.tick();
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
  
  Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
  
  // Test the unit source normalization.
  Neutronics::PostProcessor pp(NEUTRONICS_DIFFUSION);
  pp.normalize_to_unit_fission_source(&power_iterates, *matprop);
  views.show_solutions(power_iterates);
  
  SupportClasses::SourceFilter sf(power_iterates, *matprop, wf.get_fission_regions());
  Loggable::Static::info("Total fission source by normalized flux: %g.", sf.integrate());
  
  delete matprop;
  
  // Wait for the view to be closed.  
  Views::View::wait();
  
  return 0;
}
