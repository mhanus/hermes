#define HERMES_REPORT_ALL
#include "problem_data.h"
#include "definitions.h"
#include "../../../utils.h"

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 2;  // Monoenergetic (single group) problem.

const unsigned int N_EQUATIONS = N_GROUPS;

const int INIT_REF_NUM[N_EQUATIONS] = {  // Initial uniform mesh refinement for the individual solution components.
  0,0
};
const int P_INIT[N_EQUATIONS] = {        // Initial polynomial orders for the individual solution components. 
  1,1
}; 
        
const double THRESHOLD = 0.2;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
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
                                                  
                                                  
const bool HERMES_VISUALIZATION = false;  // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = true;       // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool INTERMEDIATE_VISUALIZATION = true; // Set to "true" to display coarse mesh solutions during adaptivity.

// Power iteration control.
double TOL_PIT_CM = 1e-5;   // Tolerance for eigenvalue convergence on the coarse mesh.
double TOL_PIT_FM = 1e-6;   // Tolerance for eigenvalue convergence on the fine mesh.

int main(int argc, char* argv[])
{  
    // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  //Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  MaterialProperties::MaterialPropertyMaps *matprop;
  
  matprop = new MaterialProperties::MaterialPropertyMaps(N_GROUPS, rm_map);
    
  matprop->set_D(D);
  matprop->set_Sigma_f(Sf);
  matprop->set_nu(nu);
  matprop->set_Sigma_a(Sa);
  matprop->set_Sigma_s(Ss);
  matprop->set_fission_materials(fission_materials);
  
  matprop->validate();
  
  std::cout << *matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<MeshSharedPtr > meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load(mesh_file.c_str(), meshes[0]);
  
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
  
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);
  
  // Initialize the weak formulation.
  CustomWeakForm wf(*matprop, HermesMultiArray<std::string>(bdy_vacuum));

  // Create pointers to solutions from the latest power iteration and approximation spaces with default shapeset.
  Hermes::vector<MeshFunctionSharedPtr<double> > power_iterates;  
  Hermes::vector<SpaceSharedPtr<double> > spaces_;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    spaces_.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
    power_iterates.push_back(new Solution<double>());
  }
  
  ConstantableSpacesVector spaces(&spaces_);
    
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  report_num_dof("Coarse mesh power iteration, NDOF: ", spaces);
  
  Neutronics::KeffEigenvalueIteration keff_eigenvalue_iteration(&wf, spaces);
  keff_eigenvalue_iteration.set_picard_tol(TOL_PIT_CM);
  keff_eigenvalue_iteration.set_picard_max_iter(1000);
  //keff_eigenvalue_iteration.set_matrix_E_matrix_dump_format(Hermes::Algebra::DF_HERMES_BIN);
  //keff_eigenvalue_iteration.set_matrix_filename("A");
  //keff_eigenvalue_iteration.set_matrix_number_format("%1.15f");
  //keff_eigenvalue_iteration.output_matrix(1);
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
    
    // Adaptivity loop:
    int as = 1; bool done = false; 
    Hermes::vector<SpaceSharedPtr<double> > ref_spaces_;    
    ref_spaces_.resize(N_EQUATIONS);
    std::vector<MeshSharedPtr> old_meshes(N_EQUATIONS);
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
        ref_spaces_[i] = ref_space;
      }
      
      ConstantableSpacesVector ref_spaces(&ref_spaces_);
      
      // Solve the fine mesh problem.
      report_num_dof("Fine mesh power iteration, NDOF: ", ref_spaces);     
      keff_eigenvalue_iteration.set_spaces(ref_spaces);
      keff_eigenvalue_iteration.solve(power_iterates);
      Solution<double>::vector_to_solutions(keff_eigenvalue_iteration.get_sln_vector(), ref_spaces, power_iterates);
      
      // Delete meshes dynamically created in previous adaptivity iteration (they are needed to project previous power_iterates 
      // to current reference meshes during the above call of keff_eigenvalue_iteration.solve(power_iterates) ).
      if (as > 1)
        for(unsigned int i = 0; i < N_EQUATIONS; i++)
          delete old_meshes[i];
        
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
      double kinf_err = 1e5*fabs(keff_eigenvalue_iteration.get_keff() - REF_K_INF)/REF_K_INF;
      
      report_errors("group err_est_coarse (H1): ", h1_group_errors);
      Loggable::Static::info("total err_est_coarse (H1): %g%%", h1_err_est);
      Loggable::Static::info("k_eff err: %g milli-percent", kinf_err);
      
      // Add entry to DOF convergence graph.
      int ndof_coarse = Space<double>::get_num_dofs(spaces);
      graph_dof.add_values(0, ndof_coarse, h1_err_est);
      graph_dof.add_values(1, ndof_coarse, kinf_err);
      
      // Add entry to CPU convergence graph.
      graph_cpu.add_values(0, cta, h1_err_est);
      graph_cpu.add_values(1, cta, kinf_err);
            
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
      
      if (!done)
      {
        for(unsigned int i = 0; i < N_EQUATIONS; i++)
        {
          old_meshes[i] = ref_spaces_[i]->get_mesh();
          delete ref_spaces_[i];
        }
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
      views.show_orders(spaces);
    }
    if (VTK_VISUALIZATION)
    {
      views.save_solutions_vtk("flux", "flux", power_iterates);
      views.save_orders_vtk("space", spaces);
    }
    
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double kinf_err = 1e5*fabs(keff_eigenvalue_iteration.get_keff() - REF_K_INF)/REF_K_INF;
    Loggable::Static::info("K_eff error = %g pcm", kinf_err);
  }
    
  cpu_time.tick();
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
  
  if (HERMES_VISUALIZATION)
    Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
  
  // Test the unit source normalization.
  Neutronics::PostProcessor pp(NEUTRONICS_DIFFUSION);
  pp.normalize_to_unit_fission_source(&power_iterates, *matprop);
  
  SupportClasses::SourceFilter sf(power_iterates, *matprop);
  Loggable::Static::info("Total fission source by normalized flux: %g.", sf.integrate());
  
  delete matprop;
  
  if (HERMES_VISUALIZATION)
  {
    views.show_solutions(power_iterates);
    Views::View::wait();
  }
  if (VTK_VISUALIZATION)
  {
    views.save_solutions_vtk("flux", "flux", power_iterates);
    views.save_orders_vtk("space", spaces);
  }
  
  return 0;
}
