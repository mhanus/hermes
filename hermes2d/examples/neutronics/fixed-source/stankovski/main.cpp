#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include "problem_data.h"
#include <sstream>

using namespace RefinementSelectors;

using namespace WeakFormsNeutronics::Multigroup;
using namespace MaterialProperties::SPN;
using namespace CompleteWeakForms::SPN;

const unsigned int N_GROUPS = 1;  // Monoenergetic (single group) problem.
const unsigned int SPN_ORDER = 3; // SP3 approximation

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

const int INIT_REF_NUM[N_EQUATIONS] = {  // Initial uniform mesh refinement for the individual solution components.
  1, 1//, 0, 0
};
const int P_INIT[N_EQUATIONS] = {        // Initial polynomial orders for the individual solution components. 
  1, 1//, 1, 1
};      
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
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
const double ERR_STOP = 1;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const bool HERMES_VISUALIZATION = true;  // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;     // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = true;       // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const int VIEW_COARSE_SOLUTIONS = 0, VIEW_FINE_SOLUTIONS = 1; // Solution sets for optimized visualization of SPN moments.

void report_num_dof(const std::string& msg, const Hermes::vector< Space* > spaces)
{
  std::stringstream ss;
  
  ss << msg << Space::get_num_dofs(spaces[0]);
  
  for (unsigned int i = 1; i < spaces.size(); i++)
    ss << " + " << Space::get_num_dofs(spaces[i]);
  
  if (spaces.size() > 1)
    ss << " = " << Space::get_num_dofs(spaces);
  
  info(ss.str().c_str());
}

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
   
  MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, rm_map);
  matprop.set_iso_src(src);
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
  matprop.validate();
  
  cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  H2DReader mloader;
  mloader.load(mesh_file.c_str(), meshes[0]);
  
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
  
  SupportClasses::SPN::Views views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);

  // Create pointers to final solutions (the fine mesh solutions in case of STRATEGY >= 0).
  Hermes::vector<Solution*> solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    solutions.push_back(new Solution());
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space *> spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces.push_back(new H1Space(meshes[i], P_INIT[i]));
    
  // Initialize the weak formulation.
  DefaultWeakFormFixedSource wf(matprop, SPN_ORDER);
  
  // Initialize the discrete algebraic representation of the problem.
  DiscreteProblem dp(&wf, spaces);
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  if (STRATEGY >= 0)
  {
    // DOF and CPU convergence graphs
    GnuplotGraph graph_dof_est("Error convergence", "NDOF", "log(error)");
    graph_dof_est.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_dof_est.set_log_x();
    graph_dof_est.set_log_y();
    graph_dof_est.show_legend();
    graph_dof_est.show_grid();
    
    GnuplotGraph graph_cpu_est("Error convergence", "CPU time [s]", "log(error)");
    graph_cpu_est.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_dof_est.set_log_x();
    graph_dof_est.set_log_y();
    graph_cpu_est.show_legend();
    graph_cpu_est.show_grid();
        
    // Create pointers to coarse mesh solutions used for error estimation.
    Hermes::vector<Solution*> coarse_solutions;
    
    // Initialize all the new solution variables.
    for (unsigned int i = 0; i < N_EQUATIONS; i++) 
      coarse_solutions.push_back(new Solution());
      
    // Initialize the refinement selectors.
    H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
    
    // Adaptivity loop:
    int as = 1; bool done = false;
    Hermes::vector<Space *>* fine_spaces;
    do 
    {
      info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      info("Solving on fine meshes.");
      
      fine_spaces = Space::construct_refined_spaces(spaces);
      int ndof_fine = Space::get_num_dofs(*fine_spaces);
      
      DiscreteProblem dp(&wf, *fine_spaces);

      // Initial coefficient vector for the Newton's method.  
      scalar* coeff_vec = new scalar[ndof_fine];
      memset(coeff_vec, 0, ndof_fine * sizeof(scalar));
      
      // Perform Newton's iteration.
      bool jacobian_changed = true;
      bool verbose = true;
      if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed, 1e-8, 100, verbose)) 
        error("Newton's iteration failed.");
      
      // Translate the resulting coefficient vector into solutions.
      Solution::vector_to_solutions(coeff_vec, *fine_spaces, solutions);
      
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting fine-mesh solutions onto coarse meshes.");
      OGProjection::project_global(spaces, solutions, coarse_solutions, matrix_solver);
      
      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        views.show_solutions(coarse_solutions, VIEW_COARSE_SOLUTIONS);
        views.show_orders(spaces);
        cpu_time.tick(HERMES_SKIP);
      }
      
      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt adaptivity(spaces);
      
      // Calculate error estimate for each solution component and the total error estimate.
      Hermes::vector<double> err_est_rel;
      double err_est_rel_total = adaptivity.calc_err_est(coarse_solutions, solutions, &err_est_rel) * 100;
      
      // Report results.
      info("total H1 error estimate: %g%%", err_est_rel_total);
      
      cpu_time.tick();
      graph_dof_est.add_values(Space::get_num_dofs(spaces), err_est_rel_total);
      graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
      cpu_time.tick(HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP || as == MAX_ADAPT_NUM) 
        done = true;
      else 
      {
        info("Adapting the coarse meshes.");
        done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);        
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) 
          done = true;
      }
      
      // Clean up.
      delete [] coeff_vec;
      
      if (!done)
      {
        for(unsigned int i = 0; i < N_EQUATIONS; i++)
          delete (*fine_spaces)[i]->get_mesh();
        delete fine_spaces;
        
        // Increase counter.
        as++;
      }
    }
    while (done == false);
    
    graph_dof_est.save("conv_dof_est.gp");
    graph_cpu_est.save("conv_cpu_est.gp");
    
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
    {
      // Make the fine-mesh spaces the final spaces for further analyses.
      delete spaces[i]->get_mesh(); // Alternatively "delete meshes[i]".
      delete spaces[i];
      spaces[i] = fine_spaces->at(i);
      
      // Delete the auxiliary coarse-mesh solutions.
      delete coarse_solutions[i];    
    }
    delete fine_spaces; 
  }
  else
  {
    int ndof = Space::get_num_dofs(spaces);
    
    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof];
    memset(coeff_vec, 0, ndof*sizeof(scalar));
    
    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) 
      error("Newton's iteration failed.");
    
    // Translate the resulting coefficient vector into a Solution.
    Solution::vector_to_solutions(coeff_vec, spaces, solutions);
  } 
  
  delete solver;
  delete matrix;
  delete rhs;
  
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  //
  // Analysis of the solution.
  //
  
  // Visualization.
  
  if (HERMES_VISUALIZATION)
  {
    views.show_solutions(solutions, VIEW_FINE_SOLUTIONS);
    views.show_orders(spaces);
  }
  if (VTK_VISUALIZATION)
  {
    views.save_solutions_vtk("flux", "flux", solutions);
    views.save_orders_vtk("space", spaces);
  }
  
  // Integrate absorption rates and scalar fluxes over specified edit regions and compare with results
  // from the collision probabilities code DRAGON (see problem_data.h).
    
  SupportClasses::PostProcessor pp(NEUTRONICS_SPN);
  
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
    info("Absorption rate integrated over %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, absorption_rates[i], 
         fabs(absorption_rates[i] - ref_integrated_absorption_rates[i])/ref_integrated_absorption_rates[i] * 100);
  for (int i = 0; i < 2*n_pins; i++)
    info("Scalar flux integrated over %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, integrated_fluxes[i],
         fabs(integrated_fluxes[i] - ref_integrated_fluxes[i])/ref_integrated_fluxes[i] * 100);
  for (int i = 0; i < 2*n_pins; i++)
    info("Area of %s (DRAGON reg. #%d) = %f (error %g%%)", edit_regions[i].c_str(), i+1, areas[i],
         fabs(areas[i] - ref_regions_areas[i])/ref_regions_areas[i] * 100);
        
  // Wait for the view to be closed.  
  View::wait();
  
  // Final clean up.
  for(unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    delete spaces[i]->get_mesh(); // In case of STRATEGY == -1, this is the same as "delete meshes[i]"; 
                                  // otherwise, it will delete "fine_spaces[i]->get_mesh()".
    delete spaces[i];
    delete solutions[i];
  }
  
  return 0;
}
