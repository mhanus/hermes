#define HERMES_REPORT_ALL

////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop, unsigned int N,
                                const Hermes::vector<Solution*>& iterates, 
                                const Hermes::vector<std::string>& fission_regions,
                                double init_keff, std::string bdy_vacuum )
  : DefaultWeakFormSourceIteration(matprop, N, iterates, fission_regions, init_keff)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new VacuumBoundaryCondition::Residual(m, N, g, bdy_vacuum, matprop));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(m, n, g, bdy_vacuum, matprop));    
    }
}

// Integral over the active core.
double integrate(MeshFunction* sln, const std::vector<std::string>& fission_regions)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  
  double integral = 0.0;
  Element* e;
  Mesh *mesh = sln->get_mesh();
  
  std::set<int> markers;
  std::vector<std::string>::const_iterator it = fission_regions.begin();
  for ( ; it != fission_regions.end(); ++it)
    markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it));
  
  for_all_active_elements(e, mesh)
  {
    if (markers.find(e->marker) != markers.end())
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double result = 0.0;
      h1_integrate_expression(uval[i]);
      integral += result;
    }
  }
  
  return integral;
}

int power_iteration(const Hermes2D& hermes2d, const MaterialPropertyMaps& matprop, 
                    const Hermes::vector<Space *>& spaces, DefaultWeakFormSourceIteration* wf, 
                    const Hermes::vector<Solution *>& solutions, 
                    const std::vector<std::string>& fission_regions,
                    double tol, SparseMatrix *mat, Vector* rhs, Solver *solver)
{
  // Sanity checks.
  if (spaces.size() != solutions.size()) 
    error("Spaces and solutions supplied to power_iteration do not match."); 
   
  // Initialize the discrete problem.
  DiscreteProblem dp(wf, spaces);
  int ndof = Space::get_num_dofs(spaces);
    
  // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
  // updating the eigenvalue. 
  Hermes::vector<Solution*> new_solutions;
  for (int i = 0; i < solutions.size(); i++) 
    new_solutions.push_back(new Solution(solutions[i]->get_mesh()));
  
  // This power iteration will most probably run on a different mesh than the previous one and so will be different
  // the corresponding algebraic system. We will need to factorize it anew (but then, the L and U factors may be 
  // reused until the next adaptation changes the mesh again).
  // TODO: This could be solved more elegantly by defining a function Solver::reinit().
  solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
  
  // Initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  
  // Force the Jacobian assembling in the first iteration.
  bool Jacobian_changed = true;
  
  bool eigen_done = false; int it = 0;
  do 
  {
    memset(coeff_vec, 0.0, ndof*sizeof(scalar));

    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, mat, rhs, Jacobian_changed, 1e-8, 10, true)) 
      error("Newton's iteration failed.");
    
    // The matrix doesn't change within the power iteration loop, so it does not need to be reassembled again and 
    // the first computed LU factorization may be completely reused in following iterations.
    Jacobian_changed = false;
    solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
    
    // Convert coefficients vector into a set of Solution pointers.
    Solution::vector_to_solutions(solver->get_solution(), spaces, new_solutions);

    // Update fission sources.
    SupportClasses::SPN::SourceFilter new_source(new_solutions, matprop, fission_regions);
    SupportClasses::SPN::SourceFilter old_source(solutions, matprop, fission_regions);

    // Compute the eigenvalue for current iteration.
    double k_new = wf->get_keff() * (integrate(&new_source, fission_regions) / 
                                     integrate(&old_source, fission_regions));

    info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, fabs((wf->get_keff() - k_new) / k_new));

    // Stopping criterion.
    if (fabs((wf->get_keff() - k_new) / k_new) < tol) eigen_done = true;

    // Update the final eigenvalue.
    wf->update_keff(k_new);

    it++;
        
    // Store the new eigenvector approximation in the result.
    for (int i = 0; i < solutions.size(); i++)  
      solutions[i]->copy(new_solutions[i]); 
  }
  while (!eigen_done);
  
  // Free memory.
  for (int i = 0; i < solutions.size(); i++) 
    delete new_solutions[i];
  
  return it;
}
