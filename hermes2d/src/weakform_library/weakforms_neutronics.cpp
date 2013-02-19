#include "weakforms_neutronics.h"
#include "newton_solver.h"
#include <examples/neutronics/utils.h>

namespace Hermes { namespace Hermes2D { namespace Neutronics
{
  int keff_eigenvalue_iteration(const Hermes::vector<Solution<double> *>& solutions, 
                                Common::WeakForms::KeffEigenvalueProblem* wf, const Hermes::vector<const Space<double> *>& spaces,
                                MatrixSolverType matrix_solver, double tol_keff, double tol_flux, bool output_matrix_and_rhs, EMatrixDumpFormat output_fmt)
  {
    // Sanity checks.
    if (spaces.size() != solutions.size()) 
      ErrorHandling::error_function("Spaces and solutions supplied to power_iteration do not match.");
                                  
    // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
    // updating the eigenvalue. 
    Hermes::vector<Solution<double>*> new_solutions;
    for (unsigned int i = 0; i < solutions.size(); i++) 
      new_solutions.push_back(new Solution<double>(spaces[i]->get_mesh()));
      
    // Initial coefficient vector for the Newton's method.
    int ndof = Space<double>::get_num_dofs(spaces);
    
    DiscreteProblem<double> dp(wf, spaces);
    //dp.set_do_not_use_cache();
    NewtonSolver<double> solver(&dp);

    if (output_matrix_and_rhs)
    {
      solver.output_matrix(1);
      solver.output_rhs(1);
      solver.set_matrix_number_format("%1.16f");
      solver.set_rhs_number_format("%1.16f");
      solver.set_matrix_E_matrix_dump_format(output_fmt);
      solver.set_rhs_E_matrix_dump_format(output_fmt);
    }
        
    bool meshes_changed = true;
    bool eigen_done = false; int it = 0;
    do 
    {
      // The matrix doesn't change within the power iteration loop, so we don't have to reassemble the Jacobian again.
      try
      {
        if (output_matrix_and_rhs)
        {
          std::stringstream ss; ss << it;
          solver.set_rhs_filename(std::string("rhs_") + ss.str() + std::string("_"));
          solver.set_rhs_varname(std::string("b_") + ss.str() + std::string("_"));
        }
        solver.solve_keep_jacobian();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.print_msg();
        ErrorHandling::error_function("Newton's iteration failed.");
      }
      
      // Compute the eigenvalue for current iteration.
      double lambda = 0.0;
      double *v = solver.get_sln_vector();
      for (int i = 0; i < ndof; i++)
        lambda += sqr(v[i]);
      
      double k_new = sqrt(lambda);
      
      for (int i = 0; i < ndof; i++)
        v[i] /= k_new;
        
      Solution<double>::vector_to_solutions(v, spaces, new_solutions);
      
      double diff_keff = fabs(wf->get_keff() - k_new) / k_new;
      Hermes::Mixins::Loggable::Static::info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, diff_keff);
      
      // Stopping criteria.
      if (diff_keff < tol_keff) 
        eigen_done = true;

      // cout << "Iteration: " << it << ", flux diffces: " << endl; 
     
      // Update the final eigenvalue.
      wf->update_keff(k_new);
      
      // Update the eigenvector approximation.
      // FIXME: It seems that Space::construct_refined_spaces prevents automatic determination of meshes_changed
      //        through their seq number.
      wf->update_fluxes(new_solutions, meshes_changed);
      meshes_changed = false; // FIXME: true forces a reinitialization of the scalar flux filter used by the SPN weak form to
                             // evaluate the right hand side. If meshes are not changed, the unimesh created in that
                             // filter for the first time should not change either and this should be not neccessary.
  
      it++;
    }
    while (!eigen_done && it < 1);
    
    // Free memory.
    for (unsigned int i = 0; i < solutions.size(); i++) 
      delete new_solutions[i];
    
    return it;
  }
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 
