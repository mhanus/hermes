#include "weakforms_neutronics.h"
#include "newton_solver.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{
  int keff_eigenvalue_iteration(const Hermes::vector<Solution<double> *>& solutions, 
                                Common::WeakForms::KeffEigenvalueProblem* wf, const Hermes::vector<const Space<double> *>& spaces,
                                MatrixSolverType matrix_solver, double tol_keff, double tol_flux)
  {
    // Sanity checks.
    if (spaces.size() != solutions.size()) 
      error_function("Spaces and solutions supplied to power_iteration do not match.");
                                  
    // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
    // updating the eigenvalue. 
    Hermes::vector<Solution<double>*> new_solutions;
    for (unsigned int i = 0; i < solutions.size(); i++) 
      new_solutions.push_back(new Solution<double>(spaces[i]->get_mesh()));
    
    Common::SupportClasses::SourceFilter *new_source = wf->create_source_filter(new_solutions);
    Common::SupportClasses::SourceFilter *old_source = wf->create_source_filter(solutions);
      
    // Initial coefficient vector for the Newton's method.
    int ndof = Space<double>::get_num_dofs(spaces);
    double* coeff_vec = new double[ndof];
    
    DiscreteProblem<double> dp(wf, spaces);
    NewtonSolver<double> solver(&dp, matrix_solver);
        
    // NOTE:
    //  According to the 'physical' definitions of keff as a fraction of neutron sources and sinks, keff should be 
    //  determined in an iteration like this: 
    //    k_new = k_old * src_new/src_old.  (1)
    //  This iteration can be expanded into the following form:
    //    k_new = k_0 * src_new/src_0,      (2)
    //  so if k_0 is set to src_0 at the beginning, we could potentially save one integration (of src_old) in each 
    //  iteration. Although the results are the same for both approaches (1) and (2) (up to a different eigenvector
    //  scaling), the calculation is slower surprisingly in case (2), even if compared to (1) where both old and new 
    //  sources are integrated in each iteration (this is actually not neccessary as you can see in the implementation
    //  below).
    //        
    double src_new = old_source->integrate();
    bool meshes_changed = true;
    bool eigen_done = false; int it = 0;
    do 
    {
      memset(coeff_vec, 0, ndof*sizeof(double));

      // The matrix doesn't change within the power iteration loop, so we don't have to reassemble the Jacobian again.
      try
      {
        solver.solve_keep_jacobian(coeff_vec, 1e-8, 3);
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.printMsg();
        error("Newton's iteration failed.");
      }
      
      // Convert coefficients vector into a set of Solution pointers.
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, new_solutions);

      // Compute the eigenvalue for current iteration.
      double src_old = src_new;
      src_new = new_source->integrate();
      double k_new = wf->get_keff() * src_new / src_old;
      
      double diff_keff = fabs(wf->get_keff() - k_new) / k_new;
      info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, diff_keff);
      
      // Stopping criteria.
      if (diff_keff < tol_keff) 
        eigen_done = true;

      // cout << "Iteration: " << it << ", flux diffces: " << endl; 
     
      if (tol_flux > 0)
      {
        for (unsigned int i = 0; i < solutions.size(); i++)
        {
          PostProcessor pp(wf->get_method_type(), wf->get_geom_type());
          
          // Normalize both flux iterates with the same criterion (unit integrated fission source).
          Solution<double> sln;
          sln.copy(solutions[i]);
          Solution<double> new_sln;
          new_sln.copy(new_solutions[i]);
          
          sln.multiply(1./src_old);
          new_sln.multiply(1./src_new);
          
          //cout << i << " = " << fabs(pp.integrate(new_sln) - pp.integrate(sln)) << ", ";
          
          // Compare the two solutions.
          if (fabs(pp.integrate(&new_sln) - pp.integrate(&sln)) >= tol_flux * pp.integrate(&new_sln))
            eigen_done = false;
          
          if (eigen_done == false)
            break;
        }
      }
      
      //cout << endl;

      // Update the final eigenvalue.
      wf->update_keff(k_new);

      // Update the eigenvector approximation.
      // FIXME: It seems that Space::construct_refined_spaces prevents automatic determination of meshes_changed
      //        through their seq number.
      wf->update_fluxes(new_solutions, meshes_changed);
      meshes_changed = true; // FIXME: This forces a reinitialization of the scalar flux filter used by the SPN weak form to
                             // evaluate the right hand side. If meshes are not changed, the unimesh created in that
                             // filter for the first time should not change either and this should be not neccessary. Currently,
                             // however, it doesn't work correctly.
  
      it++;
    }
    while (!eigen_done);
    
    // Free memory.
    for (unsigned int i = 0; i < solutions.size(); i++) 
      delete new_solutions[i];
    
    delete new_source;
    delete old_source;
    
    return it;
  }
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 