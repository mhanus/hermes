#include "weakforms_neutronics.h"
#include "solver/newton_solver.h"
#include <examples/neutronics/utils.h>

namespace Hermes { namespace Hermes2D {
    
void StationaryPicardSolver::solve(double *coeff_vec)
{
  int ndof = Space<double>::get_num_dofs(this->dp->get_spaces());
  
  bool _delete_coeff_vec = false;
  
  if(coeff_vec == NULL)
  {
    coeff_vec = (double*) malloc(ndof * sizeof(double));
    
    for (int i = 0; i < ndof; i++)
      coeff_vec[i] = 1.0;
    
    _delete_coeff_vec = true;
  }
  
  this->init_solving(coeff_vec);
  
  this->delete_coeff_vec = _delete_coeff_vec;

  this->init_anderson();

  unsigned int it = 1;
  unsigned int vec_in_memory = 1;   // There is already one vector in the memory.
  this->set_parameter_value(this->p_iteration, &it);
  this->set_parameter_value(this->p_vec_in_memory, &vec_in_memory);

  while (true)
  {
    this->on_step_begin();

    // Assemble the residual and also jacobian when necessary (nonconstant jacobian, not reusable, ...).
    this->conditionally_assemble(coeff_vec);

    this->process_matrix_output(this->jacobian, it); 
    this->process_vector_output(this->residual, it);

    // Solve the linear system.
    if(!this->matrix_solver->solve())
      throw Exceptions::LinearMatrixSolverException();

    memcpy(this->sln_vector, this->matrix_solver->get_sln_vector(), sizeof(double)*ndof);
    
    if (!this->on_step_end())
    {
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      return;
    }

    this->handle_previous_vectors(vec_in_memory);

    double rel_error = this->calculate_relative_error(coeff_vec);
    
    // Output for the user.
    this->info("\tPicard: iteration %d, nDOFs %d, relative error %g%%", it, ndof, rel_error * 100);

    // Find out the state with respect to all residual norms.
    PicardSolver<double>::ConvergenceState state = get_convergence_state(rel_error, it);

    switch(state)
    {
    case Converged:
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      return;
      break;

    case AboveMaxIterations:
      throw Exceptions::ValueException("iterations", it, this->max_allowed_iterations);
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      return;
      break;

    case Error:
      throw Exceptions::Exception("Unknown exception in PicardSolver.");
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      return;
      break;

    default:
      // The only state here is NotConverged which yields staying in the loop.
      break;
    }

    if(this->anderson_is_on && !this->on_step_end())
    {
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      return;
    }

    // Increase counter of iterations.
    it++;

    // Renew the last iteration vector.
    memcpy(coeff_vec, this->sln_vector, ndof*sizeof(double));
  }
}

namespace Neutronics
{
  int keff_eigenvalue_iteration(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, 
                                Common::WeakForms::KeffEigenvalueProblem* wf, const Hermes::vector<SpaceSharedPtr<double> >& spaces,
                                MatrixSolverType matrix_solver, double tol_keff, double tol_flux, bool output_matrix_and_rhs, EMatrixDumpFormat output_fmt)
  {
    // Sanity checks.
    if (spaces.size() != solutions.size()) 
      ErrorHandling::error_function("Spaces and solutions supplied to power_iteration do not match.");
                                  
    // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
    // updating the eigenvalue. 
    Hermes::vector<MeshFunctionSharedPtr<double> > new_solutions;
    for (unsigned int i = 0; i < solutions.size(); i++) 
      new_solutions.push_back(new Solution<double>(spaces[i]->get_mesh()));
    
    Common::SupportClasses::SourceFilter *new_source = wf->create_source_filter(new_solutions);
    Common::SupportClasses::SourceFilter *old_source = wf->create_source_filter(solutions);
      
    // Initial coefficient vector for the Newton's method.
    int ndof = Space<double>::get_num_dofs(spaces);
    
    DiscreteProblem<double> dp(wf, spaces);
    //dp.set_do_not_use_cache();
    NewtonSolver<double> solver(&dp);

    if (output_matrix_and_rhs)
    {
      solver.output_matrix(1);
      solver.output_rhs(1);
      solver.set_matrix_E_matrix_dump_format(output_fmt);
      solver.set_rhs_E_matrix_dump_format(output_fmt);
    }
        
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
      // The matrix doesn't change within the power iteration loop, so we don't have to reassemble the Jacobian again.
      try
      {
        if (output_matrix_and_rhs)
        {
          std::stringstream ss; ss << it;
          solver.set_rhs_filename(std::string("rhs_") + ss.str() + std::string("_"));
          solver.set_rhs_varname(std::string("b_") + ss.str() + std::string("_"));
        }
        solver.solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.print_msg();
        ErrorHandling::error_function("Newton's iteration failed.");
      }
      
      // Convert coefficients vector into a set of Solution pointers.
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, new_solutions);

      // Compute the eigenvalue for current iteration.
      double src_old = src_new;
      src_new = new_source->integrate();
      double k_new = wf->get_keff() * src_new / src_old;
      
      double diff_keff = fabs(wf->get_keff() - k_new) / k_new;
      Hermes::Mixins::Loggable::Static::info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, diff_keff);
      
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
          sln.copy(solutions[i].get());
          Solution<double> new_sln;
          new_sln.copy(new_solutions[i].get());
          
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
      meshes_changed = false; // FIXME: true forces a reinitialization of the scalar flux filter used by the SPN weak form to
                             // evaluate the right hand side. If meshes are not changed, the unimesh created in that
                             // filter for the first time should not change either and this should be not neccessary.
  
      it++;
    }
    while (!eigen_done);
    
    // Free memory.
    // Not needed when MeshFunctionSharedPtr is used instead of ordinary pointers.
    /*
    for (unsigned int i = 0; i < solutions.size(); i++) 
      delete new_solutions[i];
    */
    
    delete new_source;
    delete old_source;
    
    return it;
  }
    
  bool KeffEigenvalueIteration::on_initialization()
  {
    this->update();
    return true;
  }
  
  bool KeffEigenvalueIteration::on_step_begin()
  {
    this->old_keff = keff;
    return true;
  }
  
  bool KeffEigenvalueIteration::on_step_end()
  {
    this->update();
    
    double rel_err = abs(keff - old_keff) / keff;
    this->info("     k_eff = %g, rel. error %g%%", keff, rel_err * 100);
    
    if (keff_tol > 0)
      return (rel_err > keff_tol);
    
    return true;
  }
  
  void KeffEigenvalueIteration::update()
  {
    int ndof = Space<double>::get_num_dofs(this->dp->get_spaces());
    double lambda = 0.0;
    for (int i = 0; i < ndof; i++) 
      lambda += sqr(this->sln_vector[i]);
    
    keff = sqrt(lambda);
    
    for (int i = 0; i < ndof; i++) 
      this->sln_vector[i] /= keff;
  }
    
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 
