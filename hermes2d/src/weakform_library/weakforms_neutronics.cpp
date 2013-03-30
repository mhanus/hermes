#include "weakforms_neutronics.h"
#include "newton_solver.h"
#include <examples/neutronics/utils.h>

namespace Hermes { namespace Hermes2D {
  
void StationaryPicardSolver::solve(double *coeff_vec)
{
  this->check();
  this->tick();
      
  // Sanity check.
  if(this->num_last_vectors_used < 1)
    throw Hermes::Exceptions::Exception("PicardSolver: Bad number of last iterations to be used (must be at least one).");

  // Preliminaries.
  
  int num_spaces = static_cast<DiscreteProblem<double>*>(this->dp)->get_spaces().size();
  int ndof = static_cast<DiscreteProblem<double>*>(this->dp)->get_num_dofs();
  Hermes::vector<SpaceSharedPtr<double> > spaces = static_cast<DiscreteProblem<double>*>(this->dp)->get_spaces();
  Hermes::vector<bool> add_dir_lift;
  for(unsigned int i = 0; i < spaces.size(); i++)
    add_dir_lift.push_back(false);
  
  // Delete solution vector if there is any.
  if(this->sln_vector != NULL)
  { 
    delete [] this->sln_vector;
    this->sln_vector = NULL;
  }

  this->sln_vector = new double[ndof];
  
  bool delete_coeff_vec = false;
  bool default_initial_vector = false;
  if(coeff_vec == NULL)
  {
    coeff_vec = new double[ndof];
    for (int i = 0; i < ndof; i++) coeff_vec[i] = 1.0;
    delete_coeff_vec = true;
    default_initial_vector = true;
  }

  memcpy(this->sln_vector, coeff_vec, ndof*sizeof(double));

  // Save the coefficient vector, it will be used to calculate increment error
  // after a new coefficient vector is calculated.
  double* last_iter_vector = new double[ndof];
  for (int i = 0; i < ndof; i++)
    last_iter_vector[i] = this->sln_vector[i];

  // If Anderson is used, allocate memory for vectors and coefficients.
  double** previous_vectors = NULL;      // To store num_last_vectors_used last coefficient vectors.
  double* anderson_coeffs = NULL;        // To store num_last_vectors_used - 1 Anderson coefficients.
  if (anderson_is_on)
  {
    previous_vectors = new double*[num_last_vectors_used];
    for (int i = 0; i < num_last_vectors_used; i++) previous_vectors[i] = new double[ndof];
    anderson_coeffs = new double[num_last_vectors_used-1];
  }

  // If Anderson is used, save the initial coefficient vector in the memory.
  if (anderson_is_on)
    for (int i = 0; i < ndof; i++) previous_vectors[0][i] = this->sln_vector[i];

  int it = 1;
  int vec_in_memory = 1;   // There is already one vector in the memory.

  this->on_initialization();
  
  double rel_error;
  while (true)
  {
    this->on_step_begin();
    
    if (it == 1)
    {
      this->linear_solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
      
      this->info("Assembling the matrix and the right-hand side of the discrete problem.");
      static_cast<DiscreteProblem<double>*>(this->dp)->assemble(last_iter_vector, this->matrix, this->rhs);
      
      if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= it))
      {
        char* fileName = new char[this->matrixFilename.length() + 5];
        if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
          sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), it);
        else
          sprintf(fileName, "%s%i", this->matrixFilename.c_str(), it);
        FILE* matrix_file = fopen(fileName, "w+");

        matrix->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
        fclose(matrix_file);
        delete [] fileName;
      }
      if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
      {
        char* fileName = new char[this->RhsFilename.length() + 5];
        if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
          sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
        else
          sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
        FILE* rhs_file = fopen(fileName, "w+");
        rhs->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
        fclose(rhs_file);
        delete [] fileName;
      }
    }
    else
    { 
      this->linear_solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
      
      this->info("Assembling the right-hand side of the discrete problem.");
      static_cast<DiscreteProblem<double>*>(this->dp)->assemble(last_iter_vector, this->rhs);
      
      if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
      {
        char* fileName = new char[this->RhsFilename.length() + 5];
        if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
          sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
        else
          sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
        FILE* rhs_file = fopen(fileName, "w+");
        rhs->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
        fclose(rhs_file);
        delete [] fileName;
      }
    } 
     
    this->info("Solving the linear system.");
    if(!linear_solver->solve())
      throw Exceptions::LinearMatrixSolverException();
    
    memcpy(this->sln_vector, this->linear_solver->get_sln_vector(), sizeof(double)*ndof);

    // If Anderson is used, store the new vector in the memory.
    if (anderson_is_on)
    {
      // If memory not full, just add the vector.
      if (vec_in_memory < num_last_vectors_used)
      {
        for (int i = 0; i < ndof; i++) previous_vectors[vec_in_memory][i] = this->sln_vector[i];
        vec_in_memory++;
      }
      else
      {
        // If memory full, shift all vectors back, forgetting the oldest one.
        // Save this->sln_vector[] as the newest one.
        double* oldest_vec = previous_vectors[0];
        for (int i = 0; i < num_last_vectors_used-1; i++) previous_vectors[i] = previous_vectors[i + 1];
        previous_vectors[num_last_vectors_used-1] = oldest_vec;
        for (int j = 0; j < ndof; j++) previous_vectors[num_last_vectors_used-1][j] = this->sln_vector[j];
      }
    }

    // If there is enough vectors in the memory, calculate Anderson coeffs.
    if (anderson_is_on && vec_in_memory >= num_last_vectors_used)
    {
      // Calculate Anderson coefficients.
      calculate_anderson_coeffs(previous_vectors, anderson_coeffs, num_last_vectors_used, ndof);

      // Calculate new vector and store it in this->sln_vector[].
      for (int i = 0; i < ndof; i++)
      {
        this->sln_vector[i] = 0;
        for (int j = 1; j < num_last_vectors_used; j++)
        {
          this->sln_vector[i] += anderson_coeffs[j-1] * previous_vectors[j][i] - (1.0 - anderson_beta) * anderson_coeffs[j-1] * (previous_vectors[j][i] - previous_vectors[j-1][i]);
        }
      }
    }
    
    this->on_step_end();
    
    // Calculate relative error between last_iter_vector[] and this->sln_vector[].
    // FIXME: this is wrong in the complex case (complex conjugation must be used).
    // FIXME: This will crash if norm of last_iter_vector[] is zero.
    double last_iter_vec_norm = 0;
    for (int i = 0; i < ndof; i++)
      last_iter_vec_norm += sqr(last_iter_vector[i]);

    last_iter_vec_norm = sqrt(last_iter_vec_norm);

    double abs_error = 0;
    for (int i = 0; i < ndof; i++) abs_error += sqr(this->sln_vector[i] - last_iter_vector[i]);
    abs_error = sqrt(abs_error);

    rel_error = abs_error / last_iter_vec_norm;

    // Output for the user.
    if(default_initial_vector)
    {
      default_initial_vector = false;
      this->info("\tPicard: iteration %d, nDOFs %d, starting from vector of all ones.", it, ndof);
    }
    else
      this->info("\tPicard: iteration %d, nDOFs %d, relative error %1.16g%%", it, ndof, rel_error * 100);
        
    // Stopping because error is sufficiently low.
    if(rel_error < tol)
    {
      delete [] last_iter_vector;
      // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
      if (anderson_is_on)
      {
        for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
        delete [] previous_vectors;
        delete [] anderson_coeffs;
      }
      
      static_cast<DiscreteProblem<double>*>(this->dp)->have_matrix = false;

      this->tick();
      this->info("\tPicard: solution duration: %f s.\n", this->last());
      this->on_finish();
      return;
    }
    
    // Stopping because maximum number of iterations reached.
    if(it >= max_iter)
    {
      delete [] last_iter_vector;
      // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
      if (anderson_is_on)
      {
        for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
        delete [] previous_vectors;
        delete [] anderson_coeffs;
      }
      static_cast<DiscreteProblem<double>*>(this->dp)->have_matrix = false;

      this->tick();
      this->info("Picard: solution duration: %f s.\n", this->last());

      this->on_finish();
      throw Hermes::Exceptions::Exception("Picard: maximum allowed number of Picard iterations exceeded.");
      return;
    }
    
    // Increase counter of iterations.
    it++;
    
    // Renew the last iteration vector.
    for (int i = 0; i < ndof; i++)
      last_iter_vector[i] = this->sln_vector[i];
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
        solver.solve_keep_jacobian();
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
      meshes_changed = false; // FIXME: true forces a reinitialization of the scalar flux filter used by the SPN weak form to
                             // evaluate the right hand side. If meshes are not changed, the unimesh created in that
                             // filter for the first time should not change either and this should be not neccessary.
  
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
    
  void KeffEigenvalueIteration::on_initialization()
  {
    this->update();
  }
  
  void KeffEigenvalueIteration::on_step_begin()
  {
    this->old_keff = keff;
  }
  
  void KeffEigenvalueIteration::on_step_end()
  {
    this->update();
    
    double rel_err = abs(keff - old_keff) / keff;
    this->info("     k_eff = %g, rel. error %g%%", keff, rel_err * 100);      
  }
  
  void KeffEigenvalueIteration::update()
  {
    int ndof = static_cast<DiscreteProblem<double>*>(this->dp)->get_num_dofs();
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
