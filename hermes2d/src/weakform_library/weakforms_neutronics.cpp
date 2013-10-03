#include "weakforms_neutronics.h"

namespace Hermes { namespace Hermes2D {

void StationaryPicardSolver::set_tolerance(double tolerance_, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
{
  if (!handleMultipleTolerancesAnd)
    this->clear_tolerances();
  
  NonlinearSolver<double>::set_tolerance(tolerance_, toleranceType, handleMultipleTolerancesAnd);
}
    
void StationaryPicardSolver::solve(double *coeff_vec)
{
  this->dp->set_do_not_use_cache();
  int ndof = Space<double>::get_num_dofs(this->dp->get_spaces());
  
  bool _delete_coeff_vec = false;
  
  if(coeff_vec == NULL)
  {
    coeff_vec = (double*) malloc(ndof * sizeof(double));
    
    for (int i = 0; i < ndof; i++)
      coeff_vec[i] = 1.0;
    
    _delete_coeff_vec = true;
  }
  
  unsigned int it = 0;
  this->set_parameter_value(this->p_iteration, &it);
  
  this->init_solving(coeff_vec);
  
  this->delete_coeff_vec = _delete_coeff_vec;

  this->init_anderson();

  unsigned int vec_in_memory = 1;   // There is already one vector in the memory.
  this->set_parameter_value(this->p_vec_in_memory, &vec_in_memory);
  it++;
  
  if (this->tolerance_set[2])
  {
    // Chances are residual has already been assembled in 'init_solving'.
    if (!have_rhs)
    {
      this->conditionally_assemble(sln_vector, false, true);
      have_rhs = true;
    }
    
    initial_residual_norm = calculate_residual_norm();
  }
  
  while (true)
  {
    // Handle the event of step beginning.
    if(!this->on_step_begin())
    {
      this->info("\tPicard: aborted.");
      this->finalize_solving(coeff_vec);
      return;
    }

    // Assemble the residual and also jacobian when necessary (nonconstant jacobian, not reusable, ...).
    this->conditionally_assemble(coeff_vec, false, !have_rhs);
    if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);
    
    this->process_matrix_output(this->jacobian, it); 
    this->process_vector_output(this->residual, it);

    // Solve the linear system.
    Solvers::IterSolver<double>* iter_solver = dynamic_cast<Solvers::IterSolver<double>*>(matrix_solver);
    
    if (iter_solver)
    {
      if (dynamic_solver_tolerance)
        iter_solver->set_tolerance(
          std::min(0.1, (it==1 && this->tolerance_set[2]) ? initial_residual_norm : calculate_residual_norm())
        );
    }
    else
    {
      // Always factorize for the first time or in the case when Jacobian is updated on an algebraic level 
      // (new assembling is not needed, but new factorization is). 
      if (jacobian_change_on_algebraic_level)
        this->matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
      else
        this->matrix_solver->set_reuse_scheme(Solvers::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
    }

    this->matrix_solver->solve(this->sln_vector);
    this->handle_UMFPACK_reports();
    
    jacobian_change_on_algebraic_level = false;

    memcpy(this->sln_vector, this->matrix_solver->get_sln_vector(), sizeof(double)*ndof);
    have_rhs = false;

    if(!this->on_step_end())
    {
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      this->on_finish();
      return;
    }
    
    if (this->tolerance_set[2] && !have_rhs)
    {
      dp->assemble(sln_vector, residual);
      have_rhs = true;
    }
    
    if (this->converged())
    {
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      this->on_finish();
      return;
    }
    
    this->handle_previous_vectors(vec_in_memory);
    
    if(this->anderson_is_on && !this->on_step_end())
    {
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      this->on_finish();
      return;
    }

    // Increase counter of iterations.
    it++;
    
    // Renew the last iteration vector.
    memcpy(coeff_vec, this->sln_vector, ndof*sizeof(double));
    
    if (it >= this->max_allowed_iterations)
    {
      throw Exceptions::ValueException("iterations", it, this->max_allowed_iterations);
      this->deinit_solving(coeff_vec);
      this->deinit_anderson();
      this->on_finish();
      return;
    }
  }
}

void StationaryPicardSolver::use_dynamic_solver_tolerance(bool to_set)
{
  if (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != Hermes::SOLVER_AZTECOO)
    warn("Decreasing solver tolerance may be used only with an iterative solver (this currently means setting Hermes::matrixSolverType = Hermes::SOLVER_AZTECOO).");
  dynamic_solver_tolerance = to_set;
}

bool StationaryPicardSolver::converged()
{  
  double rel_err = std::numeric_limits<double>::max();
  double rel_res = std::numeric_limits<double>::max();
  
  if (this->tolerance_set[2])
  {
    rel_res = calculate_residual_norm() / initial_residual_norm;
    
    // Output for the user.
    this->info("\tPicard: iteration %d, nDOFs %d, relative residual %g%%", this->get_current_iteration_number(), this->ndof, rel_res * 100);

    if(rel_res < 1e-12)
      this->warn("\tPicard: a very small error threshold met, the loop should end.");
  }
  
  if (this->tolerance_set[6])
  {  
    // This is the new sln_vector.
    double* new_sln_vector = this->matrix_solver->get_sln_vector();
    
    double new_norm = get_l2_norm(new_sln_vector, this->ndof);

    // sln_vector still stores the old solution.
    // !!!! coeff_vec stores the Anderson-generated previous solution.
    double abs_error = 0.;
    for (int i = 0; i < this->ndof; i++)
      abs_error += std::abs((this->sln_vector[i] - new_sln_vector[i]) * (this->sln_vector[i] - new_sln_vector[i]));
    abs_error = std::sqrt(abs_error);
    
    rel_err = abs_error / new_norm;
    
    // Output for the user.
    this->info("\tPicard: iteration %d, nDOFs %d, relative error %g%%", this->get_current_iteration_number(), this->ndof, rel_err * 100);
  }
  
  bool conv;
  if (this->handleMultipleTolerancesAnd)
  {
    conv = true;
    
    if (this->tolerance_set[2])
      conv = conv && (rel_res < this->tolerance[2]);
    if (this->tolerance_set[6])
      conv = conv && (rel_err < this->tolerance[6]);
  }
  else
  {
    if (this->tolerance_set[2])
      conv = (rel_res < this->tolerance[2]);
    if (this->tolerance_set[6])
      conv = (rel_err < this->tolerance[6]);
  }
  
  return conv;
}

namespace Neutronics
{  
  KeffEigenvalueIteration::~KeffEigenvalueIteration()
  { 
    if (unshifted_jacobian != jacobian) 
      delete unshifted_jacobian;
    // jacobian will be deleted by the destructor of Solver
    
    if (production_matrix)
      delete production_matrix;
    if (production_dp)
      delete production_dp;
    if (Ax)
      delete [] Ax;
  }
    
  bool KeffEigenvalueIteration::on_initialization()
  {
    // NOTE: It is neccessary to call assign_dofs(spaces) before this method.
    
    if (shift_strategy == RAYLEIGH_QUOTIENT_SHIFT && !rayleigh)
    {
      warn("Using Rayleigh quotient shifting strategy requires KeffEigenvalueIteration::rayleigh==true.");
      use_rayleigh_quotient(true);
    }
    
    if (production_wf)
    {
      if (!production_dp)
      {
        production_dp = new DiscreteProblem<double>();
        production_dp->set_weak_formulation(production_wf);
        production_dp->set_linear();
        production_dp->set_do_not_use_cache();
      }
            
      if (!production_matrix)
        production_matrix = create_matrix<double>();
      
      production_dp->set_spaces(this->dp->get_spaces());
      production_dp->assemble(production_matrix);
    }
    
    // Initial assembly of Jacobian, so that the initial shift or Rayleigh quotient can be computed.
    // Residual is computed as well, since it will be needed anyway.
    conditionally_assemble(this->sln_vector); 
    have_rhs = true;
    
    if (shift_strategy == FIXED_SHIFT)
      set_shift(fixed_shift);
        
    if (rayleigh || this->tolerance_set[2])
    {
      if (Ax)
        delete [] Ax;
      
      Ax = new double [ndof];
    }
    
    this->update();
    return true;
  }
  
  bool KeffEigenvalueIteration::on_step_begin()
  {
    this->old_keff = keff;
    have_Ax = false;
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
    double lambda = 0.0;
    
    if (rayleigh)
    {
      double *resv = new double[ndof];

      if (have_rhs)
        residual->extract(resv);
      else
        production_matrix->multiply_with_vector(sln_vector, resv);
      
      double inv_norm = 1./get_l2_norm(resv, ndof);
      for (int i = 0; i < ndof; i++)
      {
        sln_vector[i] *= inv_norm;
        resv[i] *= inv_norm;
      }
      
      unshifted_jacobian->multiply_with_vector(sln_vector, Ax);
            
      for (int i = 0; i < ndof; i++) 
      {
        double Mxi = resv[i];
        lambda += Ax[i]*Mxi;
      }
      
      residual->set_vector(resv);
      delete [] resv;
      
      have_rhs = true;
      have_Ax = true;
            
      keff = 1./lambda;
      
      if (shift_strategy == RAYLEIGH_QUOTIENT_SHIFT)
        set_shift(lambda);
    }
    else
    {
      for (int i = 0; i < ndof; i++) 
        lambda += sqr(this->sln_vector[i]);
      
      lambda = 1./sqrt(lambda);
      
      for (int i = 0; i < ndof; i++) 
        this->sln_vector[i] *= lambda;
      
      keff = 1./(lambda + fixed_shift);
      
      if (production_wf)
      {
        // Save some time by using the mat-vec product instead of assembling the rhs.
        double *resv = new double[ndof];
        production_matrix->multiply_with_vector(sln_vector, resv);
        residual->set_vector(resv);
        have_rhs = true;
        delete [] resv;
      }
    }
  }
  
  bool KeffEigenvalueIteration::on_finish()
  {
    if (jacobian != unshifted_jacobian)
    {
      // If shifted jacobian has been used, return to the original one and reset the
      // matrix solver appropriately.
      delete jacobian;
      jacobian = unshifted_jacobian;
      matrix_solver = Solvers::create_linear_solver<double>(jacobian, residual);
      matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    }
    
    // Multiply the resulting dominant eigenvector by -1 if it is negative.
    double s = 0.0;
    for (int i = 0; i < ndof; i++) 
      s += this->sln_vector[i];
    
    if (s < 0)
      for (int i = 0; i < ndof; i++) 
        this->sln_vector[i] *= -1;
  }
  
  void KeffEigenvalueIteration::set_shift(double shift)
  {
    if (this->get_parameter_value(p_iteration) < num_unshifted_iterations)
      return;
    
    if (jacobian == unshifted_jacobian)
      unshifted_jacobian = jacobian->duplicate();
    else
    {
      delete jacobian;
      jacobian = unshifted_jacobian->duplicate();
      this->matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    }
    
    SparseMatrix<double> *shift_prod = production_matrix->duplicate();
    shift_prod->multiply_with_Scalar(-shift);
    jacobian->add_sparse_matrix(shift_prod);
    delete shift_prod;
    
    jacobian_change_on_algebraic_level = true;
  }
  
  void KeffEigenvalueIteration::set_spaces(const Hermes::vector<SpaceSharedPtr<double> >& spaces)
  {
    Solver<double>::set_spaces(spaces);
    if (production_dp)
      production_dp->set_spaces(spaces);
  }
  
  void KeffEigenvalueIteration::set_spectral_shift_strategy(ShiftStrategies strategy, int num_unshifted_iterations, double fixed_shift)
  {
    shift_strategy = strategy;
    this->num_unshifted_iterations = num_unshifted_iterations;
    
    if (strategy)
      assert(production_wf);
    if (strategy == FIXED_SHIFT)
    {
      assert(fixed_shift >= 0.0);    
      this->fixed_shift = fixed_shift;
    }
  }
  
  void KeffEigenvalueIteration::set_fixed_spectral_shift(double fixed_shift)
  {
    shift_strategy = FIXED_SHIFT;
    assert(production_wf);
    assert(fixed_shift >= 0.0);
    this->fixed_shift = fixed_shift;
  }
  
  void KeffEigenvalueIteration::use_rayleigh_quotient(bool to_set)
  {
    if (to_set) assert(production_wf);
    rayleigh = to_set;
  }
  
  double KeffEigenvalueIteration::calculate_residual_norm()
  {
    if (!have_rhs) // e.g. when residual norm is needed after an iteration which was monitored differently
    {
      dp->assemble(sln_vector, residual);
      have_rhs = true;
    }
    
    if (!Ax)  // the case above
      Ax = new double [ndof];  
    
    if (!have_Ax) // the case above, or when R. quot. was not used
      this->unshifted_jacobian->multiply_with_vector(this->sln_vector, Ax);
    
    double res = 0.0;
    double r = 1./keff;
    
    for (int i = 0; i < ndof; i++)
    {
      double Mxi = this->residual->get(i);
      res += (Ax[i] - r * Mxi) * (Ax[i] - r * Mxi);
    }
        
    return sqrt(res);
  }
  
  void KeffEigenvalueIteration::init_solving(double* coeff_vec)
  {
    PicardSolver<double>::init_solving(coeff_vec);
    // Copy back the scaled solution vector.
    memcpy(coeff_vec, this->sln_vector, ndof*sizeof(double));
  }
  
  double SourceIteration::calculate_residual_norm()
  {
    double *Ax = new double[ndof];
    this->jacobian->multiply_with_vector(this->sln_vector, Ax);
    
    double res = 0.0;
    
    for (int i = 0; i < ndof; i++)
    {
      double bi = this->residual->get(i);
      res += (Ax[i] - bi) * (Ax[i] - bi);
    }
    
    delete [] Ax;
    
    return sqrt(res);
  }
    
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 
