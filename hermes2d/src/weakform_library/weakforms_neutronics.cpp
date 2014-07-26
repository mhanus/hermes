#include "weakforms_neutronics.h"

namespace Hermes { namespace Hermes2D {

void StationaryPicardSolver::set_tolerance(double tolerance_, Solvers::NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
{
  if (toleranceType != Solvers::NonlinearConvergenceMeasurementType::ResidualNormRatioToInitial && toleranceType != Solvers::NonlinearConvergenceMeasurementType::SolutionChangeRelative)
    this->error("Only ResidualNormRatioToInitial or SolutionChangeRelative tolerance types are supported by StationaryPicardSolver.");

  NonlinearMatrixSolver<double>::set_tolerance(tolerance_, toleranceType, handleMultipleTolerancesAnd);
  measure_convergence_by_residual = this->tolerance_set[2];
}

Vector<double>* StationaryPicardSolver::duplicate_rhs()
{
	double *resv = new double[problem_size];
	linear_matrix_solver->get_rhs()->extract(resv);
	Vector<double>* rhs =  create_vector<double>();
	rhs->alloc(problem_size);
	rhs->set_vector(resv);
	delete [] resv;

	return rhs;
}


void StationaryPicardSolver::solve(double *coeff_vec)
{
  unsigned int it = 0;
  Hermes::vector<double> residual_norms;
  Hermes::vector<double> solution_norms;
  Hermes::vector<double> solution_change_norms;
  Hermes::vector<double> damping_factors;

  // Initial damping factor.
  damping_factors.push_back(this->manual_damping ? manual_damping_factor : initial_auto_damping_factor);

  // Link parameters.
  this->set_parameter_value(this->p_iteration, &it);
  this->set_parameter_value(this->p_residual_norms, &residual_norms);
  this->set_parameter_value(this->p_solution_norms, &solution_norms);
  this->set_parameter_value(this->p_solution_change_norms, &solution_change_norms);
  this->set_parameter_value(this->p_damping_factors, &damping_factors);

  this->init_solving(coeff_vec);
  this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(this->sln_vector, this->problem_size));

  it++;
  
  if (measure_convergence_by_residual || dynamic_solver_tolerance)
  {
    // Chances are residual has already been assembled in 'init_solving'.
    if (!have_rhs)
    {
      this->dp->assemble(this->sln_vector, this->get_jacobian(), this->get_residual());
      have_rhs = true;
      jacobian_reusable = true;
    }
    
    this->get_parameter_value(this->p_residual_norms).push_back(calculate_residual_norm());
  }
  
  while (true)
  {
    // Handle the event of step beginning.
    if(!this->on_step_begin())
    {
      this->info("\tPicard: aborted.");
      this->finalize_solving();
      return;
    }

    // Assemble the residual and also jacobian when necessary (nonconstant jacobian, not reusable, ...).
    // Always factorize for the first time or in the case when Jacobian is updated on an algebraic level 
    // (new assembling is not needed, but new factorization is). 
    if (this->jacobian_reusable)
    {
    	this->linear_matrix_solver->set_reuse_scheme(Solvers::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
    	if (!have_rhs)
    	{
    		this->dp->assemble(this->sln_vector, this->get_residual());
    		have_rhs = true;
    	}
    }
    else
    {
			this->linear_matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

			if (!have_rhs)
			{
				this->dp->assemble(this->sln_vector, this->get_jacobian(), this->get_residual());
				have_rhs = true;
			}
			else
				this->dp->assemble(this->sln_vector, this->get_jacobian());

			jacobian_reusable = true;
		}

    this->process_matrix_output(this->get_jacobian(), it);
    this->process_vector_output(this->get_residual(), it);

    // Solve the linear system.

    Solvers::IterSolver<double>* iter_solver = dynamic_cast<Solvers::IterSolver<double>*>(linear_matrix_solver);

    if (iter_solver && dynamic_solver_tolerance)
      iter_solver->set_tolerance( std::min(0.1, this->get_parameter_value(this->p_residual_norms).back()) );

    this->solve_linear_system();
    
    have_rhs = false;

    bool use_Anderson = this->anderson_is_on && (this->vec_in_memory >= this->num_last_vectors_used);
    if (use_Anderson)
    	memcpy(this->sln_vector, this->previous_Anderson_sln_vector, this->problem_size*sizeof(double));

    if(!this->on_step_end())
    {
			this->finalize_solving();
			return;
    }
    
    if (measure_convergence_by_residual || dynamic_solver_tolerance)
    	this->get_parameter_value(this->p_residual_norms).push_back(calculate_residual_norm());
    
    if (this->handle_convergence_state_return_finished(this->get_convergence_state()))
      return;

    // Increase counter of iterations.
    it++;

  }
}

void StationaryPicardSolver::use_dynamic_solver_tolerance(bool to_set)
{
  if (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != Hermes::SOLVER_AZTECOO)
    warn("Decreasing solver tolerance may be used only with an iterative solver (this currently means setting Hermes::matrixSolverType = Hermes::SOLVER_AZTECOO).");
  dynamic_solver_tolerance = to_set;
}

Solvers::NonlinearConvergenceState StationaryPicardSolver::get_convergence_state()
{
	if (this->get_current_iteration_number() >= this->max_allowed_iterations)
		return Solvers::NonlinearConvergenceState::AboveMaxIterations;
	if (this->converged())
		return Solvers::NonlinearConvergenceState::Converged;
	else
		return Solvers::NonlinearConvergenceState::NotConverged;

	return Solvers::NonlinearConvergenceState::Error;
}

bool StationaryPicardSolver::on_initialization()
{
    this->have_rhs = false;
    this->jacobian_reusable = false;
    return true;
}

bool StationaryPicardSolver::converged()
{
	double rel_err = std::numeric_limits<double>::max();
	double rel_res = std::numeric_limits<double>::max();

	this->info("\n\tPicard: iteration %d, nDOFs %d", this->get_current_iteration_number(), this->problem_size);

	if (measure_convergence_by_residual || dynamic_solver_tolerance)
	{
		double res_nrm = this->get_parameter_value(this->p_residual_norms).back();
		rel_res = res_nrm / this->get_parameter_value(this->p_residual_norms).front();

		// Output for the user.
		this->info("\t\tresidual norm: %g", res_nrm);
		this->info("\t\tresidual norm relative to initial: %g%%", rel_res * 100);

		if(rel_res < 1e-15)
			this->warn("\tPicard: a very small error threshold met, the loop should end.");
	}

	double solution_norm = this->get_parameter_value(this->p_solution_norms).back();
	double previous_solution_norm = this->get_parameter_value(this->p_solution_norms)[this->get_parameter_value(this->p_solution_norms).size() - 2];
	double solution_change_norm = this->get_parameter_value(this->p_solution_change_norms).back();
	this->info("\t\tsolution norm: %g,", solution_norm);
	this->info("\t\tsolution change norm: %g.", solution_change_norm);
	this->info("\t\trelative solution change: %g.", solution_change_norm / previous_solution_norm);

	if (!this->tolerance_set[2] && !this->tolerance_set[6])
	{
		if (this->handleMultipleTolerancesAnd)
			return true;
		else
			return false;
	}

	bool measure_convergence_by_relative_error = this->tolerance_set[6];
	bool conv = false;
	if (this->handleMultipleTolerancesAnd)
	{
		conv = true;

		if (this->measure_convergence_by_residual)
			conv = conv && (rel_res < this->tolerance[2]);
		if (measure_convergence_by_relative_error)
			conv = conv && (rel_err < this->tolerance[6]);
	}
	else
	{
		if (this->measure_convergence_by_residual)
			conv = conv || (rel_res < this->tolerance[2]);
		if (measure_convergence_by_relative_error)
			conv = conv || (rel_err < this->tolerance[6]);
	}

	return conv;
}

namespace Neutronics
{  
  KeffEigenvalueIteration::~KeffEigenvalueIteration()
  { 
    if (unshifted_jacobian != this->get_jacobian())
      delete unshifted_jacobian;
    // jacobian will be deleted by the destructor of Solver
    
    if (production_matrix)
      delete production_matrix;
    if (production_dp)
      delete production_dp;
    if (Ax)
      delete [] Ax;
  }

  void KeffEigenvalueIteration::init_solving(double* coeff_vec)
  {
    this->problem_size = Space<double>::assign_dofs(this->get_spaces());
    if (coeff_vec == nullptr)
    {
      coeff_vec = new double [this->problem_size];
      for (int i = 0; i < this->problem_size; i++)
        coeff_vec[i] = 1.0;
    }
    PicardMatrixSolver<double>::init_solving(coeff_vec);
    delete [] coeff_vec;
  }

  void KeffEigenvalueIteration::set_tolerance(double tolerance_, KeffEigenvalueIteration::ConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
  {
    Solvers::NonlinearConvergenceMeasurementType ctype;
    switch (toleranceType)
    {
      case KeffEigenvalueIteration::ConvergenceMeasurementType::ResidualNormRatioToInitial:
        ctype = Solvers::NonlinearConvergenceMeasurementType::ResidualNormRatioToInitial;
        break;
      case KeffEigenvalueIteration::ConvergenceMeasurementType::SolutionChangeRelative:
        ctype = Solvers::NonlinearConvergenceMeasurementType::SolutionChangeRelative;
        break;
      case KeffEigenvalueIteration::ConvergenceMeasurementType::EigenvalueRelative:
        // HACK: some convergence measurement must be set in NonlinearMatrixSolver; this one is ignored is StationaryPicardSolver::converged()
        ctype = Solvers::NonlinearConvergenceMeasurementType::SolutionChangeAbsolute;
        this->keff_tol = tolerance_;
        break;
    }

    NonlinearMatrixSolver<double>::set_tolerance(tolerance_, ctype, handleMultipleTolerancesAnd);
    measure_convergence_by_residual = this->tolerance_set[2];
  }

  bool KeffEigenvalueIteration::converged()
  {
    bool conv = StationaryPicardSolver::converged();

    double rel_err = abs(keff - old_keff) / keff;
    this->info("\n\t\tk_eff rel. error %g%%", rel_err * 100);

    if (keff_tol > 0)
      if (this->handleMultipleTolerancesAnd)
        conv = conv && (rel_err < keff_tol);
      else
        conv = conv || (rel_err < keff_tol);

    return conv;
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
      }
            
      if (!production_matrix)
        production_matrix = create_matrix<double>();
      
      production_dp->set_spaces(this->dp->get_spaces());
      production_dp->assemble(production_matrix);
    }
    
    // Initial assembly of Jacobian, so that the initial shift or Rayleigh quotient can be computed.
    // Residual is computed as well, since it will be needed anyway.
    this->dp->assemble(this->sln_vector, this->get_jacobian(), this->get_residual());
    have_rhs = true;
    jacobian_reusable = true;
    
    if (shift_strategy == FIXED_SHIFT)
      set_shift(fixed_shift);
        
    if (rayleigh || measure_convergence_by_residual || dynamic_solver_tolerance)
    {
      if (Ax)
        delete [] Ax;
      
      Ax = new double [problem_size];
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
    return true;
  }
  
  void KeffEigenvalueIteration::update()
  {
    double lambda = 0.0;
    
    if (rayleigh)
    {
      double *resv = new double[problem_size];

      if (have_rhs)
        this->get_residual()->extract(resv);
      else
        production_matrix->multiply_with_vector(sln_vector, resv, true);
      
      double inv_norm = 1./get_l2_norm(resv, problem_size);
      for (int i = 0; i < problem_size; i++)
      {
        sln_vector[i] *= inv_norm;
        resv[i] *= inv_norm;
      }
      
      unshifted_jacobian->multiply_with_vector(sln_vector, Ax, true);
            
      for (int i = 0; i < problem_size; i++)
      {
        double Mxi = resv[i];
        lambda += Ax[i]*Mxi;
      }
      
      this->get_residual()->set_vector(resv);
      delete [] resv;
      
      have_rhs = true;
      have_Ax = true;
            
      keff = 1./lambda;
      
      if (shift_strategy == RAYLEIGH_QUOTIENT_SHIFT)
        set_shift(lambda);
    }
    else
    {
      for (int i = 0; i < problem_size; i++)
        lambda += sqr(this->sln_vector[i]);

      lambda = 1./sqrt(lambda);

      for (int i = 0; i < problem_size; i++)
        this->sln_vector[i] *= lambda;

      keff = 1./(lambda + fixed_shift);

      if (production_wf)
      {
        // Save some time by using the mat-vec product instead of assembling the rhs.
        double *resv = new double[problem_size];
        production_matrix->multiply_with_vector(sln_vector, resv, true);
        this->get_residual()->set_vector(resv);
        have_rhs = true;
        delete [] resv;
      }
    }

    this->info("\n\t\tk_eff = %g", keff);
  }
  
  bool KeffEigenvalueIteration::on_finish()
  {
    if (this->get_jacobian() != unshifted_jacobian)
    {
      // If shifted jacobian has been used, return to the original one and reset the
      // matrix solver appropriately.
      this->get_jacobian()->zero();
      this->get_jacobian()->add_sparse_matrix(unshifted_jacobian);
      this->linear_matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    }
    
    // Multiply the resulting dominant eigenvector by -1 if it is negative.
    double s = 0.0;
    for (int i = 0; i < problem_size; i++)
      s += this->sln_vector[i];
    
    if (s < 0)
      for (int i = 0; i < problem_size; i++)
        this->sln_vector[i] *= -1;

    return true;
  }
  
  void KeffEigenvalueIteration::set_shift(double shift)
  {
    if (this->get_parameter_value(p_iteration) < num_unshifted_iterations)
      return;
    
    if (this->get_jacobian() == unshifted_jacobian)
      unshifted_jacobian = this->get_jacobian()->duplicate();
    else
    {
    	this->get_jacobian()->zero();
    	this->get_jacobian()->add_sparse_matrix(unshifted_jacobian);
      this->linear_matrix_solver->set_reuse_scheme(Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    }
    
    SparseMatrix<double> *shift_prod = production_matrix->duplicate();
    shift_prod->multiply_with_Scalar(-shift);
    this->get_jacobian()->add_sparse_matrix(shift_prod);
    delete shift_prod;
    
    jacobian_reusable = false;
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
      dp->assemble(sln_vector, this->get_residual());
      have_rhs = true;
    }
    
    if (!Ax)  // the case above
      Ax = new double [problem_size];
    
    if (!have_Ax) // the case above, or when R. quot. was not used
      this->unshifted_jacobian->multiply_with_vector(this->sln_vector, Ax, true);
    
    double res = 0.0;
    double r = 1./keff;
    
    for (int i = 0; i < problem_size; i++)
    {
      double Mxi = this->get_residual()->get(i);
      res += (Ax[i] - r * Mxi) * (Ax[i] - r * Mxi);
    }
        
    return sqrt(res);
  }
  
  double SourceIteration::calculate_residual_norm()
  {
  	if (!have_rhs)
		{
			dp->assemble(sln_vector, this->get_residual());
			have_rhs = true;
		}

    double *Ax = new double[problem_size];
    this->get_jacobian()->multiply_with_vector(this->sln_vector, Ax, true);
    
    double res = 0.0;
    
    for (int i = 0; i < problem_size; i++)
    {
      double bi = this->get_residual()->get(i);
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
