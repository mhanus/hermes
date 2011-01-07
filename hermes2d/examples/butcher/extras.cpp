bool HERMES_RESIDUAL_AS_VECTOR = false;
bool solve_newton_butcher(ButcherTable* bt, 
                          scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
                          Vector* rhs, Hermes::Tuple<MeshFunction*>stage_solutions, 
                          double newton_tol, int newton_max_iter, bool verbose = true, 
                          double damping_coeff = 1.0, double max_allowed_residual_norm = 1e6)
{
  // Get number of stages.
  int num_stages = dp->get_num_spaces();

  // Space.
  Space* space = dp->get_space(0);

  // Get ndof.
  int ndof = space->get_num_dofs();

  // Stage vector of length num_stages * ndof, initialize with zeros.
  scalar* stage_vec = new scalar[num_stages*ndof];
  memset(stage_vec, 0, num_stages * ndof * sizeof(scalar));

  // Helper vector.
  scalar* vec = new scalar[ndof];

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (true)
  {
    // Prepare external solution for each stage.
    for (int r = 0; r < num_stages; r++) {
      memset(vec, 0, ndof * sizeof(scalar));
      double increment;
      for (int i = 0; i < ndof; i++) {
        increment = 0;
        for (int s = 0; s < num_stages; s++) {
          increment += bt->get_A(r, s) * stage_vec[s*ndof + i]; 
        }
        vec[i] = coeff_vec[i] + bt->get_time_step() * increment;
      }
      Solution::vector_to_solution(vec, space, (Solution*)stage_solutions[r]);
    } 

    // Calculating weight coefficients for blocks in the 
    // Stage jacobian matrix. 
    Table block_weights(num_stages);
    for (int r = 0; r < num_stages; r++) {
      for (int s = 0; s < num_stages; s++) {
        block_weights.set_A(r, s, bt->get_A(r, s) * bt->get_time_step());
      }
    }

    // Assemble the stage Jacobian matrix and residual vector.
    // Blocks that would be zeroed will not be assembled, and 
    // all assembled blocks will be weighted according to the 
    // Butcher's table.
    bool rhs_only = false;
    dp->assemble(NULL, matrix, rhs, rhs_only, &block_weights);

    // Add -1 to each diagonal element of the matrix.
    matrix->add_to_diagonal(-1);

    // Subtract stage_vec from rhs.
    for (int i = 0; i < num_stages*ndof; i++) rhs->add(i, -stage_vec[i]);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    rhs->change_sign();
    
    // Measure the residual norm.
    if (HERMES_RESIDUAL_AS_VECTOR) {
      // Calculate the l2-norm of residual vector.
      residual_norm = get_l2_norm(rhs);
    }
    else {
      // Translate the residual vector into a residual function (or multiple functions) 
      // in the corresponding finite element space(s) and measure their norm(s) there.
      // This is more meaningful since not all components in the coefficient vector 
      // have the same weight when translated into the finite element space.
      Hermes::Tuple<Solution*> residuals;
      Hermes::Tuple<bool> add_dir_lift;
      for (int i = 0; i < num_stages; i++) {
        residuals.push_back((Solution*)stage_solutions[i]);
        add_dir_lift.push_back(false);
      }
      Solution::vector_to_solutions(rhs, dp->get_spaces(), residuals, add_dir_lift);
      residual_norm = calc_norms(residuals);
    }

    // Info for the user.
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g", it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < num_stages*ndof; i++) stage_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  // If max number of iterations was exceeded, fail. 
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    // Delete helper vector.
    delete [] vec;
    return false;
  }

  // Calculate new coefficient vector using the stage vector and the Butcher's table.
  for (int i = 0; i < ndof; i++) {
    double increment = 0;
    for (int s = 0; s < num_stages; s++) {
      increment += bt->get_B(s) * stage_vec[s*ndof + i]; 
    }
    coeff_vec[i] += bt->get_time_step() * increment;
  } 

  // Delete helper vector.
  delete [] vec;

  return true;
}
