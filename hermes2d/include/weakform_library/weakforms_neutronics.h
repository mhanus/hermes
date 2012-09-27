#ifndef ___H2D_WEAK_FORMS_NEUTRONICS_H
#define ___H2D_WEAK_FORMS_NEUTRONICS_H

#include "neutronics/common_definitions.h"
#include "neutronics/material_properties.h"
#include "neutronics/support_classes.h"
#include "neutronics/weakform_parts_implementation.h"
#include "neutronics/weakforms.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{
  /// \brief Power iteration method for finding the dominant eigenvalue. 
  ///
  /// Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
  /// has converged, also updating the global eigenvalue 'k_eff'.
  ///
  /// \param[in,out] solution     A set of Solution* pointers to solution components (neutron fluxes in each group). 
  ///                             Initial guess for the iteration on input, converged result on output.
  /// \param[in]     tol          Relative difference between two successive eigenvalue approximations that stops the iteration.
  /// \param[in]    matrix_solver Solver for the resulting matrix problem.
  ///
  /// \return  number of iterations needed for convergence within the specified tolerance.
  ///
  int keff_eigenvalue_iteration(const Hermes::vector<Solution<double> *>& solutions, 
                                Common::WeakForms::KeffEigenvalueProblem* wf, const Hermes::vector<const Space<double> *>& spaces,
                                MatrixSolverType matrix_solver = SOLVER_UMFPACK, double tol_keff = 1e-6, double tol_flux = 0,
                                bool output_matrix_and_rhs = false);

  //TODO: non-vector version
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}  
#endif
