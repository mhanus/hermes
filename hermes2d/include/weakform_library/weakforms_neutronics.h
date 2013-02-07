#ifndef ___H2D_WEAK_FORMS_NEUTRONICS_H
#define ___H2D_WEAK_FORMS_NEUTRONICS_H

#include "neutronics/common_definitions.h"
#include "neutronics/material_properties.h"
#include "neutronics/support_classes.h"
#include "neutronics/weakform_parts_implementation.h"
#include "neutronics/weakforms.h"
#include "picard_solver.h"

namespace Hermes { namespace Hermes2D {
    
  class HERMES_API StationaryPicardSolver : public PicardSolver<double>
  {
    public:
      StationaryPicardSolver(DiscreteProblem<double>* dp)
        : PicardSolver<double>(dp) {};
      StationaryPicardSolver(const WeakForm<double>* wf, const Space<double>* space) 
        : PicardSolver<double>(wf, space) {};
      StationaryPicardSolver(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces) 
        : PicardSolver<double>(wf, spaces) {};
        
      virtual void solve(double *coeff_vec = NULL);
      
      using PicardSolver<double>::solve;  // Re-expose other overloaded 'solve' methods hidden by the above override.
  };
  
  namespace Neutronics
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
                                  bool output_matrix_and_rhs = false, EMatrixDumpFormat output_fmt = DF_HERMES_BIN);
    
    //TODO: non-vector version
      
    class HERMES_API KeffEigenvalueIteration : public StationaryPicardSolver
    {
      public:
        KeffEigenvalueIteration(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces)
          : StationaryPicardSolver(wf, spaces), keff(1.0)
        {};
        
        KeffEigenvalueIteration(DiscreteProblem<double>* dp)
          : StationaryPicardSolver(dp), keff(1.0)
        {};
        
        void on_initialization();
        void on_step_begin();
        void on_step_end();
        
        double get_keff() const { return keff; }
              
      private:
        double keff;
        double old_keff;
        
        void update();
    };
    
    class HERMES_API SourceIteration : public StationaryPicardSolver
    {
      public:
        SourceIteration(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces)
          : StationaryPicardSolver(wf, spaces)
        {};
        
        SourceIteration(DiscreteProblem<double>* dp)
          : StationaryPicardSolver(dp)
        {};
    };
    
  /* Neutronics */
  }
/* Hermes2D */
}
/* Hermes */
}  
#endif
