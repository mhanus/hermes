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
      StationaryPicardSolver(DiscreteProblemLinear<double>* dp, Solution<double>* sln_prev_iter)
        : PicardSolver<double>(dp, sln_prev_iter) {};
      StationaryPicardSolver(DiscreteProblemLinear<double>* dp, Hermes::vector<Solution<double>* > slns_prev_iter) 
        : PicardSolver<double>(dp, slns_prev_iter) {};
      StationaryPicardSolver(const WeakForm<double>* wf, const Space<double>* space, Solution<double>* sln_prev_iter) 
        : PicardSolver<double>(wf, space, sln_prev_iter) {};
      StationaryPicardSolver(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces, Solution<double>* sln_prev_iter) 
        : PicardSolver<double>(wf, spaces, sln_prev_iter) {};
      StationaryPicardSolver(const WeakForm<double>* wf, const Space<double>* space, Hermes::vector<Solution<double>* > slns_prev_iter) 
        : PicardSolver<double>(wf, space, slns_prev_iter) {};
      StationaryPicardSolver(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces, Hermes::vector<Solution<double>* > slns_prev_iter) 
        : PicardSolver<double>(wf, spaces, slns_prev_iter) {};
        
      virtual void solve();
      
      virtual ~StationaryPicardSolver();
      
    protected:
      SparseMatrix<double>* mat;
      Vector<double>* rhs;
      LinearMatrixSolver<double>* linear_solver;
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
                                  MatrixSolverType matrix_solver = SOLVER_UMFPACK, double tol_keff = 1e-6, double tol_flux = 0);

    //TODO: non-vector version
    
    class HERMES_API KeffEigenvalueIteration : public StationaryPicardSolver
    {
      public:
        KeffEigenvalueIteration(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces, Hermes::vector<Solution<double>* > slns_prev_iter)
          : StationaryPicardSolver(wf, spaces, slns_prev_iter), keff_tol(0)
        {};
        
        KeffEigenvalueIteration(DiscreteProblemLinear<double>* dp, Hermes::vector<Solution<double>* > slns_prev_iter)
          : StationaryPicardSolver(dp, slns_prev_iter), keff_tol(0)
        {};
        
        void set_keff_tol(double tol);
        
        void onInitialization();
        void onStepBegin();
        void onStepEnd();
        void onFinish();
              
      private:
        double keff_tol;
        double keff;
        double old_keff;
        
        void update();
    };
    
    class HERMES_API SourceIteration : public StationaryPicardSolver
    {
      public:
        SourceIteration(const WeakForm<double>* wf, Hermes::vector<const Space<double>*> spaces, Hermes::vector<Solution<double>* > slns_prev_iter)
          : StationaryPicardSolver(wf, spaces, slns_prev_iter)
        {};
        
        SourceIteration(DiscreteProblemLinear<double>* dp, Hermes::vector<Solution<double>* > slns_prev_iter)
          : StationaryPicardSolver(dp, slns_prev_iter)
        {};
    };
    
  /* Neutronics */
  }
/* Hermes2D */
}
/* Hermes */
}  
#endif
