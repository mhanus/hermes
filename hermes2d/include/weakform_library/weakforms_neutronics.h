#ifndef ___H2D_WEAK_FORMS_NEUTRONICS_H
#define ___H2D_WEAK_FORMS_NEUTRONICS_H

#include "neutronics/common_definitions.h"
#include "neutronics/material_properties.h"
#include "neutronics/support_classes.h"
#include "neutronics/weakform_parts_implementation.h"
#include "neutronics/weakforms.h"
#include "../solver/picard_solver.h"

namespace Hermes { namespace Hermes2D {
    
  class HERMES_API StationaryPicardSolver : public PicardSolver<double>
  {
    protected:
      double initial_residual_norm;
      bool convergence_by_residual;
      bool have_rhs;
      bool dynamic_solver_tolerance;
      bool jacobian_change_on_algebraic_level;
      bool measure_convergence_by_residual;
      
    public:
      StationaryPicardSolver(DiscreteProblem<double>* dp)
        : PicardSolver<double>(dp), 
          measure_convergence_by_residual(false), have_rhs(false), dynamic_solver_tolerance(false), jacobian_change_on_algebraic_level(true) 
      {
        constant_jacobian = true;
      }
      
      StationaryPicardSolver(WeakForm<double>* wf, SpaceSharedPtr<double>& space) 
        : PicardSolver<double>(wf, space), 
          measure_convergence_by_residual(false), have_rhs(false), dynamic_solver_tolerance(false), jacobian_change_on_algebraic_level(true) 
      {
        constant_jacobian = true;
      }
      
      StationaryPicardSolver(WeakForm<double>* wf, const Hermes::vector<SpaceSharedPtr<double> >& spaces) 
        : PicardSolver<double>(wf, spaces),
          measure_convergence_by_residual(false), have_rhs(false), dynamic_solver_tolerance(false), jacobian_change_on_algebraic_level(true)
      {
        constant_jacobian = true;
      }
      
      virtual void solve(double *coeff_vec);
      void solve() { this->solve(NULL); }
      
      void use_dynamic_solver_tolerance(bool to_set = true);
      
      bool converged();
      
      virtual void set_tolerance(double tolerance_, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd = false);
      
      virtual double calculate_residual_norm() = 0;
      
      using PicardSolver<double>::solve;  // Re-expose other overloaded 'solve' methods hidden by the above override.
  };
  
  namespace Neutronics
  {   
    class HERMES_API KeffEigenvalueIteration : public StationaryPicardSolver
    {
      public:
        enum ShiftStrategies { NO_SHIFT = 0, FIXED_SHIFT, RAYLEIGH_QUOTIENT_SHIFT };
        
        KeffEigenvalueIteration(WeakForm<double>* wf, const Hermes::vector<SpaceSharedPtr<double> >& spaces, WeakForm<double>* prod_wf = NULL)
          : StationaryPicardSolver(wf, spaces), 
            keff(1.0), keff_tol(0.0), fixed_shift(0.0), rayleigh(false), 
            shift_strategy(NO_SHIFT), unshifted_jacobian(jacobian), num_unshifted_iterations(0),
            production_matrix(NULL), production_dp(NULL), production_wf(prod_wf), Ax(NULL), have_Ax(false)
        {};
        
        KeffEigenvalueIteration(DiscreteProblem<double>* dp, WeakForm<double>* prod_wf = NULL)
          : StationaryPicardSolver(dp), 
            keff(1.0), keff_tol(0.0), fixed_shift(0.0), rayleigh(false),
            shift_strategy(NO_SHIFT), unshifted_jacobian(jacobian), num_unshifted_iterations(0),
            production_matrix(NULL), production_dp(NULL), production_wf(prod_wf), Ax(NULL), have_Ax(false)
        {};
        
        virtual ~KeffEigenvalueIteration();
        
        bool on_initialization();
        bool on_step_begin();
        bool on_step_end();
        bool on_finish();
        
        double get_keff() const { return keff; }
        
        void set_keff_tol(double tol) { keff_tol = tol; }
        
        void set_production_weakform(WeakForm<double> *prod_wf) { production_wf = prod_wf; }
        void set_spectral_shift_strategy(ShiftStrategies strategy, int num_unshifted_iterations = 0, double fixed_shift = 0.0);
        void set_fixed_spectral_shift(double fixed_shift);
        
        void use_rayleigh_quotient(bool to_set = true);
        
        virtual void set_spaces(const Hermes::vector<SpaceSharedPtr<double> >& spaces);
        
        virtual void init_solving(double* coeff_vec);
        virtual double calculate_residual_norm();
        
      private:
        double keff;
        double old_keff;
        double keff_tol;
        
        ShiftStrategies shift_strategy;
        double fixed_shift;
        int num_unshifted_iterations;
        
        bool rayleigh;
        
        double *Ax;
        bool have_Ax;
        
        SparseMatrix<double>* unshifted_jacobian;
        SparseMatrix<double>* production_matrix;
        DiscreteProblem<double>* production_dp;
        WeakForm<double>* production_wf;
        
        void update();
        void set_shift(double shift);
    };
    
    class HERMES_API SourceIteration : public StationaryPicardSolver
    {
      public:
        SourceIteration(WeakForm<double>* wf, const Hermes::vector<SpaceSharedPtr<double> >& spaces)
          : StationaryPicardSolver(wf, spaces)
        {};
        
        SourceIteration(DiscreteProblem<double>* dp)
          : StationaryPicardSolver(dp)
        {};
        
        virtual double calculate_residual_norm();
    };
    
  /* Neutronics */
  }
/* Hermes2D */
}
/* Hermes */
}  
#endif
