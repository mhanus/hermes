#define HERMES_REPORT_ALL
#include "definitions.h"

// Fixed source problem with void conditions on all boundaries.
CustomWeakForm::CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N)
#ifdef USE_SPN
  : WeakForms::FixedSourceProblem(matprop, N, HERMES_PLANAR, true)
#else // DIFFUSION
  : WeakForms::FixedSourceProblem(matprop, HERMES_PLANAR, true)
#endif
{
  for (unsigned int g = 0; g < G; g++)
  {
#ifdef USE_SPN    
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual("Vacuum", m, N, g, G));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian("Vacuum", m, n, g, G));    
    }
#else // DIFFUSION
    add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual("Vacuum", g));
    add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian("Vacuum", g));    
#endif
  }
}