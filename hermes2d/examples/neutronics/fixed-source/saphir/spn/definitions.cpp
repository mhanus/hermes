#define HERMES_REPORT_ALL
#include "definitions.h"

// Fixed source problem with void conditions on all boundaries.
CustomWeakForm::CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N)
  : WeakForms::FixedSourceProblem(matprop, N, Hermes::Hermes2D::HERMES_PLANAR, true)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(m, N, g, G));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(m, n, g, G));    
    }
}

