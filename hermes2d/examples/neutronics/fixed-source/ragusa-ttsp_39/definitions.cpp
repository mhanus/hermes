#define HERMES_REPORT_ALL
#include "definitions.h"

// Fixed source problem with void conditions on all boundaries.
CustomWeakForm::CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N)
  : WeakForms::FixedSourceProblem(matprop, N)
{
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(Hermes::vector<std::string>("top","right"), m, N, 0, 1));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(Hermes::vector<std::string>("top","right"), m, n, 0, 1));
    }
}

