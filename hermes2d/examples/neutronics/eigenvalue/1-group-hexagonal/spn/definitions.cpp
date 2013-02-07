#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N,
                                const std::string& bdy_vacuum )
  : WeakForms::FixedSourceProblem(matprop, N)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, m, n, g, G));    
}

