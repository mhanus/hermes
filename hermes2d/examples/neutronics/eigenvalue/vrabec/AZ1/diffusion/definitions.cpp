#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop,
                                const Hermes::vector<std::string>& bdy_vacuum )
  : WeakForms::FixedSourceProblem(matprop)
{
  if (!bdy_vacuum.empty())
    for (unsigned int g = 0; g < matprop.get_G(); g++)
      add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, g));  
}

