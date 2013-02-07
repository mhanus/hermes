#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop,
                                const std::string& bdy_vacuum )
  : WeakForms::FixedSourceProblem(matprop)
{  
  for (unsigned int g = 0; g < matprop.get_G(); g++)
    add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, g));    
}
