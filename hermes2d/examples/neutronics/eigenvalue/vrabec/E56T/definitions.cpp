#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop,
                                const Hermes::vector<std::string>& bdy_vacuum )
  : WeakForms::FixedSourceProblem(matprop)
{
  /*for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(g, bdy_vacuum));
    add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(g, bdy_vacuum));    
  }*/
}

