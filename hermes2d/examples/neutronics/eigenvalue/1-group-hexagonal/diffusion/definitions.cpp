#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop,
                                const Hermes::vector<Solution<double>*>& iterates, 
                                const Hermes::vector<std::string>& fission_regions,
                                double init_keff, const std::string& bdy_vacuum )
  : WeakForms::KeffEigenvalueProblem(matprop, iterates, fission_regions, init_keff)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(bdy_vacuum, g));
    add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, g));    
  }
}
