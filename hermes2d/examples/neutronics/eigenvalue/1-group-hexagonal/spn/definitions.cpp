#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N,
                                const Hermes::vector<Solution<double>*>& iterates, 
                                const Hermes::vector<std::string>& fission_regions,
                                double init_keff, const std::string& bdy_vacuum )
  : WeakForms::KeffEigenvalueProblem(matprop, N, iterates, fission_regions, init_keff)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(bdy_vacuum, m, N, g, G));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, m, n, g, G));    
    }
}

