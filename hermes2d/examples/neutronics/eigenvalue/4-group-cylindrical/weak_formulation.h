////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"
using namespace Hermes::Hermes2D; 

#include "weakforms_neutronics.h"
using namespace Neutronics::Diffusion;

class CustomWeakForm : public WeakForms::KeffEigenvalueProblem
{
public:
  CustomWeakForm( const MaterialProperties::MaterialPropertyMaps& matprop, 
                  const Hermes::vector<Solution<double>*>& iterates,
                  double init_keff, const std::string& bdy_vacuum )
    : WeakForms::KeffEigenvalueProblem(matprop, iterates, init_keff, HERMES_AXISYM_Y)
  {
    for (unsigned int g = 0; g < matprop.get_G(); g++)
    {
      add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(bdy_vacuum, g, HERMES_AXISYM_Y));
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(bdy_vacuum, g, HERMES_AXISYM_Y));
    }
  }
};
