////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion; 

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
public:
  CustomWeakForm( const MaterialPropertyMaps& matprop, const Hermes::vector<Solution*>& iterates,
                  double init_keff, const std::string& bdy_vacuum )
    : DefaultWeakFormSourceIteration(matprop, iterates, init_keff, HERMES_AXISYM_Y)
  {
    for (unsigned int g = 0; g < matprop.get_G(); g++)
    {
      add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(g, bdy_vacuum, HERMES_AXISYM_Y));
      add_vector_form_surf(new VacuumBoundaryCondition::Residual(g, bdy_vacuum, HERMES_AXISYM_Y));
    }
  }
};
