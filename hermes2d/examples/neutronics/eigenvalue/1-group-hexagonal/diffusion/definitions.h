#include "hermes2d.h"
using namespace Hermes::Hermes2D; 
using namespace Hermes::Mixins;

#include "weakforms_neutronics.h"
using namespace Neutronics::Diffusion;

class CustomWeakForm : public WeakForms::FixedSourceProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop,
                   const std::string& bdy_vacuum);
};
