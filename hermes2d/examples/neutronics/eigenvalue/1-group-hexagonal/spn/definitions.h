#include "hermes2d.h"
using namespace Hermes::Hermes2D; 
using namespace Hermes::Mixins;

#include "weakforms_neutronics.h"
using namespace Neutronics::SPN;

class CustomWeakForm : public WeakForms::FixedSourceProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N,
                   const std::string& bdy_vacuum);
};
