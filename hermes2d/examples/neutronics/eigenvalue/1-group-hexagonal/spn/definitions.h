#include "hermes2d.h"
using namespace Hermes::Hermes2D; 

#include "weakforms_neutronics.h"
using namespace Neutronics::SPN;

class CustomWeakForm : public WeakForms::KeffEigenvalueProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N,
                   const Hermes::vector<Solution<double>*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, const std::string& bdy_vacuum);
};
