#include "hermes2d.h"
using namespace Hermes::Hermes2D; 

#include "weakforms_neutronics.h"
using namespace Neutronics::Diffusion;

class CustomWeakForm : public WeakForms::KeffEigenvalueProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop,
                   const Hermes::vector<Solution<double>*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, const std::string& bdy_vacuum);
};

void report_num_dof(const std::string& msg, const Hermes::vector<Space<double> *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
                    
