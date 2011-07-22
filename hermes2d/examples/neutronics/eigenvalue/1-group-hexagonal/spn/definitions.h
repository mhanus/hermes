#include "hermes2d.h"
#include <sstream>

using namespace WeakFormsNeutronics::Multigroup;
using namespace CompleteWeakForms::SPN;
using namespace ElementaryForms::SPN;

using SupportClasses::SPN::MomentFilter;

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop, unsigned int N,
                   const Hermes::vector<Solution*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, std::string bdy_vacuum);
};

void report_num_dof(const std::string& msg, const Hermes::vector<Space *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);