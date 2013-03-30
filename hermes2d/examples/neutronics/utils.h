#include "hermes2d.h"
using namespace Hermes::Hermes2D;
using namespace Hermes::Mixins;

void report_num_dof(const std::string& msg, const Hermes::vector<SpaceSharedPtr<double> > spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);
