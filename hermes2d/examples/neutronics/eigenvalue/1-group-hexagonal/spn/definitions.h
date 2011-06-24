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
                   
class Views
{
  unsigned int n_moments, n_odd_moments, n_unknowns, n_equations, n_groups;
  bool display_meshes;
  
  ScalarView** sviews;
  OrderView** oviews;
  MeshView** mviews;
  
  static const std::string  base_title_flux;
  static const std::string  base_title_order;
  static const std::string  base_title_mesh;
  static const unsigned int MAX_SOLUTIONS_SETS = 10;
  
  MomentGroupFlattener mg;
  MomentFilter::Val*** moment_filters[MAX_SOLUTIONS_SETS];
  
  std::string itos(int t)
  {
    std::stringstream ss; ss << t;
    return ss.str();
  }
  
  public:
    Views(unsigned int spn_order, unsigned int G, bool display_meshes = false);
    ~Views();
    
    void show_meshes(Hermes::vector<Mesh*> meshes);
    void show_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set = 0);
    void show_orders(Hermes::vector<Space*> spaces);
    
    void inspect_meshes(Hermes::vector<Mesh*> meshes);
    void inspect_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set = 0);
    void inspect_orders(Hermes::vector<Space*> spaces);
};