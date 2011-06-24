#include "hermes2d.h"
#include <sstream>

using namespace WeakFormsNeutronics::Multigroup;
using namespace CompleteWeakForms::Diffusion;
using namespace ElementaryForms::Diffusion;

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   const Hermes::vector<Solution*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, std::string bdy_vacuum);
};

void report_num_dof(const std::string& msg, const Hermes::vector<Space *> spaces);
                    
class Views
{
  unsigned int n_unknowns, n_equations, n_groups;
  bool display_meshes;
  
  ScalarView** sviews;
  OrderView** oviews;
  MeshView** mviews;
  
  static const std::string  base_title_flux;
  static const std::string  base_title_order;
  static const std::string  base_title_mesh;
    
  std::string itos(int t)
  {
    std::stringstream ss; ss << t;
    return ss.str();
  }
  
  public:
    Views(unsigned int G, bool display_meshes = false);
    ~Views();
    
    void show_meshes(Hermes::vector<Mesh*> meshes);
    void show_solutions(Hermes::vector< Solution* > solutions);
    void show_orders(Hermes::vector<Space*> spaces);
    
    void inspect_meshes(Hermes::vector<Mesh*> meshes);
    void inspect_solutions(Hermes::vector< Solution* > solutions);
    void inspect_orders(Hermes::vector<Space*> spaces);
};