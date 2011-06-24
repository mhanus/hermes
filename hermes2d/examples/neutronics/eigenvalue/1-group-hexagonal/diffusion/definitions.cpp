#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop,
                                const Hermes::vector<Solution*>& iterates, 
                                const Hermes::vector<std::string>& fission_regions,
                                double init_keff, std::string bdy_vacuum )
  : DefaultWeakFormSourceIteration(matprop, iterates, fission_regions, init_keff)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_vector_form_surf(new VacuumBoundaryCondition::Residual(g, bdy_vacuum));
    add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(g, bdy_vacuum));    
  }
}

void report_num_dof(const std::string& msg, const Hermes::vector< Space* > spaces)
{
  std::stringstream ss;
  
  ss << msg << Space::get_num_dofs(spaces[0]);
  
  for (unsigned int i = 1; i < spaces.size(); i++)
    ss << " + " << Space::get_num_dofs(spaces[i]);
  
  if (spaces.size() > 1)
    ss << " = " << Space::get_num_dofs(spaces);
  
  info(ss.str().c_str());
}

const std::string Views::base_title_flux = "Neutron flux: group ";
const std::string Views::base_title_order = "Polynomial orders: group ";
const std::string Views::base_title_mesh = "Core mesh for group ";

Views::Views(unsigned int G, bool display_meshes) : display_meshes(display_meshes)
{
  n_groups = G;
  n_unknowns = n_groups;
  n_equations = n_groups;
  
  sviews = new ScalarView* [n_unknowns];
  oviews = new OrderView* [n_equations];
  
  for (unsigned int g = 0; g < n_groups; g++)
  {
    std::string title_flux = base_title_flux + itos(g);
    std::string title_order = base_title_order + itos(g);
    
    sviews[g] = new ScalarView(title_flux.c_str(), new WinGeom(0, g*452, 450, 450));
    sviews[g]->show_mesh(false);
    sviews[g]->set_3d_mode(true);
    oviews[g] = new OrderView(title_order.c_str(), new WinGeom(0, n_groups*452 + g*452, 450, 450));
  }
  
  if (display_meshes)
  {
    mviews = new MeshView* [n_equations];
    
    for (unsigned int g = 0; g < n_groups; g++)
    {
      std::string title = base_title_mesh + itos(g);
      mviews[g] = new MeshView(title.c_str(), new WinGeom(0, g*352, 350, 350));
    }
  }
  else
    mviews = NULL;
}

Views::~Views()
{
  if (sviews != NULL)
  {
    for (unsigned int i = 0; i < n_unknowns; i++)
      delete sviews[i];
    delete [] sviews;
  }
  
  if (oviews != NULL)
  {
    for (unsigned int i = 0; i < n_equations; i++)
      delete oviews[i];
    delete [] oviews;
  }
  
  if (mviews != NULL)
  {
    for (unsigned int i = 0; i < n_equations; i++)
      delete mviews[i];
    delete [] mviews;
  }
}

void Views::show_meshes(Hermes::vector< Mesh* > meshes)
{
  if (display_meshes)
    for (unsigned int g = 0; g < n_groups; g++)
      mviews[g]->show(meshes[g]);
}

void Views::show_solutions(Hermes::vector< Solution* > solutions)
{
  for (unsigned int g = 0; g < n_groups; g++)
    sviews[g]->show(solutions[g]);
}

void Views::show_orders(Hermes::vector< Space* > spaces)
{
  for (unsigned int g = 0; g < n_groups; g++)
    oviews[g]->show(spaces[g]);
}

void Views::inspect_meshes(Hermes::vector< Mesh* > meshes)
{
  if (display_meshes)
  {
    show_meshes(meshes);
    View::wait();
    
    for (unsigned int i = 0; i < n_equations; i++)
      delete mviews[i];
    delete [] mviews;
    
    mviews = NULL;
  }
}

void Views::inspect_solutions(Hermes::vector< Solution* > solutions)
{
  show_solutions(solutions);
  View::wait();
  
  for (unsigned int i = 0; i < n_unknowns; i++)
    delete sviews[i];
  delete [] sviews;
  
  sviews = NULL;
}

void Views::inspect_orders(Hermes::vector< Space* > spaces)
{
  show_orders(spaces);
  View::wait();
  
  for (unsigned int i = 0; i < n_equations; i++)
    delete oviews[i];
  delete [] oviews;

  oviews = NULL;
}
