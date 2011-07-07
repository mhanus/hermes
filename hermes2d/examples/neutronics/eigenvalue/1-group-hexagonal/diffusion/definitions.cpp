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
