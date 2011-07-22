#define HERMES_REPORT_ALL
#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop, unsigned int N,
                                const Hermes::vector<Solution*>& iterates, 
                                const Hermes::vector<std::string>& fission_regions,
                                double init_keff, std::string bdy_vacuum )
  : DefaultWeakFormSourceIteration(matprop, N, iterates, fission_regions, init_keff)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new VacuumBoundaryCondition::Residual(m, N, g, bdy_vacuum, matprop));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(m, n, g, bdy_vacuum, matprop));    
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

void report_errors(const std::string& msg, const Hermes::vector< double > errors)
{
  std::stringstream ss;
  ss << msg;
  
  for (unsigned int i = 0; i < errors.size()-1; i++)
    ss << errors[i]*100 << "%%, ";
  
  ss << errors.back()*100 << "%%";
  
  info(ss.str().c_str());
}

