#define HERMES_REPORT_ALL
#include "definitions.h"

// Fixed source problem with void conditions on all boundaries.
CustomWeakForm::CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N)
  : WeakForms::FixedSourceProblem(matprop, N)
{
  for (unsigned int g = 0; g < G; g++)
    for (unsigned int m = 0; m < N_odd; m++)
    {
      add_vector_form_surf(new WeakFormParts::VacuumBoundaryCondition::Residual(m, N, g, matprop));
      for (unsigned int n = 0; n < N_odd; n++)
        add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(m, n, g, matprop));    
    }
}

// Utility functions that simplify repeated reporting of number of DOF during adaptivity.
void report_num_dof(const std::string& msg, const Hermes::vector< Space<double>* > spaces)
{
  std::stringstream ss;
  
  ss << msg << spaces[0]->get_num_dofs();
  
  for (unsigned int i = 1; i < spaces.size(); i++)
    ss << " + " << spaces[i]->get_num_dofs();
  
  if (spaces.size() > 1)
    ss << " = " << Space<double>::get_num_dofs(spaces);
  
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

std::string itos(int t)
{
  std::stringstream ss; ss << t;
  return ss.str();
}

