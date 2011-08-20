#include "definitions.h"

void report_num_dof(const std::string& msg, const Hermes::vector< Space<double>* >& spaces)
{
  std::stringstream ss;
  
  ss << msg << spaces[0]->get_num_dofs();
  
  for (unsigned int i = 1; i < spaces.size(); i++)
    ss << " + " << spaces[i]->get_num_dofs();
  
  if (spaces.size() > 1)
    ss << " = " << Space<double>::get_num_dofs(spaces);
  
  info(ss.str().c_str());
}

std::string itos(int t)
{
  std::stringstream ss; ss << t;
  return ss.str();
}