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

const std::string  Views::base_title_flux = "Neutron flux: group ";
const std::string  Views::base_title_order = "Polynomial orders: group ";
const std::string  Views::base_title_mesh = "Core mesh for group ";

Views::Views(unsigned int spn_order, unsigned int G, bool display_meshes) : display_meshes(display_meshes), mg(G)
{
  n_moments = spn_order+1;
  n_groups = G;
  n_unknowns = n_groups * n_moments;
  n_odd_moments = (n_moments+1)/2;
  n_equations = n_groups * n_odd_moments;
  
  sviews = new ScalarView* [n_unknowns];
  oviews = new OrderView* [n_equations];
  
  for (unsigned int g = 0; g < n_groups; g++)
  {
    std::string title_flux = base_title_flux + itos(g) + std::string(", moment ");
    std::string title_order = base_title_order + itos(g) + std::string(", moment ");
    for (unsigned int m = 0; m < n_moments; m++)
    {
      unsigned int i = mg.pos(m,g);
      
      sviews[i] = new ScalarView((title_flux + itos(m)).c_str(), new WinGeom(m*452, g*452, 450, 450));
      sviews[i]->show_mesh(false);
      sviews[i]->set_3d_mode(true);
      
      if (m%2) 
        oviews[mg.pos((m-1)/2,g)] = new OrderView((title_order + itos(m)).c_str(), new WinGeom(m*452, n_groups*452 + g*452, 450, 450));
    }
  }
  
  if (display_meshes)
  {
    mviews = new MeshView* [n_equations];
    
    for (unsigned int g = 0; g < n_groups; g++)
    {
      std::string title = base_title_mesh + itos(g) + std::string(", moment ");
      for (unsigned int m = 0; m < n_odd_moments; m++)
        mviews[mg.pos(m,g)] = new MeshView((title + itos(m)).c_str(), new WinGeom(m*352, g*352, 350, 350));
    }
  }
  else
    mviews = NULL;
  
  for (unsigned int i = 0; i < MAX_SOLUTIONS_SETS; i++)
  {
    moment_filters[i] = new MomentFilter::Val** [n_odd_moments];
    for (unsigned int m = 0; m < n_odd_moments; m++)
    {
      moment_filters[i][m] = new MomentFilter::Val* [n_groups];
      
      for (unsigned int g = 0; g < n_groups; g++)
        moment_filters[i][m][g] = NULL;
    }
  }
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
  
  for (unsigned int i = 0; i < MAX_SOLUTIONS_SETS; i++)
  {
    for (unsigned int m = 0; m < n_odd_moments; m++)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        if (moment_filters[i][m][g])
          delete moment_filters[i][m][g];

      delete [] moment_filters[i][m];
    }
    delete [] moment_filters[i];
  }
}

void Views::show_meshes(Hermes::vector< Mesh* > meshes)
{
  if (display_meshes)
    for (unsigned int g = 0; g < n_groups; g++)
      for (unsigned int m = 0; m < n_odd_moments; m++)
        mviews[mg.pos(m,g)]->show(meshes[mg.pos(m,g)]);
}

void Views::show_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set)
{
  if (solutions_set > MAX_SOLUTIONS_SETS)
    error("Change Views::MAX_SOLUTIONS_SETS and rebuild to allow visualizing more than 10 solutions sets.");
  
  for (unsigned int g = 0; g < n_groups; g++)
  {
    for (unsigned int m = 0; m < n_odd_moments; m++)
    {
      unsigned int i = mg.pos(m,g);
      unsigned int j = mg.pos(2*m,g);
      unsigned int k = mg.pos(2*m+1,g);
      
      if (moment_filters[solutions_set][m][g] == NULL)
        moment_filters[solutions_set][m][g] = new MomentFilter::Val(2*m, g, n_groups, solutions);
      
      sviews[j]->show(moment_filters[solutions_set][m][g]);
      sviews[k]->show(solutions[i]);
    }
  }
}

void Views::show_orders(Hermes::vector< Space* > spaces)
{
  for (unsigned int g = 0; g < n_groups; g++)
    for (unsigned int m = 0; m < n_odd_moments; m++)
      oviews[mg.pos(m,g)]->show(spaces[mg.pos(m,g)]);
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

void Views::inspect_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set)
{
  show_solutions(solutions, solutions_set);
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
