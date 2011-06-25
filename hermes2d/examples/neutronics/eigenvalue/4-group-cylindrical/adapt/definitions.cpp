#define HERMES_REPORT_ALL

////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

scalar ErrorForm::value(int n, double *wt, Func<scalar> *u_ext[],
                        Func<scalar> *u, Func<scalar> *v, Geom<double> *e,
                        ExtData<scalar> *ext) const
{
  switch (projNormType)
  {
    case HERMES_L2_NORM:
      return l2_error_form_axisym<double, scalar>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return h1_error_form_axisym<double, scalar>(n, wt, u_ext, u, v, e, ext);
    default:
      error("Only the H1 and L2 norms are currently implemented.");
      return 0.0;
  }
}

Ord ErrorForm::ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                   ExtData<Ord> *ext) const
{
  switch (projNormType)
  {
    case HERMES_L2_NORM:
      return l2_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return h1_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    default:
      error("Only the H1 and L2 norms are currently implemented.");
      return Ord();
  }
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  int n = 0;
  
  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
    
    for (int i = 0; i < np; i++)
      if (uval[i] < -1e-12)
        n++;
  }
  
  return n;
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
    ss << errors[i] << "%%, ";
  
  ss << errors.back() << "%%";
  
  info(ss.str().c_str());
}