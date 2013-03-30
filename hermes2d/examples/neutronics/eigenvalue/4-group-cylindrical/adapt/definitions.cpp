#define HERMES_REPORT_ALL

////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction<double> *sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  MeshSharedPtr mesh = sln->get_mesh();
  
  int n = 0;
  
  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o, e->get_mode());
    sln->set_quad_order(o, H2D_FN_VAL);
    double *uval = sln->get_fn_values();
    int np = quad->get_num_points(o, e->get_mode());
    
    for (int i = 0; i < np; i++)
      if (uval[i] < -1e-12)
        n++;
  }
  
  return n;
}