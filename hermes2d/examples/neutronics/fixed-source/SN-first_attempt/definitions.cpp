#include "definitions.h"
#include <cstdlib>

SNWeakForm::SNWeakForm(unsigned int N, const MaterialProperties::MaterialPropertyMaps& matprop,
                       const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries) 
  : WeakForm<double>(matprop.get_G()*N*(N+2)/2), N(N), M(N*(N+2)/2), G(matprop.get_G()), odata(SupportClasses::OrdinatesData(N, "lgvalues.txt"))
{  
    bool1 chi_nnz = matprop.get_fission_nonzero_structure();
   
    std::cout << odata;
        
    std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
    for ( ; material != matprop.get_materials_list().end(); ++material)
    {
      Hermes::vector<std::string> regions = matprop.get_regions(*material);
            
      rank1 Sigma_t = matprop.get_Sigma_t(*material);
      rank3 Sigma_sn = matprop.get_Sigma_sn(*material);
      rank1 Sigma_f = matprop.get_Sigma_f(*material);
      rank1 chi = matprop.get_chi(*material);
      rank1 nu = matprop.get_nu(*material);
      rank1 src_data = matprop.get_iso_src(*material);
     
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
      {
        for (unsigned int n = 0; n < M; n++)
        {        
          add_matrix_form(new VolumetricStreamingAndReactionsMF(regions, n, gto, G, Sigma_t[gto]));
          
          if (!Sigma_sn.empty())
            add_vector_form(new VolumetricScatteringSourceVF(odata, regions, n, gto, G, Sigma_sn));
          if (chi_nnz[gto])
            add_vector_form(new VolumetricFissionSourceVF(odata, regions, n, gto, G, chi[gto], nu, Sigma_f));
         
          if (src_data[gto] > 0) 
            add_vector_form(new VolumetricExternalSourceVF(regions, n, gto, G, src_data[gto], true));
          
          add_matrix_form_DG(new InterfaceStreamingMF(n, gto, G));
          add_matrix_form_surf(new BoundaryStreamingMF(n, gto, G));
          
          if (!inflow_boundaries.empty())
            add_vector_form_surf(new BoundaryStreamingVF(n, gto, G, inflow_boundaries));
          if (!reflective_boundaries.empty())
            add_vector_form_surf(new SpecularReflectionVF(odata, n, gto, G, reflective_boundaries));
        }
      }
    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * ( -static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[i], v->dy[i]) +
                                    Sigma_t * v->val[i] );
  
  //std::cout << "VolumetricStreamingAndReactionsMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::VolumetricStreamingAndReactionsMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                  Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0] * (v->dx[0]+v->dy[0]+v->val[0]);
}
/*
template<typename Real>
Real SNWeakForm::VolumetricStreamingAndReactionsMF::b(Real x, Real y) const
{
  return Real(1);
}
*/

template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *moment_values_at_quad_pts = new Real [n];

  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      
      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
          odata.ordinates_to_moment(l, m, gfrom, G, u_ext, n, moment_values_at_quad_pts);
          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
            
          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real group_source_moment = moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
            
            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt] * e->x[quad_pt];
            else
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt];
          }
        }
      }
      
      deg_l_result *= (2*l + 1) / (4*M_PI) * Sigma_sn[l][gto][gfrom];
      result += deg_l_result;
    }
  }
    
  delete [] moment_values_at_quad_pts;
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return result;
}


template<typename Real>
Real SNWeakForm::VolumetricFissionSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
    
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (int quad_pt = 0; quad_pt < n; quad_pt++)
    {
      // Angular integration to get scalar flux at quadrature point.
      Real group_scalar_flux(0.0);
      for (int dir = 0; dir < odata.M; dir++)
        group_scalar_flux += odata.pw[dir] * u_ext[ag.pos(dir,gfrom)]->val[quad_pt];
      //group_scalar_flux *= 2 / (4*M_PI) * 2*M_PI;
      group_scalar_flux *= 2;
      
      // Compute fission source with element-wise constant nuSigma_f.
      group_scalar_flux *= nu[gfrom] * Sigma_f[gfrom];
      
      // Contribution to spatial integral.
      if (geom_type == HERMES_AXISYM_X)
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt] * e->y[quad_pt];
      else if (geom_type == HERMES_AXISYM_Y)
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt] * e->x[quad_pt];
      else
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt];
    }
  }
  
  return result * chi_to;
}


template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  for (int i = 0; i < n; i++)
    result += wt[i] * Q * v->val[i]; //F(e->x[i], e->y[i]) * v->val[i];
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result;
}

/*
double SNWeakForm::VolumetricExternalSourceVF::Q(double x, double y) const
{
  //return 1 + (3*M_PI*(1 + x)*(1 + y)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/20. + (M_PI*std::pow(1 + y,2)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/10. + std::sin((M_PI*(1 + x)*std::pow(1 + y,2))/8.);
 
  return 6;
}

Ord SNWeakForm::VolumetricExternalSourceVF::Q(Ord x, Ord y) const
{
  return Ord(0);
}
*/

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
  {
    double x = e->x[i], y = e->y[i];
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);
    
    if (a_dot_n >= 0)
      result += wt[i] * u->val[i] * a_dot_n * v->val[i];
  }
  
  //std::cout << "BoundaryStreamingMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::BoundaryStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::InterfaceStreamingMF::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
                                               Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.0;

  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);
    double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * static_cast<SNWeakForm*>(wf)->upwind_flux(0.0, u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * static_cast<SNWeakForm*>(wf)->upwind_flux(u->val[i], 0.0, a_dot_n) * jump_v;
  }
  
  //std::cout << "InterfaceStreamingMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::InterfaceStreamingMF::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
                                          Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{ 
  Ord jump_v = (v->fn_central == NULL ? v->val_neighbor[0] : v->val[0]);
  if(u->fn_central == NULL)
    return static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), u->val_neighbor[0]) * jump_v;
  else
    return static_cast<SNWeakForm*>(wf)->upwind_flux(u->val[0], Ord(0)) * jump_v;
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double x = e->x[quad_pt], y = e->y[quad_pt];
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
    {
      double boundary_data = influx<double>(x,y);
      result += -wt[quad_pt] * boundary_data * a_dot_n * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::BoundaryStreamingVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  //std::cout << "BoundaryStreamingVF :: " << v->val[0] * influx<Ord>(e->x[0], e->y[0]) << std::endl;
  return v->val[0] * influx<Ord>(e->x[0], e->y[0]);
}

double SNWeakForm::SpecularReflectionVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double x = e->x[quad_pt], y = e->y[quad_pt];
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
   
    if (a_dot_n < 0)
    {
      double boundary_data = 0.0;
      
      if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        boundary_data = u_ext[ag.pos(odata.reflections_about_x[direction],g)]->val[quad_pt];
      else if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
        boundary_data = u_ext[ag.pos(odata.reflections_about_y[direction],g)]->val[quad_pt];
      else
        Hermes::Mixins::Loggable::Static::warn("Only horizontal or vertical boundaries are currently supported for specular reflection b. c.");

      result += -wt[quad_pt] * boundary_data * a_dot_n * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::SpecularReflectionVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  int refl_ord = std::max(u_ext[ag.pos(odata.reflections_about_x[direction],g)]->val[0].get_order(),u_ext[ag.pos(odata.reflections_about_y[direction],g)]->val[0].get_order());
  
  //std::cout << "SpecularReflectionVF :: " << static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), Ord(refl_ord), Ord(1)) * v->val[0] << std::endl;
  
  return static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), Ord(refl_ord)) * v->val[0];
}

double SNWeakForm::calculate_a_dot_v(int n, double vx, double vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return odata.xi[n]*vx + odata.eta[n]*vy;
}

double SNWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord SNWeakForm::upwind_flux(Ord u_cent, Ord u_neib) const
{
  return (u_cent + u_neib);
}
