#include "definitions.h"
#include <cstdlib>

IsotropicScatteringMatrixForm::IsotropicScatteringMatrixForm(double sigma_t, double c) : WeakForm<double>(1)
{
  add_matrix_form(new IsotropicScatteringMatrixForm::ScatteringMF(sigma_t, c));
  add_vector_form(new Diffusion::WeakFormParts::ExternalSources::LinearForm(0, -1.0));
}

template<typename Real>
Real SourceFunction::Q(Real x, Real y) const
{
  double xi = odata.xi[dir];
  double eta = odata.eta[dir];
  double kappa = c*sigma_t/4;
  
  if (xi < 0 && eta < 0)
    return eta*(x-1) + xi*(y-1) - kappa + (1-x)*(1-y)*sigma_t;
  else if (xi > 0 && eta < 0)
    return xi*(1-y) - eta*x - kappa + (1-y)*x*sigma_t;
  else if (xi < 0 && eta > 0)
    return eta*(1-x) - xi*y - kappa + (1-x)*y*sigma_t;
  else if (xi > 0 && eta > 0)
    return eta*x + xi*y - kappa + x*y*sigma_t;
  else
    return Real(kappa);
}

template<typename Real>
Real IsotropicScatteringMatrixForm::ScatteringMF::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                               Func<Real> *v, Geom<Real> *e, Func<Real> **ext  ) const
{
  Real result(0.0);
  for (int i = 0; i < n; i++)
    result += wt[i] * sigma_s * u->val[i] * v->val[i];
  return result;
}

SNWeakForm::SNWeakForm(unsigned int N, double c, double sigma_t,
                       const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries,
                       const char* out_tensor) 
  : WeakForm<double>(N*(N+2)/2), N(N), M(N*(N+2)/2), odata(SupportClasses::OrdinatesData(N, "lgvalues.txt"))
{  
    std::cout << odata;
    
    bool assemble_all = strlen(out_tensor) == 0;    
    bool assemble_L = !strcmp(out_tensor, "L");
    bool assemble_Q = !strcmp(out_tensor, "Q");

    if (!assemble_all && !strcmp(out_tensor, "LQ"))
    {
      assemble_L = true;
      assemble_Q = true;
    }
    
    SourceFunction Q(c, sigma_t, odata);
    
    for (unsigned int n = 0; n < M; n++)
    {
      if (assemble_all || assemble_L)
      {
        add_matrix_form_DG(new InterfaceStreamingMF(n));
        if (!reflective_boundaries.empty())
        {
          add_matrix_form_surf(new SpecularReflectionMF_X(odata, n, reflective_boundaries));
          add_matrix_form_surf(new SpecularReflectionMF_Y(odata, n, reflective_boundaries));
        }
        add_matrix_form_surf(new BoundaryStreamingMF(n));
        
        add_matrix_form(new VolumetricStreamingAndReactionsMF(n, sigma_t));
      }
      
      if (assemble_all || assemble_Q) 
      {
        Q.set_direction(n);
        add_vector_form(new VolumetricExternalSourceVF(n, Q));
        add_vector_form_surf(new BoundaryStreamingVF(n, inflow_boundaries));
      }
      
      if (assemble_all)
        add_vector_form(new VolumetricScatteringSourceVF(odata, n, c*sigma_t));
      
    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * ( -static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[i], v->dy[i]) +
                                    sigma_t * v->val[i] );
  
  //std::cout << "VolumetricStreamingAndReactionsMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::VolumetricStreamingAndReactionsMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                  Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0] * (v->dx[0]+v->dy[0]+v->val[0]+sigma_t);
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

  odata.ordinates_to_moment(0, 0, 0, 1, u_ext, n, moment_values_at_quad_pts);
  SupportClasses::SphericalHarmonic Rlm(0, 0);
  double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
    
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    Real source = sigma_s * moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
    
    // Contribution to spatial integral.
    if (geom_type == HERMES_AXISYM_X)
      result += source * wt[quad_pt] * v->val[quad_pt] * e->y[quad_pt];
    else if (geom_type == HERMES_AXISYM_Y)
      result += source * wt[quad_pt] * v->val[quad_pt] * e->x[quad_pt];
    else
      result += source * wt[quad_pt] * v->val[quad_pt];
  }
  
  result *= 1 / (4*M_PI);
    
  delete [] moment_values_at_quad_pts;
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return result;
}

template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  for (int i = 0; i < n; i++)
    result += wt[i] * Q.value(e->x[i], e->y[i]) * v->val[i]; //F(e->x[i], e->y[i]) * v->val[i];
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result;
}

double SNWeakForm::SpecularReflectionMF_X::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
  
    if (a_dot_n < 0)
      if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        result += wt[quad_pt] * u->val[quad_pt] * a_dot_n * v->val[quad_pt];
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_X::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::SpecularReflectionMF_Y::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
      if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
        result += wt[quad_pt] * u->val[quad_pt] * a_dot_n * v->val[quad_pt];
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_Y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::InterfaceStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                        Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;

  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);
    double jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
    result += wt[i] * static_cast<SNWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
  }
  
  //std::cout << "InterfaceStreamingMF :: " << result << std::endl;
  
  return result;
}

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int i = 0; i < n; i++)
  {
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

Ord SNWeakForm::InterfaceStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return static_cast<SNWeakForm*>(wf)->upwind_flux(u->get_val_central(0), u->get_val_neighbor(0)) * (v->get_val_central(0) - v->get_val_neighbor(0));
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
    {
      double boundary_data = influx<double>(e->x[quad_pt],e->y[quad_pt]);
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

double SNWeakForm::calculate_a_dot_v(int n, double vx, double vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return odata.xi[n]*vx + odata.eta[n]*vy;
}

double SNWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  if (fabs(a_dot_n) < 1e-14)
    return 0;
  
  return a_dot_n * (a_dot_n > 0 ? u_cent : u_neib);
}

Ord SNWeakForm::upwind_flux(Ord u_cent, Ord u_neib) const
{
  return (u_cent + u_neib);
}
