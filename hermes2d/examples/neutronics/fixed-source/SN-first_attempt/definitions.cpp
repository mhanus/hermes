#include "definitions.h"

CustomWeakForm::CustomWeakForm(const std::string& inflow_boundary, Mesh* mesh, int N) : WeakForm<double>(N), mesh(mesh), N(N)
{
  for (int n = 0; n < N; n++)
  {
    add_matrix_form(new CustomMatrixFormVol(n, n));
    add_vector_form(new CustomVectorFormVol(n));
    add_matrix_form_surf(new CustomMatrixFormSurface(n, n));
    add_matrix_form_surf(new CustomMatrixFormInterface(n, n));
    add_vector_form_surf(new CustomVectorFormSurface(n, inflow_boundary));
  }
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * ( -static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(direction, e->x[i], e->y[i], v->dx[i], v->dy[i]) + b<Real,Scalar>(e->x[i],e->y[i]) * v->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::b(Real x, Real y) const
{
  return Real(1);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * F(e->x[i], e->y[i]) * v->val[i];
  return result;
}

double CustomWeakForm::CustomVectorFormVol::F(double x, double y) const
{
  return 1 + (3*M_PI*(1 + x)*(1 + y)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/20. + (M_PI*std::pow(1 + y,2)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/10. + std::sin((M_PI*(1 + x)*std::pow(1 + y,2))/8.);
}

Ord CustomWeakForm::CustomVectorFormVol::F(Ord x, Ord y) const
{
  return Ord(0);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormSurface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                      Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
  {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = Real(static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(direction, x, y, e->nx[i], e->ny[i]));
    result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Scalar(0), a_dot_n) * v->val[i];
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                        Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = Scalar(0);

  for (int i = 0; i < n; i++) {
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(direction, e->x[i], e->y[i], e->nx[i], e->ny[i]);
    Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
    result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
  }
  return result;
}

double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i], y = e->y[i];
    double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(direction, x, y, e->nx[i], e->ny[i]);
    // Function values for Dirichlet boundary conditions.
    result += -wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0, g<double,double>(static_cast<CustomWeakForm*>(wf)->mesh->get_boundary_markers_conversion().get_user_marker(e->edge_marker).marker, x, y), a_dot_n) * v->val[i];
  }
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += -wt[i] * v->val[i];
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormSurface::g(const std::string& bdy_marker, Real x, Real y) const
{
  if(bdy_marker == inflow_boundary) return 1; else return 0;
}

double CustomWeakForm::calculate_a_dot_v(int n, double x, double y, double vx, double vy) const
{
  //double norm = std::max<double>(1e-12, std::sqrt(sqr(x) + sqr(y)));
  //return -y/norm*vx + x/norm*vy;
  return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
}

Ord CustomWeakForm::calculate_a_dot_v(int n, Ord x, Ord y, Ord vx, Ord vy) const
{
  return Ord(1);
}

double CustomWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}