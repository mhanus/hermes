////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion;

int get_num_of_neg(MeshFunction *sln);
void report_num_dof(const std::string& msg, const Hermes::vector<Space *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);

// Jacobian matrix (same as stiffness matrix since projections are linear).
class H1AxisymProjectionJacobian : public WeakForm::MatrixFormVol
{
public:
  H1AxisymProjectionJacobian(int i) : WeakForm::MatrixFormVol(i, i, HERMES_ANY, HERMES_SYM) {};

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                Geom<double> *e, ExtData<scalar> *ext) const
  {
    return h1_axisym_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
  
  template<typename Real, typename Scalar>
  static Scalar h1_axisym_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }
};

// Residual.
class H1AxisymProjectionResidual : public WeakForm::VectorFormVol
{
public:
  H1AxisymProjectionResidual(int i, MeshFunction* ext) : WeakForm::VectorFormVol(i)
  {
    this->ext.push_back(ext);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                Geom<double> *e, ExtData<scalar> *ext) const
  {
    return h1_axisym_projection_liform<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_liform<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

private:
  template<typename Real, typename Scalar>
  Scalar h1_axisym_projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                     Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * ( (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i] 
                                  + (u_ext[this->i]->dx[i]  - ext->fn[0]->dx[i])  * v->dx[i] 
                                  + (u_ext[this->i]->dy[i]  - ext->fn[0]->dy[i])  * v->dy[i]  );
    return result;
  }
};

// Matrix forms for error calculation.
class ErrorForm : public Adapt::MatrixFormVolError
{
public:
  ErrorForm(ProjNormType type) : Adapt::MatrixFormVolError(type) {};

  /// Error bilinear form.
  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                        Func<scalar> *u, Func<scalar> *v, Geom<double> *e,
                        ExtData<scalar> *ext) const;

  /// Error bilinear form to estimate order of a function.
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                  ExtData<Ord> *ext) const;

private:
  template<typename Real, typename Scalar>
  static Scalar l2_error_form_axisym(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u,
                                     Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * conj(v->val[i]));
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar h1_error_form_axisym(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u,
                                     Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i])
                                 + u->dy[i] * conj(v->dy[i]));
    return result;
  }  
};