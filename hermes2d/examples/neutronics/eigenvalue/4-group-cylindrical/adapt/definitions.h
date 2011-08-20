////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "../weak_formulation.h"
using Hermes::Ord;

int get_num_of_neg(MeshFunction<double> *sln);
void report_num_dof(const std::string& msg, const Hermes::vector<Space<double> *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);

// Jacobian matrix (same as stiffness matrix since projections are linear).
template <typename Scalar>
class H1AxisymProjectionJacobian : public MatrixFormVol<Scalar>
{
public:
  H1AxisymProjectionJacobian(int i) : MatrixFormVol<Scalar>(i, i, Hermes::HERMES_ANY, HERMES_SYM) {};

  Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
               Geom<double> *e, ExtData<Scalar> *ext) const
  {
    return h1_axisym_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
  
  template<typename TestFunctionDomain, typename SolFunctionDomain>
  static SolFunctionDomain h1_axisym_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
                                                       Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
  {
    _F_
    SolFunctionDomain result(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }
};

// Residual.
template <typename Scalar>
class H1AxisymProjectionResidual : public VectorFormVol<Scalar>
{
public:
  H1AxisymProjectionResidual(int i, MeshFunction<Scalar>* ext) : VectorFormVol<Scalar>(i)
  {
    this->ext.push_back(ext);
  }

  Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                Geom<double> *e, ExtData<Scalar> *ext) const
  {
    return h1_axisym_projection_liform<double, Scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_liform<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

private:
  template<typename TestFunctionDomain, typename SolFunctionDomain>
  SolFunctionDomain h1_axisym_projection_liform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
                                                Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
  {
    _F_
    SolFunctionDomain result(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * ( (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i] 
                                  + (u_ext[this->i]->dx[i]  - ext->fn[0]->dx[i])  * v->dx[i] 
                                  + (u_ext[this->i]->dy[i]  - ext->fn[0]->dy[i])  * v->dy[i]  );
    return result;
  }
};

// Matrix forms for error calculation.
template <typename Scalar>
class ErrorForm : public Adapt<Scalar>::MatrixFormVolError
{
public:
  ErrorForm(ProjNormType type) : Adapt<Scalar>::MatrixFormVolError(type) {};

  /// Error bilinear form.
  virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
                        Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
                        ExtData<Scalar> *ext) const;

  /// Error bilinear form to estimate order of a function.
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                  ExtData<Ord> *ext) const;

private:
  template<typename TestFunctionDomain, typename SolFunctionDomain>
  static SolFunctionDomain l2_error_form_axisym(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
                                                Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
  {
    SolFunctionDomain result(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * v->val[i]);
    return result;
  }

  template<typename TestFunctionDomain, typename SolFunctionDomain>
  static SolFunctionDomain h1_error_form_axisym(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
                                                Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
  {
    SolFunctionDomain result(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }  
};

template<typename Scalar>
Scalar ErrorForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
                                ExtData<Scalar> *ext) const
{
  switch (this->projNormType)
  {
    case HERMES_L2_NORM:
      return ErrorForm<Scalar>::l2_error_form_axisym<double, Scalar>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return ErrorForm<Scalar>::h1_error_form_axisym<double, Scalar>(n, wt, u_ext, u, v, e, ext);
    default:
      error_function("Only the H1 and L2 norms are currently implemented.");
      return 0.0;
  }
}

template<typename Scalar>
Ord ErrorForm<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[],
                        Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                        ExtData<Ord> *ext) const
{
  switch (this->projNormType)
  {
    case HERMES_L2_NORM:
      return ErrorForm<Scalar>::l2_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return ErrorForm<Scalar>::h1_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    default:
      error_function("Only the H1 and L2 norms are currently implemented.");
      return Ord();
  }
}
