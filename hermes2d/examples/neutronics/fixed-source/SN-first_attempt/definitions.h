#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(const std::string& inflow_boundary, Mesh* mesh);

private:
  class CustomMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVol(int i, int j) : MatrixFormVol<double>(i, j)
    {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone()
    {
      return new CustomWeakForm::CustomMatrixFormVol(*this);
    }
    
    template<typename Real, typename Scalar>
    Scalar b(Real x, Real y) const;
  };

  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i) : VectorFormVol<double>(i)
    {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }    
    
    VectorFormVol<double>* clone()
    {
      return new CustomWeakForm::CustomVectorFormVol(*this);
    }

    double F(double x, double y) const;
    Ord F(Ord x, Ord y) const;
  };

  class CustomMatrixFormSurface : public MatrixFormSurf<double>
  {
  public:
    CustomMatrixFormSurface(int i, int j) : MatrixFormSurf<double>(i, j, HERMES_ANY)
    {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }    

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
    
    MatrixFormSurf<double>* clone()
    {
      return new CustomWeakForm::CustomMatrixFormSurface(*this);
    }
  };

  class CustomMatrixFormInterface : public MatrixFormSurf<double>
  {
  public:
    CustomMatrixFormInterface(int i, int j) : MatrixFormSurf<double>(i, j, H2D_DG_INNER_EDGE) 
    {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormSurf<double>* clone()
    {
      return new CustomWeakForm::CustomMatrixFormInterface(*this);
    }
  };

  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface(int i, const std::string& inflow_boundary) : VectorFormSurf<double>(i, HERMES_ANY),
      inflow_boundary(inflow_boundary) 
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    VectorFormSurf<double>* clone()
    {
      return new CustomWeakForm::CustomVectorFormSurface(*this);
    }

    template<typename Real, typename Scalar>
    Scalar g(const std::string& bdy_marker, Real x, Real y) const;
    
    // Member.
    std::string inflow_boundary;
  };
  
  double calculate_a_dot_v(double x, double y, double vx, double vy) const;

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const;

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;

  Mesh* mesh;
};
