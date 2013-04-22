#include "hermes2d.h"
#include "../../utils.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Mixins;
using namespace Neutronics;
using namespace Neutronics::SN;

class SourceFunction : public Hermes2DFunction<double>
{
  double c;
  double sigma_t;
  const SupportClasses::OrdinatesData& odata;
  unsigned int dir;
  
  public:   
    SourceFunction(double c, double sigma_t, const SupportClasses::OrdinatesData& odata)
      : c(c), sigma_t(sigma_t), dir(0), odata(odata)
    {};
    
    template<typename Real>
    Real Q(Real x, Real y) const;
    
    virtual double value(double x, double y) const { return Q(x,y); }
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const { return Q(x,y); }
    
    void set_direction(unsigned int n) { dir = n; }
};

class SNWeakForm : public WeakForm<double>
{
public:
  SNWeakForm(unsigned int N, double c, double sigma_t,
             const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries = Hermes::vector<std::string>(),
             const char* out_tensor = "");
  
  const SupportClasses::OrdinatesData& get_ordinates_data() const { return odata; }
  
  virtual WeakForm<double>* clone() const { return new SNWeakForm(*this); }
  
private:
  class GenericForm
  {
    protected:
      GeomType geom_type;
      unsigned int direction;
      
      GenericForm(unsigned int direction, GeomType geom_type = HERMES_PLANAR)
        : geom_type(geom_type), direction(direction)
      {};
  };
    
  class VolumetricStreamingAndReactionsMF : protected GenericForm, public MatrixFormVol<double>
  {
    double sigma_t;
    
  public:
    VolumetricStreamingAndReactionsMF(unsigned int n,
                                      double sigma_t) 
      : GenericForm(n), MatrixFormVol<double>(n, n), 
        sigma_t(sigma_t)
    {};
    
    VolumetricStreamingAndReactionsMF(const Hermes::vector<std::string>& areas, 
                                      unsigned int n,
                                      double sigma_t) 
      : GenericForm(n), MatrixFormVol<double>(n, n), 
        sigma_t(sigma_t)
    {
      set_areas(areas);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const
    {
      return new SNWeakForm::VolumetricStreamingAndReactionsMF(*this);
    }
  };
  
  class VolumetricScatteringSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    const SupportClasses::OrdinatesData& odata;
    double sigma_s;
    unsigned int L;
    
  public:
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 unsigned int n,
                                 double sigma_s) 
      : GenericForm(n), VectorFormVol<double>(n), 
        odata(odata), sigma_s(sigma_s), L(0)//, L(sigma_s.size()-1)
    {};
    
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 const Hermes::vector<std::string>& areas, 
                                 unsigned int n,
                                 double sigma_s) 
      : GenericForm(n), VectorFormVol<double>(n), 
        odata(odata), sigma_s(sigma_s), L(0)//L(sigma_s.size()-1)
    {
      set_areas(areas);
    };

    template<typename Real>
    Real vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Real> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
    {
      return vector_form<double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
    {
      return vector_form<Ord>(n, wt, u_ext, v, e, ext);
    }  

    VectorFormVol<double>* clone() const
    {
      return new SNWeakForm::VolumetricScatteringSourceVF(*this);
    }
  };

  class VolumetricExternalSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    SourceFunction Q;
    
  public:
    VolumetricExternalSourceVF(unsigned int n, 
                               const SourceFunction& Q) 
      : GenericForm(n), VectorFormVol<double>(n), Q(Q)
    {};
    
    VolumetricExternalSourceVF(const Hermes::vector<std::string>& areas,
                               unsigned int n,
                               const SourceFunction& Q) 
      : GenericForm(n), VectorFormVol<double>(n), Q(Q)
    {
      set_areas(areas);
    };

    template<typename Real>
    Real vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Real> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
    {
      return vector_form<double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
    {
      return vector_form<Ord>(n, wt, u_ext, v, e, ext);
    }    
    
    VectorFormVol<double>* clone() const
    {
      return new SNWeakForm::VolumetricExternalSourceVF(*this);
    }
  };
  
  class InterfaceStreamingMF : protected GenericForm, public MatrixFormDG<double>
  {
  public:
    InterfaceStreamingMF(unsigned int n) 
      : GenericForm(n), MatrixFormDG<double>(n, n)
    {};

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;
    
    MatrixFormDG<double>* clone() const
    {
      return new SNWeakForm::InterfaceStreamingMF(*this);
    }
  };
  
  class BoundaryStreamingMF : protected GenericForm, public MatrixFormSurf<double>
  { 
  public:
    BoundaryStreamingMF(unsigned int n) 
      : GenericForm(n), MatrixFormSurf<double>(n, n)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormSurf<double>* clone() const
    {
      return new SNWeakForm::BoundaryStreamingMF(*this);
    }
  };

  class SpecularReflectionMF_X : protected GenericForm, public MatrixFormSurf<double>
  { 
  public:
    SpecularReflectionMF_X(const SupportClasses::OrdinatesData& odata, 
                           unsigned int n,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n), MatrixFormSurf<double>(n, odata.reflections_about_x[n])
    {
      set_areas(reflective_boundaries);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormSurf<double>* clone() const
    {
      return new SNWeakForm::SpecularReflectionMF_X(*this);
    }
  };
  
  class SpecularReflectionMF_Y : protected GenericForm, public MatrixFormSurf<double>
  { 
  public:
    SpecularReflectionMF_Y(const SupportClasses::OrdinatesData& odata,
                           unsigned int n,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n), MatrixFormSurf<double>(n, odata.reflections_about_y[n])
    {
      set_areas(reflective_boundaries);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormSurf<double>* clone() const
    {
      return new SNWeakForm::SpecularReflectionMF_Y(*this);
    }
  };
  
  class BoundaryStreamingVF : protected GenericForm, public VectorFormSurf<double>
  {    
  public:
    BoundaryStreamingVF(unsigned int n, 
                        const Hermes::vector<std::string>& inflow_boundaries) 
      : GenericForm(n), VectorFormSurf<double>(n)
    {
      set_areas(inflow_boundaries);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const
    {
      return new SNWeakForm::BoundaryStreamingVF(*this);
    }

    template<typename Real>
    Real influx(Real x, Real y) const { return Real(0.0); } // vacuum inflow boundary
  };
    
  double calculate_a_dot_v(int n, double vx, double vy) const;
  
  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;
  Ord upwind_flux(Ord u_cent, Ord u_neib) const;
  
  unsigned int N, M;
  SupportClasses::OrdinatesData odata;
};

class IsotropicScatteringMatrixForm : public WeakForm<double>
{
  private:
    class ScatteringMF : public MatrixFormVol<double>
    {
      double sigma_s;
      
      public:
        ScatteringMF(double c, double sigma_t) : MatrixFormVol< double >(0,0), sigma_s(c*sigma_t) {};
        
        template<typename Real>
        Real matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, Func<Real> **ext  ) const;
        
        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const 
        { 
          return matrix_form<double> (n, wt, u_ext, u, v, e, ext);
        }
        
        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
        { 
          return matrix_form<Ord> (n, wt, u_ext, u, v, e, ext);
        }
    
        MatrixFormVol<double>* clone() const
        {
          return new IsotropicScatteringMatrixForm::ScatteringMF(*this);
        }
    };
    
  public:
    IsotropicScatteringMatrixForm(double c, double sigma_t);
    virtual WeakForm<double>* clone() const { return new IsotropicScatteringMatrixForm(*this); }
};