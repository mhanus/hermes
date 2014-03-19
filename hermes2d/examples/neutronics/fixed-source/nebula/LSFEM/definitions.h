#include "hermes2d.h"
#include "../../../utils.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Mixins;
using namespace Neutronics;
using namespace Neutronics::SN;

class ExtinctionFunction : public Hermes2DFunction<double>
{
  double chi0;
  
  public:   
    ExtinctionFunction(double base_extinction_coeff) 
      : chi0(base_extinction_coeff)
    {};
    
    virtual double value(double x, double y) const;
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const;
};

class ScatteringFunction : public Hermes2DFunction<double>
{
  double scattering_ratio;
  ExtinctionFunction chi;
  
  public:   
    ScatteringFunction(double thermalization_coeff, const ExtinctionFunction& chi)
      : scattering_ratio(1-thermalization_coeff), chi(chi)
    {};
    
    virtual double value(double x, double y) const { return scattering_ratio*chi.value(x,y); }
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const { return scattering_ratio*chi.value(x,y); }
};

class SourceFunction : public Hermes2DFunction<double>
{
  double eps;
  ExtinctionFunction chi;
  
  public:   
    SourceFunction(double thermalization_coeff, const ExtinctionFunction& chi)
      : eps(thermalization_coeff), chi(chi)
    {};
    
    virtual double value(double x, double y) const { return /*1.0*/eps*chi.value(x,y)/*/(4*M_PI)*/; }
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const { return /*Hermes::Ord(1.0)*/ eps*chi.value(x,y); }
};

class SNWeakForm : public WeakForm<double>
{
public:
  SNWeakForm(unsigned int N, double extinction, double thermalization,
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
      unsigned int direction_from;
      
      GenericForm(unsigned int direction, GeomType geom_type = HERMES_PLANAR)
				: geom_type(geom_type), direction(direction), direction_from(direction)
			{};
      GenericForm(unsigned int direction, unsigned int direction_from, GeomType geom_type = HERMES_PLANAR)
        : geom_type(geom_type), direction(direction), direction_from(direction_from)
      {};
  };
    
  class VolumetricStreamingAndReactionsMF : protected GenericForm, public MatrixFormVol<double>
  {
    ExtinctionFunction Sigma_t;
    const SupportClasses::OrdinatesData& odata;
    
  public:
    VolumetricStreamingAndReactionsMF(const SupportClasses::OrdinatesData& odata, unsigned int n,
                                      const ExtinctionFunction& Sigma_t) 
      : GenericForm(n), MatrixFormVol<double>(n, n),
        Sigma_t(Sigma_t), odata(odata)
    {
    	setSymFlag(HERMES_SYM);
    };
    
    VolumetricStreamingAndReactionsMF(const SupportClasses::OrdinatesData& odata, const Hermes::vector<std::string>& areas,
    																	unsigned int n,
                                      const ExtinctionFunction& Sigma_t) 
      : GenericForm(n), MatrixFormVol<double>(n, n),
        Sigma_t(Sigma_t), odata(odata)
    {
    	setSymFlag(HERMES_SYM);
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
    ScatteringFunction Sigma_s;
    ExtinctionFunction Sigma_t;
    unsigned int L;
    
  public:
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 unsigned int n,
                                 const ScatteringFunction& Sigma_s,
                                 const ExtinctionFunction& Sigma_t)
      : GenericForm(n), VectorFormVol<double>(n), 
        odata(odata), Sigma_s(Sigma_s), Sigma_t(Sigma_t), L(0)//, L(Sigma_s.size()-1)
    {};
    
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 const Hermes::vector<std::string>& areas, 
                                 unsigned int n,
                                 const ScatteringFunction& Sigma_s,
                                 const ExtinctionFunction& Sigma_t)
      : GenericForm(n), VectorFormVol<double>(n), 
        odata(odata), Sigma_s(Sigma_s), Sigma_t(Sigma_t), L(0)//L(Sigma_s.size()-1)
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

  class VolumetricScatteringMF : protected GenericForm, public MatrixFormVol<double>
	{
	    const SupportClasses::OrdinatesData& odata;
	    ScatteringFunction Sigma_s;
	    ExtinctionFunction Sigma_t;
	    unsigned int L;

	public:
	    VolumetricScatteringMF(const SupportClasses::OrdinatesData& odata,
              unsigned int m, unsigned int n,
              const ScatteringFunction& Sigma_s, const ExtinctionFunction& Sigma_t)
				: GenericForm(m,n), MatrixFormVol<double>(m, n),
					odata(odata), Sigma_s(Sigma_s), Sigma_t(Sigma_t), L(0)
			{
	    	setSymFlag(HERMES_SYM);
			};

	    VolumetricScatteringMF(const SupportClasses::OrdinatesData& odata,
			  const Hermes::vector<std::string>& areas,
	                unsigned int m, unsigned int n,
	                const ScatteringFunction& Sigma_s, const ExtinctionFunction& Sigma_t)
	  		: GenericForm(m,n), MatrixFormVol<double>(m, n),
	  		  odata(odata), Sigma_s(Sigma_s), Sigma_t(Sigma_t), L(0)
			{
	    	setSymFlag(HERMES_SYM);
	    	set_areas(areas);
			};

		template<typename Real>
		Real matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Real> **ext) const;

	  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const
	  {
		  return matrix_form<double>(n, wt, u_ext, u, v, e, ext);
	  }

	  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
	  {
		  return matrix_form<Ord>(n, wt, u_ext, u, v, e, ext);
	  }

	  MatrixFormVol<double>* clone() const
	  {
	  	return new SNWeakForm::VolumetricScatteringMF(*this);
	  }

	  //template<typename Real>
	  //Real b(Real x, Real y) const;
	};

  class VolumetricExternalSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    SourceFunction Q;
    ExtinctionFunction Sigma_t;
    ScatteringFunction Sigma_s;
    SupportClasses::OrdinatesData odata;
    
  public:
    VolumetricExternalSourceVF(const SupportClasses::OrdinatesData& odata,
    													 unsigned int n,
                               const SourceFunction& Q,
                               const ScatteringFunction& Sigma_s,
                               const ExtinctionFunction& Sigma_t)
      : GenericForm(n), VectorFormVol<double>(n), Q(Q), Sigma_t(Sigma_t), Sigma_s(Sigma_s), odata(odata)
    {};
    
    VolumetricExternalSourceVF(const SupportClasses::OrdinatesData& odata,
															 const Hermes::vector<std::string>& areas,
                               unsigned int n,
                               const SourceFunction& Q,
                               const ScatteringFunction& Sigma_s,
                               const ExtinctionFunction& Sigma_t)
      : GenericForm(n), VectorFormVol<double>(n), Q(Q), Sigma_t(Sigma_t), Sigma_s(Sigma_s), odata(odata)
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
  
  class BoundaryStreamingMF : protected GenericForm, public MatrixFormSurf<double>
  { 
  	const SupportClasses::OrdinatesData& odata;
  	double alpha;

  public:
    BoundaryStreamingMF(const SupportClasses::OrdinatesData& odata, unsigned int n, double alpha=1.0)
      : GenericForm(n), MatrixFormSurf<double>(n, n), odata(odata), alpha(alpha)
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
  	const SupportClasses::OrdinatesData& odata;
  public:
    SpecularReflectionMF_X(const SupportClasses::OrdinatesData& odata, 
                           unsigned int n,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n), MatrixFormSurf<double>(n, odata.reflections_about_x[n]), odata(odata)
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
  	const SupportClasses::OrdinatesData& odata;
  public:
    SpecularReflectionMF_Y(const SupportClasses::OrdinatesData& odata,
                           unsigned int n,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n), MatrixFormSurf<double>(n, odata.reflections_about_y[n]), odata(odata)
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
  	double alpha;

  public:
    BoundaryStreamingVF(unsigned int n, 
                        const Hermes::vector<std::string>& inflow_boundaries,
                        double alpha = 1.0)
      : GenericForm(n), VectorFormSurf<double>(n), alpha(alpha)
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
  Ord calculate_a_dot_v(int n, Ord vx, Ord vy) const;
  
  unsigned int N, M;
  SupportClasses::OrdinatesData odata;
};

class IsotropicScatteringMatrixForm : public WeakForm<double>
{
  private:
    class ScatteringMF : public MatrixFormVol<double>
    {
      ScatteringFunction Sigma_s;
      
      public:
        ScatteringMF(const ScatteringFunction& Sigma_s) : MatrixFormVol< double >(0,0), Sigma_s(Sigma_s) {};
        
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
    IsotropicScatteringMatrixForm(double extinction, double thermalization);
    virtual WeakForm<double>* clone() const { return new IsotropicScatteringMatrixForm(*this); }
};
