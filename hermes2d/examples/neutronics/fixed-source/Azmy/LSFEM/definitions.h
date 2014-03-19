#include "hermes2d.h"
#include "../../../utils.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Mixins;
using namespace Neutronics;
using namespace Neutronics::SN;

//#define VER1
//#define VER2
#define VER3

class SNWeakForm : public WeakForm<double>
{
public:
  SNWeakForm(unsigned int N, const MaterialProperties::MaterialPropertyMaps& matprop,
             const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries = Hermes::vector<std::string>(),
             const char* out_tensor = "");
  
  const SupportClasses::OrdinatesData& get_ordinates_data() const { return odata; }
  
  virtual WeakForm<double>* clone() const { return new SNWeakForm(*this); }
  
private:
  class GenericForm
  {
    protected:
      GeomType geom_type;
      unsigned int G;
      unsigned int direction;
      unsigned int direction_from;
      SupportClasses::AngleGroupFlattener ag;
      
      GenericForm(unsigned int direction, unsigned int G, GeomType geom_type = HERMES_PLANAR)
        : geom_type(geom_type), direction(direction), direction_from(direction), G(G), ag(G)
      {};
      GenericForm(unsigned int direction, unsigned int direction_from, unsigned int G, GeomType geom_type = HERMES_PLANAR)
				: geom_type(geom_type), direction(direction), direction_from(direction_from), G(G), ag(G)
			{};
  };
    
  class VolumetricStreamingAndReactionsMF : protected GenericForm, public MatrixFormVol<double>
  {
    double Sigma_t;
    const SupportClasses::OrdinatesData& odata;
    
  public:
    VolumetricStreamingAndReactionsMF(const SupportClasses::OrdinatesData& odata, unsigned int n, unsigned int g, unsigned int G, double Sigma_t)
      : GenericForm(n, G), MatrixFormVol<double>(ag.pos(n,g), ag.pos(n,g)), 
        Sigma_t(Sigma_t), odata(odata)
    {
      this->setSymFlag(HERMES_SYM);
    };
    
    VolumetricStreamingAndReactionsMF(const SupportClasses::OrdinatesData& odata, const Hermes::vector<std::string>& areas,
                                      unsigned int n, unsigned int g, unsigned int G, double Sigma_t) 
      : GenericForm(n, G), MatrixFormVol<double>(ag.pos(n,g), ag.pos(n,g)), 
        Sigma_t(Sigma_t), odata(odata)
    {
      this->setSymFlag(HERMES_SYM);
      set_areas(areas);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const
    {
      return new SNWeakForm::VolumetricStreamingAndReactionsMF(*this);
    }
    
    //template<typename Real>
    //Real b(Real x, Real y) const;
  };
  
  class VolumetricScatteringSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    const SupportClasses::OrdinatesData& odata;
    rank3 Sigma_sn;
    double Sigma_t;
    unsigned int L;
    unsigned int gto;
    
  public:
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 unsigned int n, unsigned int g, unsigned int G,
                                 const rank3& Sigma_sn, double Sigma_t) 
      : GenericForm(n, G), VectorFormVol<double>(ag.pos(n,g)), 
        odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1), gto(g)
    {};
    
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 const Hermes::vector<std::string>& areas, 
                                 unsigned int n, unsigned int g, unsigned int G,
                                 const rank3& Sigma_sn, double Sigma_t) 
      : GenericForm(n, G), VectorFormVol<double>(ag.pos(n,g)), 
        odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1), gto(g)
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
    
    //template<typename Real>
    //Real b(Real x, Real y) const;
  };
  
  class VolumetricScatteringMF : protected GenericForm, public MatrixFormVol<double>
	{
	    const SupportClasses::OrdinatesData& odata;
	    rank3 Sigma_sn;
	    double Sigma_t;
	    unsigned int L;
	    unsigned int gto, gfrom;

	public:
	    VolumetricScatteringMF(const SupportClasses::OrdinatesData& odata,
              unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom, unsigned int G,
              const rank3& Sigma_sn, double Sigma_t)
		: GenericForm(m, n, G), MatrixFormVol<double>(ag.pos(m,gto), ag.pos(n,gfrom)),
		  odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1), gto(gto), gfrom(gfrom)
	  {
#if defined(VER2) or defined(VER3)
	    	setSymFlag(HERMES_SYM);
#endif
	  };

	    VolumetricScatteringMF(const SupportClasses::OrdinatesData& odata,
			  const Hermes::vector<std::string>& areas,
	                unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom, unsigned int G,
	                const rank3& Sigma_sn, double Sigma_t)
	  		: GenericForm(m, n, G), MatrixFormVol<double>(ag.pos(m,gto), ag.pos(n,gfrom)),
	  		  odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1), gto(gto), gfrom(gfrom)
	  {
#if defined(VER2) or defined(VER3)
	    	setSymFlag(HERMES_SYM);
#endif
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
    const SupportClasses::OrdinatesData& odata;
    rank3 Sigma_sn;
    double Sigma_t;
    unsigned int L;
    double Q;
    unsigned int gto;
    
  public:
    VolumetricExternalSourceVF(const SupportClasses::OrdinatesData& odata,
                               unsigned int n, unsigned int g, unsigned int G, 
                               const rank3& Sigma_sn, double Sigma_t, 
                               double Q, bool isotropic = false) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g)), 
        odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1),
        Q(isotropic ? Q/(4*M_PI) : Q), gto(g)
    {
    };
    
    VolumetricExternalSourceVF(const SupportClasses::OrdinatesData& odata,
                               const Hermes::vector<std::string>& areas,
                               unsigned int n, unsigned int g, unsigned int G, 
                               const rank3& Sigma_sn, double Sigma_t,
                               double Q, bool isotropic = false) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g)),
        odata(odata), Sigma_sn(Sigma_sn), Sigma_t(Sigma_t), L(Sigma_sn.size()-1),
        Q(isotropic ? Q/(4*M_PI) : Q), gto(g)
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
  public:
    BoundaryStreamingMF(const SupportClasses::OrdinatesData& odata, unsigned int n, unsigned int g, unsigned int G)
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(n,g)), odata(odata)
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
                           unsigned int n, unsigned int g, unsigned int G,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(odata.reflections_about_x[n],g)), odata(odata)
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
                           unsigned int n, unsigned int g, unsigned int G,
                           const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(odata.reflections_about_y[n],g)), odata(odata)
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
    BoundaryStreamingVF(unsigned int n, unsigned int g, unsigned int G, 
                        const Hermes::vector<std::string>& inflow_boundaries) 
      : GenericForm(n, G), VectorFormSurf<double>(ag.pos(n,g))
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
  
  class SpecularReflectionVF : protected GenericForm, public VectorFormSurf<double>
  {
    const SupportClasses::OrdinatesData& odata;
    unsigned int g;
    
  public:
    SpecularReflectionVF(const SupportClasses::OrdinatesData& odata,
                         unsigned int n, unsigned int g, unsigned int G, 
                         const Hermes::vector<std::string>& reflective_boundaries) 
      : GenericForm(n, G), VectorFormSurf<double>(ag.pos(n,g)), 
        odata(odata), g(g)
    {
      set_areas(reflective_boundaries);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const
    {
      return new SNWeakForm::SpecularReflectionVF(*this);
    }
  };
    
  double calculate_a_dot_v(int n, double vx, double vy) const;
  Ord calculate_a_dot_v(int n, Ord vx, Ord vy) const;
    
  unsigned int N, M, G;
  SupportClasses::OrdinatesData odata;
};

class IsotropicScatteringAndFissionMatrixForms : public WeakForm<double>
{
public:
  IsotropicScatteringAndFissionMatrixForms(const MaterialProperties::MaterialPropertyMaps& matprop, const char* out_tensor);
  virtual WeakForm<double>* clone() const { return new IsotropicScatteringAndFissionMatrixForms(*this); }
};
