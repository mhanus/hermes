#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Mixins;
using namespace Neutronics;
using namespace Neutronics::SN;

class SNWeakForm : public WeakForm<double>
{
public:
  SNWeakForm(unsigned int N, const MaterialProperties::MaterialPropertyMaps& matprop, const Hermes::vector<Solution<double>*>& iterates,
             const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries = Hermes::vector<std::string>());
  
  const SupportClasses::OrdinatesData& get_ordinates_data() const { return odata; }
  
private:
  class GenericForm
  {
    protected:
      GeomType geom_type;
      unsigned int G;
      unsigned int direction;
      SupportClasses::AngleGroupFlattener ag;
      
      GenericForm(unsigned int direction, unsigned int G, GeomType geom_type = HERMES_PLANAR)
        : geom_type(geom_type), direction(direction), G(G), ag(G)
      {};
  };
    
  class VolumetricStreamingAndReactionsMF : protected GenericForm, public MatrixFormVol<double>
  {
    double Sigma_t;
    
  public:
    VolumetricStreamingAndReactionsMF(unsigned int n, unsigned int g, unsigned int G, double Sigma_t) 
      : GenericForm(n, G), MatrixFormVol<double>(ag.pos(n,g), ag.pos(n,g)), 
        Sigma_t(Sigma_t)
    {};
    
    VolumetricStreamingAndReactionsMF(const Hermes::vector<std::string>& areas, 
                                      unsigned int n, unsigned int g, unsigned int G, double Sigma_t) 
      : GenericForm(n, G), MatrixFormVol<double>(ag.pos(n,g), ag.pos(n,g), areas), 
        Sigma_t(Sigma_t)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone()
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
    unsigned int L;
    unsigned int gto;
    
  public:
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 unsigned int n, unsigned int g, unsigned int G,
                                 const rank3& Sigma_sn, 
                                 const Hermes::vector<MeshFunction<double>*>& iterates) 
      : GenericForm(n, G), VectorFormVol<double>(ag.pos(n,g), HERMES_ANY, iterates), 
        odata(odata), Sigma_sn(Sigma_sn), L(Sigma_sn.size()-1), gto(g)
    {};
    
    VolumetricScatteringSourceVF(const SupportClasses::OrdinatesData& odata,
                                 const Hermes::vector<std::string>& areas, 
                                 unsigned int n, unsigned int g, unsigned int G,
                                 const rank3& Sigma_sn, 
                                 const Hermes::vector<MeshFunction<double>*>& iterates) 
      : GenericForm(n, G), VectorFormVol<double>(ag.pos(n,g), areas, iterates), 
        odata(odata), Sigma_sn(Sigma_sn), L(Sigma_sn.size()-1), gto(g)
    {};

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

    VectorFormVol<double>* clone()
    {
      return new SNWeakForm::VolumetricScatteringSourceVF(*this);
    }
    
    //template<typename Real>
    //Real b(Real x, Real y) const;
  };
  
  class VolumetricFissionSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    const SupportClasses::OrdinatesData& odata;
    unsigned int g;
    double chi_to;
    rank1 nu;
    rank1 Sigma_f;
    
  public:
    VolumetricFissionSourceVF(const SupportClasses::OrdinatesData& odata,
                              unsigned int n, unsigned int g, unsigned int G, 
                              double chi_to, const rank1& nu, const rank1& Sigma_f,
                              const Hermes::vector<MeshFunction<double>*>& iterates) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g), HERMES_ANY, iterates),
        odata(odata), g(g), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
    {};
    
    VolumetricFissionSourceVF(const SupportClasses::OrdinatesData& odata,
                              const Hermes::vector<std::string>& areas, 
                              unsigned int n, unsigned int g, unsigned int G,
                              double chi_to, const rank1& nu, const rank1& Sigma_f,
                              const Hermes::vector<MeshFunction<double>*>& iterates) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g), areas, iterates),
        odata(odata), g(g), chi_to(chi_to), nu(nu), Sigma_f(Sigma_f)
    {};

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

    VectorFormVol<double>* clone()
    {
      return new SNWeakForm::VolumetricFissionSourceVF(*this);
    }
    
    //template<typename Real>
    //Real b(Real x, Real y) const;
  };

  class VolumetricExternalSourceVF : protected GenericForm, public VectorFormVol<double>
  {
    double Q;
    
  public:
    VolumetricExternalSourceVF(unsigned int n, unsigned int g, unsigned int G, double Q, bool isotropic = false) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g)), Q(isotropic ? Q/(4*M_PI) : Q)
    {};
    
    VolumetricExternalSourceVF(const Hermes::vector<std::string>& areas,
                               unsigned int n, unsigned int g, unsigned int G, double Q, bool isotropic = false) 
      : GenericForm(n, G), VectorFormVol<double>(ag(n,g), areas), Q(isotropic ? Q/(4*M_PI) : Q)
    {};

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
    
    VectorFormVol<double>* clone()
    {
      return new SNWeakForm::VolumetricExternalSourceVF(*this);
    }
  };
  
  class InterfaceStreamingMF : protected GenericForm, public MatrixFormSurf<double>
  {
  public:
    InterfaceStreamingMF(unsigned int n, unsigned int g, unsigned int G) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(n,g), H2D_DG_INNER_EDGE)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormSurf<double>* clone()
    {
      return new SNWeakForm::InterfaceStreamingMF(*this);
    }
  };

  class BoundaryStreamingMF : protected GenericForm, public MatrixFormSurf<double>
  { 
  public:
    BoundaryStreamingMF(unsigned int n, unsigned int g, unsigned int G) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(n,g), HERMES_ANY)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormSurf<double>* clone()
    {
      return new SNWeakForm::BoundaryStreamingMF(*this);
    }
  };
  
  class BoundaryStreamingVF : protected GenericForm, public VectorFormSurf<double>
  {    
  public:
    BoundaryStreamingVF(unsigned int n, unsigned int g, unsigned int G, 
                        const Hermes::vector<std::string>& inflow_boundaries) 
      : GenericForm(n, G), VectorFormSurf<double>(ag.pos(n,g), inflow_boundaries)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone()
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
                         const Hermes::vector<std::string>& reflective_boundaries,
                         const Hermes::vector<MeshFunction<double>*>& iterates) 
      : GenericForm(n, G), VectorFormSurf<double>(ag.pos(n,g), reflective_boundaries, iterates), 
        odata(odata), g(g)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone()
    {
      return new SNWeakForm::SpecularReflectionVF(*this);
    }
  };
    
  double calculate_a_dot_v(int n, double x, double y, double vx, double vy) const;
  
  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;
  Ord upwind_flux(Ord u_cent, Ord u_neib) const;

  const Mesh* mesh;
  
  unsigned int N, M, G;
  SupportClasses::OrdinatesData odata;
};
