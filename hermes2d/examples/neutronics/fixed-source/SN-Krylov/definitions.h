#include "hermes2d.h"
#include "../../utils.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Mixins;
using namespace Neutronics;
using namespace Neutronics::SN;







class DiffusionWeakForm : public Common::WeakForms::NeutronicsProblem
{
  public:
    DiffusionWeakForm(const Mesh* mesh, bool DG, const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop,
                      const Hermes::vector<std::string>& void_boundaries = Hermes::vector<std::string>());
    
    virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }                      
    
  private:
    template<typename Real>
    static Real dot2(Real x1, Real y1, Real x2, Real y2) {
      return x1*x2 + y1*y2;
    }
    
    class HomogeneousPart : public Common::WeakForms::HomogeneousPart
    {
      public:
        HomogeneousPart(const Common::MaterialProperties::MaterialPropertyMaps* matprop,
                        GeomType geom_type, bool include_fission);
                  
        friend class DiffusionWeakForm;
    };
    
    class BoundaryJacobian : public MatrixFormSurf<double>
    {
      const Mesh* mesh;
      const Diffusion::MaterialProperties::MaterialPropertyMaps* matprop;
      int theta, C_W;
      
      public:
        BoundaryJacobian(const Mesh* mesh, const Diffusion::MaterialProperties::MaterialPropertyMaps* matprop,
                         int theta, int C_W) : MatrixFormSurf<double>(0,0),
          mesh(mesh), matprop(matprop), theta(theta), C_W(C_W)
        {};
           
        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                            Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;
        
        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const
        {
          return Ord(10);
        } 

        virtual MatrixFormSurf<double>* clone() {
          return new BoundaryJacobian(*this);
        }
    };
    
    class BoundaryResidual : public VectorFormSurf<double>
    {
      double D;
      int theta, C_W;
      
      public:
        BoundaryResidual(double D, int theta, int C_W) 
          : VectorFormSurf<double>(0), D(D), theta(theta), C_W(C_W)
        {};

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                            Geom<double> *e, ExtData<double> *ext) const;                   

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const
        {
          return Ord(10);
        }

        virtual VectorFormSurf<double>* clone() {
          return new BoundaryResidual(*this);
        }
    };
    
    class InterfaceJacobian : public MatrixFormSurf<double>
    {
      const Mesh* mesh;
      const Diffusion::MaterialProperties::MaterialPropertyMaps* matprop;
      int theta, C_W;
      
      public:
        InterfaceJacobian(const Mesh* mesh, const Diffusion::MaterialProperties::MaterialPropertyMaps* matprop, 
                          int theta, int C_W) 
          : MatrixFormSurf<double>(0,0,H2D_DG_INNER_EDGE), mesh(mesh), matprop(matprop), theta(theta), C_W(C_W) 
        {};
        
        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                            Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;
        
        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const
        {
          return Ord(10);
        } 

        virtual MatrixFormSurf<double>* clone() {
          return new InterfaceJacobian(*this);
        }
    };
};















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

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

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
    Real vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return vector_form<double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
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
    Real vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return vector_form<double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
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
    Real vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      return vector_form<double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
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

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
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

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    MatrixFormSurf<double>* clone()
    {
      return new SNWeakForm::BoundaryStreamingMF(*this);
    }
  };

  class SpecularReflectionMF_X : protected GenericForm, public MatrixFormSurf<double>
  { 
    const SupportClasses::OrdinatesData& odata;
  public:
    SpecularReflectionMF_X(const SupportClasses::OrdinatesData& odata, 
                           unsigned int n, unsigned int g, unsigned int G) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(odata.reflections_about_x[n],g), HERMES_ANY), odata(odata)
    {
      //this->scaling_factor = -1;
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    MatrixFormSurf<double>* clone()
    {
      return new SNWeakForm::SpecularReflectionMF_X(*this);
    }
  };
  
  class SpecularReflectionMF_Y : protected GenericForm, public MatrixFormSurf<double>
  { 
    const SupportClasses::OrdinatesData& odata;
  public:
    SpecularReflectionMF_Y(const SupportClasses::OrdinatesData& odata,
                           unsigned int n, unsigned int g, unsigned int G) 
      : GenericForm(n, G), MatrixFormSurf<double>(ag.pos(n,g), ag.pos(odata.reflections_about_y[n],g), HERMES_ANY), odata(odata)
    {
      //this->scaling_factor = -1;
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    MatrixFormSurf<double>* clone()
    {
      return new SNWeakForm::SpecularReflectionMF_Y(*this);
    }
  };
  
  class BoundaryStreamingVF : protected GenericForm, public VectorFormSurf<double>
  {    
  public:
    BoundaryStreamingVF(unsigned int n, unsigned int g, unsigned int G, 
                        const Hermes::vector<std::string>& inflow_boundaries) 
      : GenericForm(n, G), VectorFormSurf<double>(ag.pos(n,g), inflow_boundaries)
    {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

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

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    VectorFormSurf<double>* clone()
    {
      return new SNWeakForm::SpecularReflectionVF(*this);
    }
  };
    
  double calculate_a_dot_v(int n, double vx, double vy) const;
  
  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;
  Ord upwind_flux(Ord u_cent, Ord u_neib) const;

  const Mesh* mesh;
  
  unsigned int N, M, G;
  SupportClasses::OrdinatesData odata;
};
