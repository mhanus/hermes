#ifndef ___H2D_NEUTRONICS_WEAK_FORMS_H
#define ___H2D_NEUTRONICS_WEAK_FORMS_H

#include "weakform_parts_implementation.h"
#include "weakforms_h1.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics 
{
  namespace SimpleMonoenergeticDiffusionWeakForms
  {
    /* 
    Simple monoenergetic neutron diffusion, with the following weak formulation within each
    homogeneous region:
    
    \int_{region} D \nabla\phi \cdot \nabla\psi d\bfx + \int_{region} \Sigma_a \phi\psi d\bfx
    + \int_{region} Q_{ext}\psi d\bfx = 0
    
    where 
    
    D         ... diffusion coefficient, 
    \Sigma_a  ... absorption cross-section, 
    Q_{ext}   ... external neutron sources (multiplied with -1)
    
    are region-wise constant physical parameters of the problem. Each region has one entry in vector
    'regions', which is the marker used for all elements it is composed of (usually specified in the
    mesh file). A corresponding entry in the *_map arguments is the value of the particular physical 
    parameter for that marker.
    
    Dirichlet and/or zero Neumann BC are assumed - nonzero Neumann or Newton boundary conditions can 
    be enabled by creating a descendant and adding surface forms to it.
    */
    class FixedSourceProblem : public WeakForm<double>
    {        
      public:
        FixedSourceProblem(Hermes::vector<std::string> regions, 
                           Hermes::vector<double> D_map, 
                           Hermes::vector<double> Sigma_a_map, 
                           Hermes::vector<double> Q_map );
                           
        virtual WeakForm<double>* clone() const { return new FixedSourceProblem(*this); }
    };
  }
  
  namespace Common { namespace WeakForms 
  {
    using namespace MaterialProperties;
    
    class NeutronicsProblem : public WeakForm<double>
    {
      protected:
        const MaterialPropertyMaps* matprop;
        GeomType geom_type;
        unsigned int G;
        
        enum FissionTreatment { NONE=0, IMPLICIT, EXPLICIT };
        FissionTreatment include_fission;
        
        NeutronicsProblem(unsigned int n_eq,
                          const MaterialPropertyMaps* matprop,
                          GeomType geom_type)
          : WeakForm<double>(n_eq), 
            matprop(matprop), geom_type(geom_type), G(matprop->get_G()), include_fission(IMPLICIT) 
        { };
        
      public:
        virtual ~NeutronicsProblem() { };
        
        virtual NeutronicsMethod get_method_type() const = 0;
        GeomType get_geom_type() const { return geom_type; }
        
        void get_source_part(WeakForm<double>* wf);
    };    
  /* WeakForms */  
  }
  /* Common */
  }
        
  namespace Diffusion { namespace WeakForms 
  {      
    using namespace WeakFormParts;
    using namespace MaterialProperties;
    
    class DiffusionWeakForm : public Common::WeakForms::NeutronicsProblem
    {
      protected:
        DiffusionWeakForm(const MaterialPropertyMaps *mp, GeomType geom_type) 
          : Common::WeakForms::NeutronicsProblem(mp->get_G(), mp, geom_type)
        {};
        
        void add_forms_nonlinear();
        void add_forms();
                      
      public:
        void add_fission_sparse_structure();
        
        void get_fission_yield_part(WeakForm<double>* wf);
        void get_diffusion_reaction_part(WeakForm<double>* wf, 
                                         const Hermes::vector<std::string>& vacuum_boundaries = Hermes::vector<std::string>());
        void get_scattering_part(WeakForm<double>* wf);
    };
          
    class FixedSourceProblem : public DiffusionWeakForm
    {  
      protected:
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton) { 
          if (solve_by_newton) 
            add_forms_nonlinear();
          else
            add_forms();
        }
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }
        
        virtual WeakForm<double>* clone() const { return new FixedSourceProblem(*this); }
    };
  
  /* WeakForms */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakForms 
  {
    using namespace WeakFormParts; 
    using namespace MaterialProperties;
                                         
    class SPNWeakForm : public Common::WeakForms::NeutronicsProblem
    {
      protected:        
        unsigned int N, N_odd;
        MomentGroupFlattener mg;
        bool2 sym;
        bool2 present;
        
        SPNWeakForm(unsigned int N, const MaterialPropertyMaps *mp, GeomType geom_type) 
          : Common::WeakForms::NeutronicsProblem(mp->get_G()*(N+1)/2, mp, geom_type), 
            N(N), N_odd((N+1)/2), mg(mp->get_G()) 
        {
          determine_symmetries();
        };
        
        void determine_symmetries();
        void add_forms_nonlinear();
        void add_forms();
                       
      public:
        void add_fission_sparse_structure();
        
        void get_fission_yield_part(WeakForm<double>* wf);
        void get_symmetric_diffusion_reaction_part(WeakForm<double>* wf, 
                                                   const Hermes::vector<std::string>& vacuum_boundaries = Hermes::vector<std::string>());
        void get_nonsymmetric_diffusion_reaction_part(WeakForm<double>* wf);
    };
            
    class FixedSourceProblem : public SPNWeakForm
    {
      protected:
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton) { 
          if (solve_by_newton) 
            add_forms_nonlinear();
          else
            add_forms(); 
        }
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *iso_src,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *iso_src,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& iso_src,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& iso_src,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_SPN; }
        
        virtual WeakForm<double>* clone() const { return new FixedSourceProblem(*this); }
    };
  /* WeakForms */
  }
  /* SPN */ 
  }
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}      

#endif
