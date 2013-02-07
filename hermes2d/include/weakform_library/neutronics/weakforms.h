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
        
        NeutronicsProblem(unsigned int n_eq,
                          const MaterialPropertyMaps* matprop,
                          GeomType geom_type)
          : WeakForm<double>(n_eq), 
            matprop(matprop), geom_type(geom_type), G(matprop->get_G()) 
        { };
        
      public:
        virtual ~NeutronicsProblem() { };
        
        virtual NeutronicsMethod get_method_type() const = 0;
        GeomType get_geom_type() const { return geom_type; }
    };
    
    class KeffEigenvalueProblem : public NeutronicsProblem
    {
      protected:
        double keff;

        Hermes::vector<Solution<double>*> stored_flux_solutions;
        Hermes::vector<MeshFunction<double>*> scalar_flux_iterates;
        
        /// \param[in] fission_regions  Strings specifiying the parts of the solution domain where fission occurs.
        KeffEigenvalueProblem(unsigned int n_eq,
                              const MaterialPropertyMaps* matprop,
                              GeomType geom_type, 
                              double initial_keff_guess) 
          : NeutronicsProblem(n_eq, matprop, geom_type),
            keff(initial_keff_guess)
        { };
        
      public:        
        virtual void update_keff(double new_keff) = 0; //TODO: Define a common FissionYield::OuterIteration class,
                                                      // so that this method may be defined here instead of in both
                                                      // SPN and Diffusion KeffEigenvalueProblem.
        virtual void update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed) = 0;
        
        virtual SupportClasses::SourceFilter* create_source_filter() = 0;
        virtual SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<Solution<double>*>& solutions) = 0;
        virtual SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions) = 0;
                                      
        double get_keff() const { return keff; } 
        const Hermes::vector<MeshFunction<double>*>& get_scalar_flux_iterates() const { return scalar_flux_iterates; }
    };
    
  /* WeakForms */  
  }
  /* Common */
  }
        
  namespace Diffusion { namespace WeakForms 
  {      
    using namespace WeakFormParts;
    using namespace MaterialProperties;
    
    class DiffusionWeakForm
    {
      protected:
        enum FissionTreatment { IMPLICIT, EXPLICIT, NONE };
        
        DiffusionWeakForm() {};
        
        void add_forms_nonlinear(Common::WeakForms::NeutronicsProblem* wf, 
                                const MaterialPropertyMaps *mp, GeomType geom_type, FissionTreatment include_fission);
        void add_forms(Common::WeakForms::NeutronicsProblem* wf, 
                      const MaterialPropertyMaps *mp, GeomType geom_type, FissionTreatment include_fission);
    };
          
    class FixedSourceProblem : public Common::WeakForms::NeutronicsProblem, protected DiffusionWeakForm
    {  
      protected:
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton, FissionTreatment include_fission = IMPLICIT) { 
          if (solve_by_newton) 
            add_forms_nonlinear(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
          else
            add_forms(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
        }
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }
        
        virtual WeakForm<double>* clone() const { return new FixedSourceProblem(*this); }
    };
            
    class KeffEigenvalueProblem : public Common::WeakForms::KeffEigenvalueProblem, protected DiffusionWeakForm
    {
      protected:    
        void init_rhs(const Hermes::vector<Solution<double>*>& iterates);
        
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton, FissionTreatment include_fission = NONE) { 
          if (solve_by_newton) 
            add_forms_nonlinear(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
          else
            add_forms(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
        }
        
      public:                                        
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<Solution<double>*>& iterates,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR,
                              bool solve_by_newton = false);                                
                                        
        void update_keff(double new_keff);
        void update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed);
        
        Common::SupportClasses::SourceFilter* create_source_filter() {
          return new SupportClasses::SourceFilter(*matprop, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<Solution<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, geom_type);
        }
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }
        
        // FIXME: Should be OK for multithreaded assembling; may not be usable in keff
        // eigenvalue iteration, since scalar_flux_iterates that are being updated during 
        // that iteration won't point to correct OuterIterationForm ext functions after
        // cloning (brand new ext functions will be created by WeakForm::cloneMembers, but
        // they can't be assigned to scalar_flux_iterates because cloneMembers is called 
        // *after* clone() and is private to WeakForm. Could be fixed by updating directly
        // ext functions of all forms of type OuterIterationForm.
        virtual WeakForm<double>* clone() const { return new KeffEigenvalueProblem(*this); }
    };  
  
  /* WeakForms */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakForms 
  {
    using namespace WeakFormParts; 
    using namespace MaterialProperties;
                                         
    class SPNWeakForm
    {
      protected:
        enum FissionTreatment { IMPLICIT, EXPLICIT, NONE };
        
        unsigned int N, N_odd;
        MomentGroupFlattener mg;
        
        SPNWeakForm(unsigned int N, unsigned int G) : N(N), N_odd((N+1)/2), mg(G) { };
        
        void add_forms_nonlinear(Common::WeakForms::NeutronicsProblem* wf, 
                                 const MaterialPropertyMaps *mp, GeomType geom_type, FissionTreatment include_fission);
        
        void add_forms(Common::WeakForms::NeutronicsProblem* wf, 
                       const MaterialPropertyMaps *mp, GeomType geom_type, FissionTreatment include_fission);
    };
            
    class FixedSourceProblem : public Common::WeakForms::NeutronicsProblem, protected SPNWeakForm
    {
      protected:
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton, FissionTreatment include_fission = IMPLICIT) { 
          if (solve_by_newton) 
            add_forms_nonlinear(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
          else
            add_forms(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission); 
        }
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *minus_isotropic_source,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *minus_isotropic_source,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR,
                           bool solve_by_newton = false);
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_SPN; }
        
        virtual WeakForm<double>* clone() const { return new FixedSourceProblem(*this); }
    };
    
    class KeffEigenvalueProblem : public Common::WeakForms::KeffEigenvalueProblem, protected SPNWeakForm
    {
      protected:        
        virtual void add_forms_from_homogeneous_part(bool solve_by_newton, FissionTreatment include_fission = NONE) { 
          if (solve_by_newton) 
            add_forms_nonlinear(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission);
          else
            add_forms(this, static_cast<const MaterialPropertyMaps*>(matprop), geom_type, include_fission); 
        }
        
      public:
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                              const Hermes::vector<Solution<double>*>& iterates,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR,
                              bool solve_by_newton = false);
        
        virtual ~KeffEigenvalueProblem();
        
        void update_keff(double new_keff);
        void update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed);
                
        Common::SupportClasses::SourceFilter* create_source_filter() {
          return new SupportClasses::SourceFilter(*matprop, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<Solution<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, geom_type);
        }
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_SPN; }
        
        // FIXME: Should be OK for multithreaded assembling; may not be usable in keff
        // eigenvalue iteration, since scalar_flux_iterates that are being updated during 
        // that iteration won't point to correct OuterIterationForm ext functions after
        // cloning (brand new ext functions will be created by WeakForm::cloneMembers, but
        // they can't be assigned to scalar_flux_iterates because cloneMembers is called 
        // *after* clone() and is private to WeakForm. Could be fixed by updating directly
        // ext functions of all forms of type OuterIterationForm.
        virtual WeakForm<double>* clone() const { return new KeffEigenvalueProblem(*this); }
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
