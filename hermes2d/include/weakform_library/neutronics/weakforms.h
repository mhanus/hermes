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
    };
  }
  
  namespace Common { namespace WeakForms 
  {
    using namespace MaterialProperties;
    
    class HomogeneousPart
    {
      protected:
        std::vector<MatrixFormVol<double>*> matrix_forms;
        std::vector<VectorFormVol<double>*> vector_forms;

      public:                          
        virtual ~HomogeneousPart();
        
        friend class NeutronicsProblem;
    };
    
    class NeutronicsProblem : public WeakForm<double>
    {
      protected:
        const MaterialPropertyMaps* matprop;
        GeomType geom_type;
        unsigned int G;
        
        HomogeneousPart *homogeneous_part;
        
        NeutronicsProblem(unsigned int n_eq,
                          const MaterialPropertyMaps* matprop,
                          GeomType geom_type)
          : WeakForm<double>(n_eq), 
            matprop(matprop), geom_type(geom_type), G(matprop->get_G())
        { };
        
        
      public:
        virtual ~NeutronicsProblem() { delete homogeneous_part; };
        
        void add_forms_from_homogeneous_part();
        
        virtual NeutronicsMethod get_method_type() const = 0;
        GeomType get_geom_type() const { return geom_type; }
    };
    
    class KeffEigenvalueProblem : public NeutronicsProblem
    {
      protected:
        double keff;
        const Hermes::vector<std::string>& fission_materials;
        Hermes::vector<std::string> fission_regions;

        Hermes::vector<Solution<double>*> stored_flux_solutions;
        Hermes::vector<MeshFunction<double>*> scalar_flux_iterates;
        
        // TODO: keff_iteration_forms initialization also belongs here. 
        /// \param[in] fission_regions  Strings specifiying the parts of the solution domain where fission occurs.
        KeffEigenvalueProblem(unsigned int n_eq,
                              const MaterialPropertyMaps* matprop,
                              GeomType geom_type, 
                              double initial_keff_guess,
                              const Hermes::vector<std::string>& fission_materials = Hermes::vector<std::string>()) 
          : NeutronicsProblem(n_eq, matprop, geom_type),
            keff(initial_keff_guess), fission_materials(fission_materials)
        { 
          Hermes::vector<std::string>::const_iterator it = fission_materials.begin();
          Hermes::vector<std::string>::iterator insert_it = fission_regions.begin();
          for ( ; it != fission_materials.end(); ++it)
          {
            Hermes::vector<std::string> regs = matprop->get_regions(*it);
            fission_regions.insert(insert_it, regs.begin(), regs.end());
            insert_it += regs.size();
          }
        };
        
      public:
        virtual ~KeffEigenvalueProblem() { };
        
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
    
    class HomogeneousPart : public Common::WeakForms::HomogeneousPart
    {
      public:
        HomogeneousPart(const Common::MaterialProperties::MaterialPropertyMaps* matprop,
                        GeomType geom_type, bool include_fission);
                  
        friend class FixedSourceProblem;
        friend class KeffEigenvalueProblem;
    };
          
    class FixedSourceProblem : public Common::WeakForms::NeutronicsProblem
    {
      protected:
        std::vector<VectorFormVol<double>*> source_terms;
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           Hermes2DFunction<double> *minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR);
                           
        ~FixedSourceProblem();
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }
    };
            
    class KeffEigenvalueProblem : public Common::WeakForms::KeffEigenvalueProblem
    {
      protected:  
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        
        void init_rhs(const Hermes::vector<Solution<double>*>& iterates);
        
      public:                                        
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<Solution<double>*>& iterates,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );    
                                                                                
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<Solution<double>*>& iterates, 
                              const Hermes::vector<std::string>& fission_materials,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );                                            
                                        
        void update_keff(double new_keff);
        void update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed);
        
        Common::SupportClasses::SourceFilter* create_source_filter() {
          return new SupportClasses::SourceFilter(*matprop, fission_regions, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<Solution<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, fission_regions, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, fission_regions, geom_type);
        }
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_DIFFUSION; }
    };  
  
  /* WeakForms */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakForms 
  {
    using namespace WeakFormParts; 
    using namespace MaterialProperties;
    
    class HomogeneousPart : public Common::WeakForms::HomogeneousPart
    {
      public:
        HomogeneousPart(unsigned int N_odd, 
                        const Common::MaterialProperties::MaterialPropertyMaps* matprop,
                        GeomType geom_type, bool include_fission);
                  
        friend class FixedSourceProblem;
        friend class KeffEigenvalueProblem;
    };
    
    class SPNWeakForm
    {
      protected:
        unsigned int N, N_odd;
        MomentGroupFlattener mg;
        
        SPNWeakForm(unsigned int N, unsigned int G) : N(N), N_odd((N+1)/2), mg(G) { };
    };
            
    class FixedSourceProblem : public Common::WeakForms::NeutronicsProblem, protected SPNWeakForm
    {
      protected:
        std::vector<VectorFormVol<double>*> source_terms;
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *minus_isotropic_source,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           Hermes2DFunction<double> *minus_isotropic_source,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR);
                                  
        virtual ~FixedSourceProblem();
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_SPN; }
    };
    
    class KeffEigenvalueProblem : public Common::WeakForms::KeffEigenvalueProblem, protected SPNWeakForm
    {
      protected:
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        
      public:
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                              const Hermes::vector<Solution<double>*>& iterates, 
                              const Hermes::vector<std::string>& fission_regions,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );
        
        virtual ~KeffEigenvalueProblem();
        
        void update_keff(double new_keff);
        void update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed);
                
        Common::SupportClasses::SourceFilter* create_source_filter() {
          return new SupportClasses::SourceFilter(*matprop, fission_regions, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<Solution<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, fission_regions, geom_type);
        }
        Common::SupportClasses::SourceFilter* create_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions) {
          return new SupportClasses::SourceFilter(solutions, *matprop, fission_regions, geom_type);
        }
        
        virtual NeutronicsMethod get_method_type() const { return NEUTRONICS_SPN; }
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
