#ifndef ___H2D_NEUTRONICS_WEAK_FORMS_H
#define ___H2D_NEUTRONICS_WEAK_FORMS_H

#include "common_definitions.h"
#include "material_properties.h"
#include "support_classes.h"
#include "weakform_parts.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace WeakForms 
{
  namespace Simple
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
  
  namespace Common
  {
    class KeffEigenvalueProblem
    {
      protected:
        double keff;
        
        // TODO: source_areas and keff_iteration_forms initialization also belongs here. 
        KeffEigenvalueProblem(double initial_keff_guess) : keff(initial_keff_guess) {};
        
      public:
        virtual void update_keff(double new_keff) = 0; //TODO: Define a common FissionYield::OuterIteration class,
                                                      // so that this method may be defined here instead of in both
                                                      // SPN and Diffusion KeffEigenvalueProblem.
        double get_keff() const { return keff; }
    };
  }
        
  namespace Diffusion
  {      
    using namespace MaterialProperties::Diffusion;   // Contained in WeakFormParts::Diffusion; 
                                                     // explicitly needed only for KDevelop's Intellisense.
    using namespace WeakFormParts::Diffusion;    
    using DataStructures::bool1;
    using DataStructures::bool2;
          
    class FixedSourceProblem : public WeakForm<double>
    {
      protected:
        void lhs_init(unsigned int G, const MaterialPropertyMaps& matprop, GeomType geom_type);
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           HermesFunction<double> *minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           HermesFunction<double> *minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<HermesFunction<double>*>& minus_f_src,
                           const std::string& src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                           const std::vector<HermesFunction<double>*>& minus_f_src,
                           const Hermes::vector<std::string>& src_areas,
                           GeomType geom_type = HERMES_PLANAR);
    };
            
    class KeffEigenvalueProblem : public WeakForm<double>, public Common::KeffEigenvalueProblem
    {
      protected:            
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        
        void init(const MaterialPropertyMaps& matprop, const Hermes::vector<MeshFunction<double>*>& iterates, 
                  GeomType geom_type, const Hermes::vector<std::string>& areas = Hermes::vector<std::string>());
        
      public:
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );
                                        
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<Solution<double>*>& iterates,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );    
                                        
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<MeshFunction<double>*>& iterates, 
                              const Hermes::vector<std::string>& src_areas,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );
                                        
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                              const Hermes::vector<Solution<double>*>& iterates, 
                              const Hermes::vector<std::string>& src_areas,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );                                            
                                        
        void update_keff(double new_keff);
    };
    
  }

  namespace SPN
  {
    using namespace MaterialProperties::SPN;  // Contained in WeakFormParts::SPN;  
                                              // explicitly needed only for KDevelop's Intellisense.
    using namespace WeakFormParts::SPN; 
    using DataStructures::bool1;
    using DataStructures::bool2;
    
    class WeakFormHomogeneous : public WeakForm<double>
    {
      protected:
        MomentGroupFlattener mg;
        unsigned int G, N_odd;
        
        WeakFormHomogeneous(unsigned int N, const MaterialPropertyMaps& matprop,
                            GeomType geom_type, bool include_fission);
    };
            
    class FixedSourceProblem : public WeakFormHomogeneous
    {
      protected:
        std::vector<VectorFormVol<double>*> source_terms;
        
      public:
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           HermesFunction<double> *minus_isotropic_source,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           HermesFunction<double> *minus_isotropic_source,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<HermesFunction<double>*>& minus_isotropic_sources,
                           std::string src_area = HERMES_ANY,
                           GeomType geom_type = HERMES_PLANAR);
        
        FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                           const std::vector<HermesFunction<double>*>& minus_isotropic_sources,
                           Hermes::vector<std::string> src_areas,
                           GeomType geom_type = HERMES_PLANAR);
                                  
        virtual ~FixedSourceProblem();
    };
    
    class KeffEigenvalueProblem : public WeakFormHomogeneous, public Common::KeffEigenvalueProblem
    {
      protected:
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        Hermes::vector<MeshFunction<double>*> scalar_flux_iterates;
        
      public:
        KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                              const Hermes::vector<Solution<double>*>& iterates, 
                              const Hermes::vector<std::string>& src_areas,
                              double initial_keff_guess,
                              GeomType geom_type = HERMES_PLANAR );
        
        virtual ~KeffEigenvalueProblem();
        
        void update_keff(double new_keff);
        
        const Hermes::vector<MeshFunction<double>*>& get_scalar_flux_iterates() const { 
          return scalar_flux_iterates; 
        }
    };
  }
        
/* WeakForms */
}
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}      

#endif