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
    class DefaultWeakFormFixedSource : public WeakForm
    {        
      public:
        DefaultWeakFormFixedSource( Hermes::vector<std::string> regions, 
                                    Hermes::vector<double> D_map, 
                                    Hermes::vector<double> Sigma_a_map, 
                                    Hermes::vector<double> Q_map );
    };
  }
  
  namespace Common
  {
    class WeakFormSourceIteration
    {
      protected:
        double keff;
        
        // TODO: source_areas and keff_iteration_forms initialization also belongs here. 
        WeakFormSourceIteration(double initial_keff_guess) : keff(initial_keff_guess) {};
        
      public:
        virtual void update_keff(double new_keff) = 0; //TODO: Define a common FissionYield::OuterIteration class,
                                                      // so that this method may be defined here instead of in both
                                                      // SPN and Diffusion DefaultWeakFormSourceIteration.
        double get_keff() const { return keff; }
    };
  }
        
  namespace Diffusion
  {      
    using namespace MaterialProperties::Diffusion;
    using namespace WeakFormParts::Diffusion;
          
    class DefaultWeakFormFixedSource : public WeakForm
    {
      protected:
        void lhs_init(unsigned int G, const MaterialPropertyMaps& matprop, GeomType geom_type);
        
      public:
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                  HermesFunction *minus_f_src,
                                  const std::string& src_area = HERMES_ANY,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                  HermesFunction *minus_f_src,
                                  const Hermes::vector<std::string>& src_areas,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                  const std::vector<HermesFunction*>& minus_f_src,
                                  const std::string& src_area = HERMES_ANY,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                  const std::vector<HermesFunction*>& minus_f_src,
                                  const Hermes::vector<std::string>& src_areas,
                                  GeomType geom_type = HERMES_PLANAR);
    };
            
    class DefaultWeakFormSourceIteration : public WeakForm, public Common::WeakFormSourceIteration
    {
      protected:            
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        
        void init(const MaterialPropertyMaps& matprop, const Hermes::vector<MeshFunction*>& iterates, 
                  GeomType geom_type, const Hermes::vector<std::string>& areas = Hermes::vector<std::string>());
        
      public:
        DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                        const Hermes::vector<MeshFunction*>& iterates,
                                        double initial_keff_guess,
                                        GeomType geom_type = HERMES_PLANAR );
                                        
        DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                        const Hermes::vector<Solution*>& iterates,
                                        double initial_keff_guess,
                                        GeomType geom_type = HERMES_PLANAR );    
                                        
        DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                        const Hermes::vector<MeshFunction*>& iterates, 
                                        const Hermes::vector<std::string>& src_areas,
                                        double initial_keff_guess,
                                        GeomType geom_type = HERMES_PLANAR );
                                        
        DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                        const Hermes::vector<Solution*>& iterates, 
                                        const Hermes::vector<std::string>& src_areas,
                                        double initial_keff_guess,
                                        GeomType geom_type = HERMES_PLANAR );                                            
                                        
        void update_keff(double new_keff);
    };
    
  }

  namespace SPN
  {
    using namespace MaterialProperties::SPN;
    using namespace WeakFormParts::SPN;
    
    class WeakFormHomogeneous : public WeakForm
    {
      protected:
        MomentGroupFlattener mg;
        unsigned int G, N_odd;
        
        WeakFormHomogeneous(unsigned int N, const MaterialPropertyMaps& matprop,
                            GeomType geom_type, bool include_fission);
    };
            
    class DefaultWeakFormFixedSource : public WeakFormHomogeneous
    {
      protected:
        std::vector<VectorFormVol*> source_terms;
        
      public:
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                  HermesFunction *minus_isotropic_source,
                                  std::string src_area = HERMES_ANY,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                  HermesFunction *minus_isotropic_source,
                                  Hermes::vector<std::string> src_areas,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                  const std::vector<HermesFunction*>& minus_isotropic_sources,
                                  std::string src_area = HERMES_ANY,
                                  GeomType geom_type = HERMES_PLANAR);
        
        DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                  const std::vector<HermesFunction*>& minus_isotropic_sources,
                                  Hermes::vector<std::string> src_areas,
                                  GeomType geom_type = HERMES_PLANAR);
                                  
        virtual ~DefaultWeakFormFixedSource();
    };
    
    class DefaultWeakFormSourceIteration : public WeakFormHomogeneous, public Common::WeakFormSourceIteration
    {
      protected:
        std::vector<FissionYield::OuterIterationForm*> keff_iteration_forms;
        Hermes::vector<MeshFunction*> scalar_flux_iterates;
        
      public:
        DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop, unsigned int N,
                                        const Hermes::vector<Solution*>& iterates, 
                                        const Hermes::vector<std::string>& src_areas,
                                        double initial_keff_guess,
                                        GeomType geom_type = HERMES_PLANAR );
        
        virtual ~DefaultWeakFormSourceIteration();
        
        void update_keff(double new_keff);
        
        const Hermes::vector<MeshFunction*>& get_scalar_flux_iterates() const { 
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