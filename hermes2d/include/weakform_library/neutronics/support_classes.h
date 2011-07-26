#ifndef ___H2D_NEUTRONICS_SUPPORT_CLASSES_H
#define ___H2D_NEUTRONICS_SUPPORT_CLASSES_H

#include "material_properties.h"

#include "views/mesh_view.h"
#include "views/scalar_view.h"
#include "views/order_view.h"

#include "function/filter.h"
#include "integrals/h1.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{   
  namespace Common { namespace SupportClasses
  {    
    class SourceFilter : public SimpleFilter<double>
    {
      public: 
        
        /* Lazy constructors: the vector of solutions to be filtered will be added by 'assign_solutions'. */
        
        SourceFilter(const MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(), matprop(matprop), geom_type(geom_type),
            source_regions(source_regions.begin(), source_regions.end())
        {
          pre_init();
        };
        SourceFilter(const MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(), matprop(matprop), geom_type(geom_type)
        { 
          source_regions.insert(source_region);
          pre_init();
        }
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions, 
                     const MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop), geom_type(geom_type),
            source_regions(source_regions.begin(), source_regions.end())
        {
          post_init();
        };
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                     const MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop), geom_type(geom_type),
            source_regions(source_regions.begin(), source_regions.end())
        {
          post_init();
        };
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions,
                    const MaterialProperties::MaterialPropertyMaps& matprop,
                    const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop), geom_type(geom_type)
        { 
          source_regions.insert(source_region); 
          post_init();
        }
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                    const MaterialProperties::MaterialPropertyMaps& matprop,
                    const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop), geom_type(geom_type)
        { 
          source_regions.insert(source_region); 
          post_init();
        }
        
        /// \brief Empty virtual destructor.
        /// Required in order to properly delete derived classes accessed through a pointer to this class.
        virtual ~SourceFilter() {}
        
        virtual void assign_solutions(const Hermes::vector<Solution<double>*>& solutions);
        virtual void assign_solutions(const Hermes::vector<MeshFunction<double>*>& solutions);
        
        virtual void set_active_element(Element* e);
        
        double integrate();
                    
      protected:
        const MaterialProperties::MaterialPropertyMaps& matprop;
        GeomType geom_type;
        std::set<std::string> source_regions;
        std::set<int> markers;
        bool have_solutions;
        
        virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
        
        virtual void pre_init();
        virtual void post_init();
    };
  
    class Visualization
    {
      protected:
        unsigned int n_unknowns, n_equations, n_groups;
        bool display_meshes;
        
        Views::ScalarView<double>** sviews;
        Views::OrderView<double>** oviews;
        Views::MeshView** mviews;
        
        static const std::string  base_title_flux;
        static const std::string  base_title_order;
        static const std::string  base_title_mesh;
        
        std::string itos(int t)
        {
          std::stringstream ss; ss << t;
          return ss.str();
        }
        
        void init(unsigned int nu, unsigned int ne, unsigned int ng);
        
      public:
        Visualization(bool display_meshes = false) : display_meshes(display_meshes) { init(0,0,0); }
        Visualization(unsigned int n_unknowns, unsigned int n_equations, unsigned int n_groups, 
                      bool display_meshes = false) : display_meshes(display_meshes) { init(n_unknowns, n_equations, n_groups); }
        
        virtual ~Visualization();
        
        virtual void show_meshes(Hermes::vector<Mesh*> meshes) = 0;
        virtual void show_solutions(Hermes::vector< Solution<double>* > solutions) = 0;
        virtual void show_orders(Hermes::vector<Space<double>*> spaces) = 0;
        
        void inspect_meshes(Hermes::vector<Mesh*> meshes);
        void inspect_solutions(Hermes::vector< Solution<double>* > solutions);
        void inspect_orders(Hermes::vector<Space<double>*> spaces);
    };
  
  /* SupportClasses */
  }
  /* Common */
  }
 
  namespace Diffusion { namespace SupportClasses
  {
    using Common::SupportClasses::SourceFilter;
    
    class Visualization : public Common::SupportClasses::Visualization
    {
      public:
        Visualization(unsigned int G, bool display_meshes = false);
        
        void show_meshes(Hermes::vector<Mesh*> meshes);
        void show_solutions(Hermes::vector< Solution<double>* > solutions);
        void show_orders(Hermes::vector<Space<double>*> spaces);     
        
        void save_solutions_vtk(const std::string& base_filename, const std::string& base_varname,
                                Hermes::vector< Solution<double>* > solutions,  bool mode_3D = false);
        void save_orders_vtk(const std::string& base_filename, Hermes::vector<Space<double>*> spaces);       
    };
    
  /* SupportClasses */
  }
  /* Diffusion */
  }
  
  namespace SPN { namespace SupportClasses
  {    
    class Coeffs
    {
      private:
        
        static const unsigned int N_MAX = 5;
        
        static const double SYSTEM_MATRIX[N_MAX][N_MAX][N_MAX];
        static const double D_GRAD_F[N_MAX][N_MAX];
        static const double EVEN_MOMENTS[N_MAX][N_MAX];
                    
      public:
        
        static double system_matrix(unsigned int m, unsigned int n, unsigned int k_of_Sigma_t2k);
        static double D_grad_F(unsigned int m, unsigned int n);
        static double even_moment(unsigned int m, unsigned int n);
        
        static double D(unsigned int m) { 
          return 1./(4*m+3); 
        }
        
        static int max_order() { 
          return N_MAX; 
        }
    };
    
    class MomentGroupFlattener
    {
      unsigned int G;
      
      public:
        MomentGroupFlattener() : G(0) {};
        MomentGroupFlattener(unsigned int G) : G(G) {};
        
        void set_G(unsigned int G) { this->G = G; }
        unsigned int get_G() const { return G; }
        
        unsigned int pos(unsigned int angular_moment, unsigned int group) const {
          return angular_moment * G + group;
        }
    };
    
    struct MomentFilter
    {
      class Common 
      {
        protected:
          Common(unsigned int angular_moment, unsigned int group, unsigned int G);  
          
          unsigned int odd_req_mom, req_mom_idx, g;
          MomentGroupFlattener mg;
      };
      
      class Val : protected Common, public SimpleFilter<double>
      {
        public:       
          Val(unsigned int angular_moment, unsigned int group, unsigned int G, 
              const Hermes::vector<MeshFunction<double>*>& solutions)
            : Common(angular_moment, group, G), SimpleFilter<double>(solutions, Hermes::vector<int>())
          {};
          Val(unsigned int angular_moment, unsigned int group, unsigned int G,
              const Hermes::vector<Solution<double>*>& solutions)
            : Common(angular_moment, group, G), SimpleFilter<double>(solutions, Hermes::vector<int>())
          {};
          
          virtual void set_active_element(Element* e);
          
        protected:             
          void filter_fn(int n, Hermes::vector<double*> values, double* result);
      };
      
      class ValDxDy : protected Common, public DXDYFilter<double>
      {
        public:
          ValDxDy(unsigned int angular_moment, unsigned int group, unsigned int G, 
                  const Hermes::vector<MeshFunction<double>*>& solutions)
            : Common(angular_moment, group, G), DXDYFilter<double>(solutions)
          {};
          ValDxDy(unsigned int angular_moment, unsigned int group, unsigned int G,
                  const Hermes::vector<Solution<double>*>& solutions)
            : Common(angular_moment, group, G), DXDYFilter<double>(solutions)
          {};
          
          virtual void set_active_element(Element* e);
          
        protected:
          void filter_fn(int n, 
                        Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, 
                        double* rslt, double* rslt_dx, double* rslt_dy);
      };
      
      static void get_scalar_fluxes(const Hermes::vector<Solution<double>*>& angular_fluxes,
                                    Hermes::vector<MeshFunction<double>*>* scalar_fluxes,
                                    unsigned int G);               
      static void get_scalar_fluxes_with_derivatives(const Hermes::vector<Solution<double>*>& angular_fluxes,
                                                    Hermes::vector<MeshFunction<double>*>* scalar_fluxes,
                                                    unsigned int G);                                                       
      static void clear_scalar_fluxes(Hermes::vector<MeshFunction<double>*>* scalar_fluxes);
    };
    
    class SourceFilter : public Common::SupportClasses::SourceFilter
    {
      public: 
        
        /* Lazy constructors: the vector of solutions to be filtered will be added by 'assign_solutions'. */
        
        SourceFilter(const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(matprop, source_regions, geom_type),
            G(matprop.get_G()), mg(G)
        {};
        SourceFilter(const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(matprop, source_region, geom_type),
            G(matprop.get_G()), mg(G) 
        {};
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions, 
                     const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(solutions, matprop, source_regions, geom_type), 
            G(matprop.get_G()), mg(G)
        {};
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                     const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::vector<std::string>& source_regions = std::vector<std::string>(),
                     GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(solutions, matprop, source_regions, geom_type), 
            G(matprop.get_G()), mg(G) 
        {};
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions,
                     const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(solutions, matprop, source_region, geom_type),
            G(matprop.get_G()), mg(G) 
        {};
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                     const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     const std::string& source_region, GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(solutions, matprop, source_region, geom_type), 
            G(matprop.get_G()), mg(G) 
        {};
                    
        virtual void assign_solutions(const Hermes::vector<Solution<double>*>& solutions) {
          this->num = solutions.size();
          Common::SupportClasses::SourceFilter::assign_solutions(solutions);
        }
        virtual void assign_solutions(const Hermes::vector<MeshFunction<double>*>& solutions) {
          this->num = solutions.size();
          Common::SupportClasses::SourceFilter::assign_solutions(solutions);
        }
        
      protected:
        unsigned int G;
        MomentGroupFlattener mg;
        
        virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
    };
  
    class Visualization : public Common::SupportClasses::Visualization
    {
      unsigned int n_moments, n_odd_moments;
      MomentGroupFlattener mg;
      
      public:
        Visualization(unsigned int spn_order, unsigned int G, bool display_meshes = false);
        
        void show_meshes(Hermes::vector<Mesh*> meshes);
        void show_solutions(Hermes::vector< Solution<double>* > solutions);
        void show_orders(Hermes::vector<Space<double>*> spaces);
        
        void save_solutions_vtk(const std::string& base_filename, const std::string& base_varname,
                                Hermes::vector< Solution<double>* > solutions, bool mode_3D = false);                   
        void save_orders_vtk(const std::string& base_filename, Hermes::vector<Space<double>*> spaces);
        
        void inspect_solutions(Hermes::vector< Solution<double>* > solutions);
    };       
    
  /* SupportClasses */
  }
  /* SPN */
  }
  
  // FIXME: Ad-hoc class. Replace accordingly as soon asarithmetic operations with solutions will be ready.
  template <typename Scalar>
  class MultipliableSolution : public Solution<Scalar>
  {
    public:
      MultipliableSolution(const Solution<Scalar>* solution) { this->copy(solution); }
      
      void multiply(double coeff)
      {
        if (this->sln_type == HERMES_SLN)
          for (int i = 0; i < this->num_coefs; i++)
            this->mono_coefs[i] *= coeff;
        else
          error("Not implemented.");
      }
  };
  
  class PostProcessor
  {
    NeutronicsMethod method;
    GeomType geom_type;
    
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunction<double>* solution,
                                                        const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                        const Hermes::vector<std::string>& regions,
                                                        unsigned int this_group, int other_group = -1) const;
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunction<double>* solution,
                                                        const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                        const std::string& region,
                                                        unsigned int this_group, int other_group = -1) const 
    {
      Hermes::vector<std::string> regions; 
      regions.push_back(region);
      return get_integrated_group_reaction_rates_internal(reaction, solution, matprop, regions, this_group, other_group);
    }
    
    public:
      PostProcessor(NeutronicsMethod method, GeomType geom_type = HERMES_PLANAR) : method(method), geom_type(geom_type) {};
      
      double integrate(MeshFunction<double>* solution, const Hermes::vector<std::string>& areas = Hermes::vector<std::string>()) const;
      double integrate(MeshFunction<double>* solution, const std::string& area) const 
      {
        Hermes::vector<std::string> areas; 
        areas.push_back(area);
        return integrate(solution, areas);
      }
      
      
      void normalize_to_unit_fission_source(Hermes::vector<Solution<double>*>* solutions, 
                                            double integrated_fission_source) const;
                                            
      void normalize_to_unit_fission_source(Hermes::vector<Solution<double>*>* solutions,
                                            const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                            const Hermes::vector<std::string>& src_areas = Hermes::vector<std::string>()) const;
                                            
      void normalize_to_unit_power(Hermes::vector<Solution<double>*>* solutions,
                                  const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                  double power_per_fission,
                                  const Hermes::vector<std::string>& src_areas = Hermes::vector<std::string>()) const;
      
                                  
                                  
      void get_integrated_group_reaction_rates( ReactionType reaction, 
                                                const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results,
                                                const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                unsigned int group, const Hermes::vector<std::string>& regions) const;
                                                
      void get_integrated_group_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results, 
                                              unsigned int group, unsigned int G, 
                                              const Hermes::vector<std::string>& regions) const;
                                              
      void get_integrated_reaction_rates( ReactionType reaction, 
                                          const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results,
                                          const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                          const Hermes::vector<std::string>& regions) const;
                                          
      void get_integrated_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results, 
                                        unsigned int G, const Hermes::vector<std::string>& regions) const;                                                                             
      
      void get_areas(Mesh *mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const;
                                        
                                        
      double get_integrated_group_reaction_rates( ReactionType reaction, const Hermes::vector<Solution<double>*>& solutions,
                                                  const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                  unsigned int group,
                                                  const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                  
      double get_integrated_group_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, unsigned int group, unsigned int G,
                                                const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                
      double get_integrated_reaction_rates( ReactionType reaction, const Hermes::vector<Solution<double>*>& solutions,
                                            const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                            const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                            
      double get_integrated_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions,
                                          unsigned int G, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                          
      double get_area(Mesh *mesh, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
  };
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}      

#endif