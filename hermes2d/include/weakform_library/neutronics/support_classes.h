#ifndef ___H2D_NEUTRONICS_SUPPORT_CLASSES_H
#define ___H2D_NEUTRONICS_SUPPORT_CLASSES_H

#include "material_properties.h"

#include "views/mesh_view.h"
#include "views/scalar_view.h"
#include "views/order_view.h"

#include "function/filter.h"
#include "integrals/h1.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace SupportClasses
{   
  namespace Common
  {    
    class SourceFilter : public SimpleFilter<double>
    {
      public: 
        
        /* Lazy constructors: the vector of solutions to be filtered will be added by 'assign_solutions'. */
        
        SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : SimpleFilter<double>(), matprop(matprop),
            source_regions(source_regions.begin(), source_regions.end())
        {
          pre_init();
        };
        SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::string& source_region)
          : SimpleFilter<double>(), matprop(matprop)
        { 
          source_regions.insert(source_region);
          pre_init();
        }
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions, 
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop),
            source_regions(source_regions.begin(), source_regions.end())
        {
          post_init();
        };
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop),
            source_regions(source_regions.begin(), source_regions.end())
        {
          post_init();
        };
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::string& source_region)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop)
        { 
          source_regions.insert(source_region); 
          post_init();
        }
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::string& source_region)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop)
        { 
          source_regions.insert(source_region); 
          post_init();
        }
        
        virtual void assign_solutions(const Hermes::vector<Solution<double>*>& solutions);
        virtual void assign_solutions(const Hermes::vector<MeshFunction<double>*>& solutions);
        
        virtual void set_active_element(Element* e);
        
        double integrate(GeomType geom_type = HERMES_PLANAR);
                    
      protected:
        const MaterialProperties::Common::MaterialPropertyMaps& matprop;
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
  }
 
  namespace Diffusion
  {
    class Visualization : public Common::Visualization
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
  }
  
  namespace SPN
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
    
    class SourceFilter : public Common::SourceFilter
    {
      public: 
        
        /* Lazy constructors: the vector of solutions to be filtered will be added by 'assign_solutions'. */
        
        SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : Common::SourceFilter(matprop, source_regions), G(matprop.get_G()), mg(G)
        {};
        SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop, const std::string& source_region)
          : Common::SourceFilter(matprop, source_region), G(matprop.get_G()), mg(G) 
        {};
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions, 
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : Common::SourceFilter(solutions, matprop, source_regions), G(matprop.get_G()), mg(G)
        {};
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::vector<std::string>& source_regions = std::vector<std::string>())
          : Common::SourceFilter(solutions, matprop, source_regions), G(matprop.get_G()), mg(G) 
        {};
        SourceFilter(Hermes::vector<MeshFunction<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::string& source_region)
          : Common::SourceFilter(solutions, matprop, source_region), G(matprop.get_G()), mg(G) 
        {};
        SourceFilter(Hermes::vector<Solution<double>*> solutions,
                    const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                    const std::string& source_region)
          : Common::SourceFilter(solutions, matprop, source_region), G(matprop.get_G()), mg(G) 
        {};
                    
        virtual void assign_solutions(const Hermes::vector<Solution<double>*>& solutions) {
          num = solutions.size();
          Common::SourceFilter::assign_solutions(solutions);
        }
        virtual void assign_solutions(const Hermes::vector<MeshFunction<double>*>& solutions) {
          num = solutions.size();
          Common::SourceFilter::assign_solutions(solutions);
        }
        
      protected:
        unsigned int G;
        MomentGroupFlattener mg;
        
        virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
    };
  
    class Visualization : public Common::Visualization
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
  }
  
  class PostProcessor
  {
    NeutronicsMethod method;
    GeomType geom_type;
    
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunction<double>* solution,
                                                        const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                                        const Hermes::vector<std::string>& regions,
                                                        unsigned int this_group, int other_group = -1) const;
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunction<double>* solution,
                                                        const MaterialProperties::Common::MaterialPropertyMaps& matprop,
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
                                            const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                            const Hermes::vector<std::string>& src_areas = Hermes::vector<std::string>()) const;
                                            
      void normalize_to_unit_power(Hermes::vector<Solution<double>*>* solutions,
                                  const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                  double power_per_fission,
                                  const Hermes::vector<std::string>& src_areas = Hermes::vector<std::string>()) const;
      
                                  
                                  
      void get_integrated_group_reaction_rates( ReactionType reaction, 
                                                const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results,
                                                const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                                unsigned int group, const Hermes::vector<std::string>& regions) const;
                                                
      void get_integrated_group_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results, 
                                              unsigned int group, unsigned int G, 
                                              const Hermes::vector<std::string>& regions) const;
                                              
      void get_integrated_reaction_rates( ReactionType reaction, 
                                          const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results,
                                          const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                          const Hermes::vector<std::string>& regions) const;
                                          
      void get_integrated_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, Hermes::vector<double>* results, 
                                        unsigned int G, const Hermes::vector<std::string>& regions) const;                                                                             
      
      void get_areas(Mesh *mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const;
                                        
                                        
      double get_integrated_group_reaction_rates( ReactionType reaction, const Hermes::vector<Solution<double>*>& solutions,
                                                  const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                                  unsigned int group,
                                                  const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                  
      double get_integrated_group_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions, unsigned int group, unsigned int G,
                                                const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                
      double get_integrated_reaction_rates( ReactionType reaction, const Hermes::vector<Solution<double>*>& solutions,
                                            const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                                            const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                            
      double get_integrated_scalar_fluxes(const Hermes::vector<Solution<double>*>& solutions,
                                          unsigned int G, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                          
      double get_area(Mesh *mesh, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
  };
  
  class SourceIteration
  {
    NeutronicsMethod method;
    const std::vector<std::string>& fission_regions;
    GeomType geom_type;
            
    Common::SourceFilter *new_source, *old_source;
    
    public:
      /// \param[in] fission_regions  Strings specifiying the parts of the solution domain where fission occurs.
      /// \param[in]     hermes2d     Class encapsulating global Hermes2D functions.
      /// \param[in]     spaces       Pointers to spaces on which the solutions are defined (one space for each energy group).
      /// \param[in]     wf           Pointer to the weak form of the problem.
      SourceIteration(NeutronicsMethod method, const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                      const std::vector<std::string>& fission_regions = std::vector<std::string>(),
                      GeomType geom_type = HERMES_PLANAR);
                      
      ~SourceIteration() { delete new_source; delete old_source; }
                      
      // \brief Power iteration method for finding the dominant eigenvalue. 
      ///
      /// Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
      /// has converged, also updating the global eigenvalue 'k_eff'.
      ///
      /// \param[in,out] solution     A set of Solution* pointers to solution components (neutron fluxes in each group). 
      ///                             Initial guess for the iteration on input, converged result on output.
      /// \param[in]     tol          Relative difference between two successive eigenvalue approximations that stops the iteration.
      /// \param[in]    matrix_solver Solver for the resulting matrix problem.
      ///
      /// \return  number of iterations needed for convergence within the specified tolerance.
      ///
      int eigenvalue_iteration(const Hermes::vector<Solution<double> *>& solutions, DiscreteProblem<double>& dp,
                              double tol_keff = 1e-6, double tol_flux = 0, MatrixSolverType matrix_solver = SOLVER_UMFPACK);
  };
  
/* SupportClasses */
}
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}      

#endif