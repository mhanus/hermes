#ifndef ___H2D_NEUTRONICS_SUPPORT_CLASSES_H
#define ___H2D_NEUTRONICS_SUPPORT_CLASSES_H

#include "material_properties.h"

#include "views/mesh_view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
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
                     GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(), matprop(matprop), geom_type(geom_type),
            source_regions(matprop.get_fission_regions().begin(), matprop.get_fission_regions().end())
        {
          pre_init();
        };
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, 
                     const MaterialProperties::MaterialPropertyMaps& matprop,
                     GeomType geom_type = HERMES_PLANAR)
          : SimpleFilter<double>(solutions, Hermes::vector<int>()), matprop(matprop), geom_type(geom_type),
            source_regions(matprop.get_fission_regions().begin(), matprop.get_fission_regions().end())
        {
          // We need to setup the array 'item' manually, since 'solutions' may be a vector of
          // freshly created Solution's, which have unset num_components.
          for (int i = 0; i < 10; i++)
            item[i] = H2D_FN_VAL & H2D_FN_COMPONENT_0;
          
          post_init();
        };
        
        /// \brief Empty virtual destructor.
        /// Required in order to properly delete derived classes accessed through a pointer to this class.
        virtual ~SourceFilter() {}
        
        virtual void assign_solutions(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions);
        
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
  
    class Visualization : public Hermes::Mixins::Loggable
    {
      protected:
        unsigned int n_equations, n_groups;
        unsigned int width, height;
        bool display_meshes;
        
        Views::ScalarView** sviews;
        Views::OrderView** oviews;
        Views::MeshView** mviews;
        
        static const std::string  base_title_flux;
        static const std::string  base_title_order;
        static const std::string  base_title_mesh;
        
        std::string itos(int t)
        {
          std::stringstream ss; ss << t;
          return ss.str();
        }
        
        void init(unsigned int ne, unsigned int ng);
        
      public:
        Visualization(unsigned int width = 450, unsigned int height = 450, bool display_meshes = false) 
          : width(width), height(height), display_meshes(display_meshes) { init(0,0); }
        Visualization(unsigned int n_equations, unsigned int n_groups, unsigned int width = 450, unsigned int height = 450, bool display_meshes = false) 
          : width(width), height(height), display_meshes(display_meshes) { init(n_equations, n_groups); }
        
        virtual ~Visualization();
        
        virtual void show_meshes(Hermes::vector<MeshSharedPtr> meshes) = 0;
        virtual void show_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions) = 0;
        virtual void show_orders(Hermes::vector<SpaceSharedPtr<double> > spaces) = 0;

#ifndef NOGLUT
        void inspect_meshes(Hermes::vector<MeshSharedPtr> meshes);
        void inspect_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions);
        void inspect_orders(Hermes::vector<SpaceSharedPtr<double> > spaces);
#else
        void inspect_meshes(Hermes::vector<MeshSharedPtr> meshes) {};
        void inspect_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions) {};
        void inspect_orders(Hermes::vector<SpaceSharedPtr<double> > spaces) {};
#endif
        
        Views::ScalarView** get_solution_views(unsigned int* num) { *num = n_equations; return sviews; }
        Views::OrderView** get_order_views(unsigned int* num)     { *num = n_equations; return oviews; }
        Views::MeshView** get_mesh_views(unsigned int* num)       { *num = n_equations; return mviews; }
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
      protected:
#ifndef NOGLUT
        void init(unsigned int G, unsigned int width, unsigned int height, bool display_meshes);
#else
        void init(unsigned int G, unsigned int width, unsigned int height, bool display_meshes) {};
#endif
      public:
        Visualization(unsigned int G, bool display_meshes = false)
          : Common::SupportClasses::Visualization(G, G, 450, 450, display_meshes)
        {
          init(G, 450, 450, display_meshes); 
        }
        Visualization(unsigned int G, unsigned int width, unsigned int height, bool display_meshes = false)
          : Common::SupportClasses::Visualization(G, G, width, height, display_meshes)
        {
          init(G, width, height, display_meshes); 
        }
        
#ifndef NOGLUT
        void show_meshes(Hermes::vector<MeshSharedPtr> meshes);
        void show_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions);
        void show_orders(Hermes::vector<SpaceSharedPtr<double> > spaces);     
#else
        void show_meshes(Hermes::vector<MeshSharedPtr> meshes) {};
        void show_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions) {};
        void show_orders(Hermes::vector<SpaceSharedPtr<double> > spaces) {};     
#endif
        
        void save_solutions_vtk(const std::string& base_filename, const std::string& base_varname,
                                Hermes::vector< MeshFunctionSharedPtr<double> > solutions,  bool mode_3D = false);
        void save_orders_vtk(const std::string& base_filename, Hermes::vector<SpaceSharedPtr<double> > spaces);       
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
          Common(unsigned int angular_moment, unsigned int group, unsigned int G) 
            : odd_req_mom((angular_moment%2) == 1), req_mom_idx(angular_moment/2),  g(group), mg(G)
          {
            if (group >= G) ErrorHandling::error_function("MomentFilter::Common > %s", Messages::E_INVALID_GROUP_INDEX);
          }
          
          unsigned int odd_req_mom, req_mom_idx, g;
          MomentGroupFlattener mg;
      };
      
      class EvenMomentVal : protected Common, public SimpleFilter<double>
      {
        public:       
          EvenMomentVal(unsigned int angular_moment, unsigned int group, unsigned int G, 
                        const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions)
            : Common(angular_moment, group, G), SimpleFilter<double>(solutions, Hermes::vector<int>()),
              angular_moment(angular_moment), group(group), G(G)
          {
            if (odd_req_mom) ErrorHandling::error_function("MomentFilter::EvenMomentVal constructor > %s", Messages::E_EVEN_MOMENT_EXPECTED);
          };
          
          virtual MeshFunction<double>* clone() const;
          
          virtual void set_active_element(Element* e);
          
        protected:             
          void filter_fn(int n, Hermes::vector<double*> values, double* result);
          
          unsigned int angular_moment;
          unsigned int group;
          unsigned int G;
      };
      
      class EvenMomentValDxDy : protected Common, public DXDYFilter<double>
      {
        public:
          EvenMomentValDxDy(unsigned int angular_moment, unsigned int group, unsigned int G, 
                  const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions)
            : Common(angular_moment, group, G), DXDYFilter<double>(solutions),
              angular_moment(angular_moment), group(group), G(G)
          {
            if (odd_req_mom) ErrorHandling::error_function("MomentFilter::EvenMomentVal constructor > %s", Messages::E_EVEN_MOMENT_EXPECTED);
          };
          
          virtual MeshFunction<double>* clone() const;
          
          virtual void set_active_element(Element* e);
          
        protected:
          void filter_fn(int n, 
                         Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, 
                         double* rslt, double* rslt_dx, double* rslt_dy);
                         
          unsigned int angular_moment;
          unsigned int group;
          unsigned int G;
      };
      
      class OddMomentVal : protected Common, public Filter<double>
      {
        public:
          OddMomentVal(unsigned int component, unsigned int angular_moment, unsigned int group, unsigned int G, 
                  const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
                  const MaterialProperties::MaterialPropertyMaps *matprop)
            : Common(angular_moment, group, G), Filter<double>(solutions), 
              component(component), matprop(matprop), angular_moment(angular_moment), group(group), G(G)
          {
            if (!odd_req_mom) ErrorHandling::error_function("MomentFilter::OddMomentVal constructor > %s", Messages::E_ODD_MOMENT_EXPECTED);
            if (component >= 2) ErrorHandling::error_function("MomentFilter::OddMomentVal > %s", Messages::E_INVALID_COMPONENT);
          };
          
          virtual MeshFunction<double>* clone() const;
          
          virtual void set_active_element(Element* e);
          
        protected:
          virtual void precalculate(int order, int mask);
          virtual Func<double>* get_pt_value(double x, double y, Element* e = NULL)
          { 
            ErrorHandling::error_function("Not implemented yet"); 
            return NULL; 
          }
                         
          const MaterialProperties::MaterialPropertyMaps *matprop;
          unsigned int component;
          unsigned int angular_moment;
          unsigned int group;
          unsigned int G;
      };
      
      static void get_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                    Hermes::vector<MeshFunctionSharedPtr<double> >* scalar_fluxes,
                                    unsigned int G);               
      static void get_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                    Hermes::vector<Filter<double>*>* scalar_fluxes,
                                    unsigned int G);  
      static void get_scalar_fluxes_with_derivatives(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                                    Hermes::vector<MeshFunctionSharedPtr<double> >* scalar_fluxes,
                                                    unsigned int G);
      static void get_scalar_fluxes_with_derivatives(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                                    Hermes::vector<Filter<double>*>* scalar_fluxes,
                                                    unsigned int G);                          
      // DEPRECATED (not needed when shared pointers are used for scalar_fluxes)                                              
      static void clear_scalar_fluxes(Hermes::vector< MeshFunction<double>* >* scalar_fluxes);
      static void clear_scalar_fluxes(Hermes::vector<Filter<double>*>* scalar_fluxes);
    };
    
    class SourceFilter : public Common::SupportClasses::SourceFilter
    {
      public: 
        
        /* Lazy constructors: the vector of solutions to be filtered will be added by 'assign_solutions'. */
        
        SourceFilter(const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(matprop, geom_type),
            G(matprop.get_G()), mg(matprop.get_G())
        {};
        
        /* Immediate constructors: the vector of solutions to be filtered is given by the first argument.  */
        
        SourceFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, 
                     const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                     GeomType geom_type = HERMES_PLANAR)
          : Common::SupportClasses::SourceFilter(solutions, matprop, geom_type), 
            G(matprop.get_G()), mg(matprop.get_G())
        {};
                    
        virtual void assign_solutions(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions) {
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
      
      Views::ScalarView** sviews_app;
      Views::VectorView** vviews;
      
      protected:
        void init(unsigned int spn_order, unsigned int G, unsigned int width, unsigned int height, bool display_meshes);
      
      public:
        Visualization(unsigned int spn_order, unsigned int G, bool display_meshes = false)
          : Common::SupportClasses::Visualization(450, 450, display_meshes), mg(G), sviews_app(NULL), vviews(NULL)
        {
          init(spn_order, G, 450, 450, display_meshes);
        }
        Visualization(unsigned int spn_order, unsigned int G, unsigned int width, unsigned int height, bool display_meshes = false)
          : Common::SupportClasses::Visualization(width, height, display_meshes), mg(G), sviews_app(NULL), vviews(NULL)
        {
          init(spn_order, G, width, height, display_meshes);
        }
        virtual ~Visualization();

#ifndef NOGLUT 
        void show_meshes(Hermes::vector<MeshSharedPtr> meshes);
        void show_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions);
        void show_orders(Hermes::vector<SpaceSharedPtr<double> > spaces);
        
        void inspect_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions);

        void show_even_flux_moment(unsigned int moment, unsigned int group, Views::ScalarView* sview,
                                   Hermes::vector< MeshFunctionSharedPtr<double> > solutions);
        void show_odd_flux_moment(unsigned int moment, unsigned int group, Views::VectorView* vview,
                                  Hermes::vector< MeshFunctionSharedPtr<double> > solutions, const MaterialProperties::MaterialPropertyMaps& matprop);
        void show_all_flux_moments(Hermes::vector< MeshFunctionSharedPtr<double> > solutions, const MaterialProperties::MaterialPropertyMaps& matprop);
#else        
        void show_meshes(Hermes::vector<MeshSharedPtr> meshes) {};
        void show_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions) {};
        void show_orders(Hermes::vector<SpaceSharedPtr<double> > spaces) {};
        
        void inspect_solutions(Hermes::vector< MeshFunctionSharedPtr<double> > solutions) {};
        
        void show_even_flux_moment(unsigned int moment, unsigned int group, Views::ScalarView* sview,
                                   Hermes::vector< MeshFunctionSharedPtr<double> > solutions) {};
        void show_odd_flux_moment(unsigned int moment, unsigned int group, Views::VectorView* vview,
                                  Hermes::vector< MeshFunctionSharedPtr<double> > solutions, const MaterialProperties::MaterialPropertyMaps& matprop) {};
        void show_all_flux_moments(Hermes::vector< MeshFunctionSharedPtr<double> > solutions, const MaterialProperties::MaterialPropertyMaps& matprop) {};
#endif
        
        void save_solutions_vtk(const std::string& base_filename, const std::string& base_varname,
                                Hermes::vector< MeshFunctionSharedPtr<double> > solutions, bool mode_3D = false);                   
        void save_orders_vtk(const std::string& base_filename, Hermes::vector<SpaceSharedPtr<double> > spaces);
    };       
    
  /* SupportClasses */
  }
  /* SPN */
  }
  
  namespace SN { namespace SupportClasses 
  {
    class DegreeOrderGroupFlattener
    {
      unsigned int G;
      
      public:
        DegreeOrderGroupFlattener() : G(1) {};
        DegreeOrderGroupFlattener(unsigned int G) : G(G) {};
        
        void set_G(unsigned int G) { this->G = G; }
        unsigned int get_G() const { return G; }
        
        // Accessed array is supposed to contain only elements with even deg+ord. 
        // Progression A002620 used to define the first term (starting position in the linear array for orders and groups 
        // belonging to given degree).
        
        unsigned int pos(unsigned int degree, int order, unsigned int group = 0) const {
          return G*(((degree+1)*(degree+1))/2) + G*(order/2) + group;
        }
        unsigned int operator() (unsigned int degree, int order, unsigned int group = 0) const {
          return G*(((degree+1)*(degree+1))/2) + G*(order/2) + group;
        }
    };
    
    class AngleGroupFlattener
    {
      unsigned int G;
      
      public:
        AngleGroupFlattener() : G(1) {};
        AngleGroupFlattener(unsigned int G) : G(G) {};
        
        void set_G(unsigned int G) { this->G = G; }
        unsigned int get_G() const { return G; }
        
        unsigned int pos(unsigned int angle, unsigned int group) const {
          return angle * G + group;
        }
        unsigned int operator() (unsigned int angle, unsigned int group) const {
          return angle * G + group;
        }
    };
    
    struct OrdinatesData
    {
      OrdinatesData(unsigned int N, const std::string& filename);
      
      template<typename Real>
      void ordinates_to_moment(unsigned int l, int m, unsigned int g, unsigned int G,
                               Func< Real >* *const solution_fns, int num_quad_pts, Real *moment_values_at_quad_pts) const;
      
      template<typename Real>
      void ordinates_to_moment(unsigned int l, int m, unsigned int g, unsigned int G,
                               const Hermes::vector< Real* >& solution_values_at_quad_pts, int num_quad_pts, Real *moment_values_at_quad_pts) const;
      
      unsigned int N, M;
      std::vector<double> xi, eta, mu, pw;
      
      std::vector<int> reflections_about_x;
      std::vector<int> reflections_about_y;
      
      friend std::ostream & operator<< (std::ostream& os, const OrdinatesData& odata);
    };
    
    class SphericalHarmonic
    {
      unsigned int l;
      int m;
      
      double plgndr(double x) const;
      
      int factorial(unsigned int n) const
      {
        int fact=1;
        for (int i=1; i<=n; i++)
          fact*=i;
        return fact;
      }
      
      public:
        SphericalHarmonic(unsigned int l, int m) : l(l), m(m)
        {
          if (m > l)
            Neutronics::ErrorHandling::error_function("Invalid arguments for SphericalHarmonic::operator(): m > l");
        }
        
        double operator() (double xi, double eta, double mu) const;
    };
    
    struct MomentFilter
    {
      class Common 
      {
        protected:
          Common(unsigned int l, int m, unsigned int group, unsigned int G,
                 const OrdinatesData& odata) 
            : l(l), m(m),  g(group), G(G), ag(G), odata(odata)
          {
            if (group >= G) ErrorHandling::error_function("MomentFilter::Common > %s", Messages::E_INVALID_GROUP_INDEX);
          }
          
          unsigned int l, g, G;
          int m;
          AngleGroupFlattener ag;
          const OrdinatesData& odata;
      };
      
      class Val : protected Common, public SimpleFilter<double>
      {
        public:       
          Val(unsigned int l, int m, unsigned int group, unsigned int G, 
              const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
              const OrdinatesData& odata)
            : Common(l, m, group, G, odata), SimpleFilter<double>(solutions, Hermes::vector<int>())
          {
          };
          
          virtual MeshFunction<double>* clone() const;
          
          virtual void set_active_element(Element* e);
          
        protected:             
          void filter_fn(int n, Hermes::vector<double*> values, double* result) {
            odata.ordinates_to_moment<double>(l, m, g, G, values, n, result);
          }
          
      };
      
      class ValDxDy : protected Common, public DXDYFilter<double>
      {
        public:
          ValDxDy(unsigned int l, int m, unsigned int group, unsigned int G,
                  const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
                  const OrdinatesData& odata)
            : Common(l, m, group, G, odata), DXDYFilter<double>(solutions)
          {
          };
          
          virtual MeshFunction<double>* clone() const;
          
          virtual void set_active_element(Element* e);
          
        protected:
          void filter_fn(int n, 
                         Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, 
                         double* rslt, double* rslt_dx, double* rslt_dy) 
          {
            odata.ordinates_to_moment<double>(l, m, g, G, values, n, rslt);
            odata.ordinates_to_moment<double>(l, m, g, G, dx, n, rslt_dx);
            odata.ordinates_to_moment<double>(l, m, g, G, dy, n, rslt_dy);
          };
      };
      
      static void get_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                    Hermes::vector<MeshFunctionSharedPtr<double> >* scalar_fluxes,
                                    unsigned int G, const OrdinatesData& odata);               
      static void get_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                    Hermes::vector<Filter<double>*>* scalar_fluxes,
                                    unsigned int G, const OrdinatesData& odata);  
      static void get_scalar_fluxes_with_derivatives(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                                    Hermes::vector<MeshFunctionSharedPtr<double> >* scalar_fluxes,
                                                    unsigned int G, const OrdinatesData& odata);
      static void get_scalar_fluxes_with_derivatives(const Hermes::vector<MeshFunctionSharedPtr<double> >& angular_fluxes,
                                                    Hermes::vector<Filter<double>*>* scalar_fluxes,
                                                    unsigned int G, const OrdinatesData& odata);
      
      // DEPRECATED (not needed when shared pointers are used for scalar_fluxes) 
      static void clear_scalar_fluxes(Hermes::vector<MeshFunction<double>* >* scalar_fluxes);
      static void clear_scalar_fluxes(Hermes::vector<Filter<double>*>* scalar_fluxes);
    };
    
  /* SupportClasses */
  }
  /* SN */
  }
    
  class PostProcessor
  {
    NeutronicsMethod method;
    GeomType geom_type;
    const SN::SupportClasses::OrdinatesData& odata;
    
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunctionSharedPtr<double> solution,
                                                        const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                        const Hermes::vector<std::string>& regions,
                                                        unsigned int this_group, int other_group = -1) const;
    double get_integrated_group_reaction_rates_internal(ReactionType reaction, MeshFunctionSharedPtr<double> solution,
                                                        const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                        const std::string& region,
                                                        unsigned int this_group, int other_group = -1) const 
    {
      Hermes::vector<std::string> regions; 
      regions.push_back(region);
      return get_integrated_group_reaction_rates_internal(reaction, solution, matprop, regions, this_group, other_group);
    }
    
    public:
      PostProcessor(NeutronicsMethod method, GeomType geom_type = HERMES_PLANAR, const SN::SupportClasses::OrdinatesData* odata = NULL)
        : method(method), geom_type(geom_type), odata(*odata) {};
      
      double integrate(MeshFunctionSharedPtr<double> solution, const Hermes::vector<std::string>& areas = Hermes::vector<std::string>()) const;
      double integrate(MeshFunctionSharedPtr<double> solution, const std::string& area) const 
      {
        Hermes::vector<std::string> areas; 
        areas.push_back(area);
        return integrate(solution, areas);
      }
      
      
      void normalize_to_unit_fission_source(Hermes::vector<MeshFunctionSharedPtr<double> >* solutions, 
                                            double integrated_fission_source) const;
                                            
      void normalize_to_unit_fission_source(Hermes::vector<MeshFunctionSharedPtr<double> >* solutions,
                                            const Common::MaterialProperties::MaterialPropertyMaps& matprop) const;
                                            
      void normalize_to_unit_power(Hermes::vector<MeshFunctionSharedPtr<double> >* solutions,
                                  const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                  double power_per_fission) const;
      
                                  
                                  
      void get_integrated_group_reaction_rates( ReactionType reaction, 
                                                const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, Hermes::vector<double>* results,
                                                const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                unsigned int group, const Hermes::vector<std::string>& regions) const;
                                                
      void get_integrated_group_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, Hermes::vector<double>* results, 
                                              unsigned int group, unsigned int G, 
                                              const Hermes::vector<std::string>& regions) const;
                                              
      void get_integrated_reaction_rates( ReactionType reaction, 
                                          const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, Hermes::vector<double>* results,
                                          const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                          const Hermes::vector<std::string>& regions) const;
                                          
      void get_integrated_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, Hermes::vector<double>* results, 
                                        unsigned int G, const Hermes::vector<std::string>& regions) const;                                                                             
      
      void get_areas(MeshSharedPtr mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const;
                                        
                                        
      double get_integrated_group_reaction_rates( ReactionType reaction, const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
                                                  const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                                  unsigned int group,
                                                  const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                  
      double get_integrated_group_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions, unsigned int group, unsigned int G,
                                                const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                                
      double get_integrated_reaction_rates( ReactionType reaction, const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
                                            const Common::MaterialProperties::MaterialPropertyMaps& matprop,
                                            const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                            
      double get_integrated_scalar_fluxes(const Hermes::vector<MeshFunctionSharedPtr<double> >& solutions,
                                          unsigned int G, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
                                          
      double get_area(MeshSharedPtr mesh, const Hermes::vector<std::string>& regions = Hermes::vector<std::string>()) const;
  };
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}      

#endif
