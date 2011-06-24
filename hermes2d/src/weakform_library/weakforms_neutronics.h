#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "weakforms_h1.h"
#include "../function/forms.h"
#include "../function/filter.h"

namespace WeakFormsNeutronics
{
  namespace Monoenergetic
  {    
    namespace Diffusion 
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
  }
    
  namespace Multigroup
  {
    enum NeutronicsMethod { NEUTRONICS_DIFFUSION, NEUTRONICS_SPN };
    
    namespace MaterialProperties
    {
      namespace Definitions
      {
        typedef double rank0;
        typedef std::vector<double> rank1;
        typedef std::vector<std::vector<double > > rank2;
        typedef std::vector<std::vector<std::vector<double > > > rank3;
        
        typedef std::map<std::string, rank0> MaterialPropertyMap0;
        typedef std::map<std::string, rank1> MaterialPropertyMap1;
        typedef std::map<std::string, rank2> MaterialPropertyMap2;
        typedef std::map<std::string, rank3> MaterialPropertyMap3;
        
        typedef std::map<std::string, std::string> RegionMaterialMap;
        
        typedef std::vector<bool > bool1;
        typedef std::vector<std::vector<bool > > bool2;
        typedef std::vector<std::vector<std::vector<bool > > > bool3;
      }
      
      namespace Messages
      {
        static const char* E_INF_VALUE = 
          "Attempt to set an infinite material property.";
        static const char* W_NEG_VALUE =
          "Entered material data lead to some negative properties.";
        static const char* W_SA_LT_SF =
          "Possible unphysical situation detected: Sigma_a < Sigma_f.";  
        static const char* E_MR_EXTENSION = 
          "Cannot create a multiregion material-property map: no regions specified.";
        static const char* E_INSUFFICIENT_DATA =
          "Not all required material properties have been set.";
        static const char* W_NO_FISSION =
          "Not all required fission properties have been set or could be determined automatically."
          "Assuming a non-fissioning system.";
        static const char* W_NO_SCATTERING =
          "Not all required scattering properties have been set or could be determined automatically."
          "Assuming a purely absorbing system.";
        static const char* E_INVALID_COMBINATION =
          "Invalid combination of entered material properties.";
        static const char* E_NONMATCHING_PROPERTIES =
          "All properties must be defined for a single given number of materials.";
        static const char* E_INVALID_SIZE =
          "Material property defined for an unexpected number of groups.";
        static const char* E_INVALID_GROUP_INDEX =
          "Attempted to access an out-of-range group.";
        static const char* E_INVALID_MARKER =
          "Material data undefined for the given element marker.";
        static const char* E_SG_SIGMA_R = 
          "Group-reduction cross-section (Sigma_r) is not defined for one-group (i.e. monoenergetic) problems."
          "Set Sigma_a instead.";
        static const char* E_LAPACK_ERROR = 
          "Failure of a LAPACK routine with error %d";
        static const char* W_SINGULAR_MATRIX =
          "Entered material properties make the Sigma_rn matrix almost singular. Its inversion will most likely fail.";
        static const char* E_EVEN_SPN =
          "The \"simplified P_N\" model may be used only with odd N.";
        static const char* W_SCATTERING_TRUNCATION =
          "Using SP_%d model - scattering data entered for N > %d will be ignored.";
        static const char* E_MISMATCHED_ORDER_OF_ANISOTROPY =
          "SP_%d model requires Sigma_tn for n = 0,...,%d.";
        static const char* E_SIGMA_T_REQUIRED =
          "SP_%d model requires Sigma_tn for n = 0,...,%d.";
      }
      
      namespace ValidationFunctors
      {
        using namespace Definitions;
        using namespace Messages;
        
        struct ensure_trivial { 
          void operator() (MaterialPropertyMap1::value_type x) { 
            MaterialPropertyMap1::mapped_type::iterator it;
            for (it = x.second.begin(); it != x.second.end(); ++it) 
              if (fabs(*it) > 1e-14)
                error(E_INVALID_COMBINATION);
          }
        };
                
        struct ensure_size { 
          ensure_size(unsigned int nrows, unsigned int ncols = 0, unsigned int npages = 0) 
            : nrows(nrows), ncols(ncols), npages(npages) {};
          
          void operator() (MaterialPropertyMap1::value_type x) { 
            if (x.second.size() != nrows)
              error(E_INVALID_SIZE);
          }
          
          void operator() (MaterialPropertyMap2::value_type x) {
            if (x.second.size() != nrows)
              error(E_INVALID_SIZE);
            
            MaterialPropertyMap2::mapped_type::iterator it;
            for (it = x.second.begin(); it != x.second.end(); ++it) 
              if (it->size() != ncols)
                error(E_INVALID_SIZE);
          }
          
          void operator() (MaterialPropertyMap3::value_type x) {
            if (x.second.size() != npages)
              error(E_MISMATCHED_ORDER_OF_ANISOTROPY, npages);
            
            MaterialPropertyMap3::mapped_type::iterator matrix;
            for (matrix = x.second.begin(); matrix != x.second.end(); ++matrix) 
            {
              if (matrix->size() != nrows)
                error(E_INVALID_SIZE);
              
              rank2::iterator row;
              for (row = matrix->begin(); row != matrix->end(); ++row) 
                if (row->size() != nrows)
                  error(E_INVALID_SIZE);
            }
          }
          
          private:
            unsigned int nrows, ncols, npages;
        };
      }
      
      namespace Common
      {
        using namespace Definitions;
        using namespace Messages;
        
        class NDArrayMapOp
        {
          //
          // NOTE: Could be perhaps combined with the classes material_property_map and MultiArray below
          // and moved to hermes_common as a general way of handling maps with multidimensional mapped types.
          //
          
          template <typename NDArrayType>
          static rank0 divide(rank0 x, rank0 y) {
            if (x == 0 && y == 0) 
              return 0.0;
            else if (y == 0)
            {
              error(E_INF_VALUE);
              return -1.0;
            }
            else
              return x/y;
          }
          
          template <typename NDArrayType>
          static rank0 multiply(rank0 x, rank0 y) {
            return x*y;
          }
          
          template <typename NDArrayType>
          static rank0 add(rank0 x, rank0 y) {
            return x + y;
          }
          
          template <typename NDArrayType>
          static rank0 subtract(rank0 x, rank0 y) {
            rank0 ret = x - y;
            if(ret < 0)
              warning(W_NEG_VALUE);
            return ret;
          }
          
          template <typename NDArrayType>
          static rank0 subtract_nowarn(rank0 x, rank0 y) {
            return x - y;
          }
          
          public: 
            
            #define for_each_element_in_dimension \
                      typedef typename NDArrayType::value_type dim_type;                      \
                      typename NDArrayType::const_iterator dim_iterator_x = x.begin();        \
                      typename NDArrayType::const_iterator dim_iterator_y = y.begin();        \
                      for ( ; dim_iterator_x != x.end(); ++dim_iterator_x, ++dim_iterator_y ) 
            
            template <typename NDArrayType>
            static NDArrayType divide(const NDArrayType& x, const NDArrayType& y)
            { 
              NDArrayType res; res.reserve(x.size());
              
              for_each_element_in_dimension
                res.push_back( divide<dim_type>(*dim_iterator_x, *dim_iterator_y) );
              
              return res;
            }
                        
            template <typename NDArrayType>
            static NDArrayType multiply(const NDArrayType& x, const NDArrayType& y)
            { 
              NDArrayType res; res.reserve(x.size());
              
              for_each_element_in_dimension
                res.push_back( multiply<dim_type>(*dim_iterator_x, *dim_iterator_y) );
              
              return res;
            }
            
            template <typename NDArrayType>
            static NDArrayType add(const NDArrayType& x, const NDArrayType& y)
            { 
              NDArrayType res; res.reserve(x.size());
              
              for_each_element_in_dimension
                res.push_back( add<dim_type>(*dim_iterator_x, *dim_iterator_y) );
              
              return res;
            }

            template <typename NDArrayType>
            static NDArrayType subtract(const NDArrayType& x, const NDArrayType& y)
            { 
              NDArrayType res; res.reserve(x.size());
              
              for_each_element_in_dimension
                res.push_back( subtract<dim_type>(*dim_iterator_x, *dim_iterator_y) );
              
              return res;
            }
            
            template <typename NDArrayType>
            static NDArrayType subtract_nowarn(const NDArrayType& x, const NDArrayType& y)
            { 
              NDArrayType res; res.reserve(x.size());
              
              for_each_element_in_dimension
                res.push_back( subtract_nowarn<dim_type>(*dim_iterator_x, *dim_iterator_y) );
              
              return res;
            }
            
            #undef for_each_element_in_dimension
            
            #define for_each_element_in_map \
                      typename std::map<std::string, T>::iterator iterator_ret = ret.begin();   \
                      typename std::map<std::string, T>::const_iterator iterator_x = x.begin(); \
                      typename std::map<std::string, T>::const_iterator iterator_y = y.begin(); \
                      for ( ; iterator_x != x.end(); ++iterator_x, ++iterator_y, ++iterator_ret ) 
          
            template <typename T>
            static std::map<std::string, T> divide(const std::map<std::string, T>& x, 
                                                   const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = divide<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
            
            template <typename T>
            static std::map<std::string, T> multiply(const std::map<std::string, T>& x, 
                                                     const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = multiply<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
                        
            template <typename T>
            static std::map<std::string, T> add(const std::map<std::string, T>& x, 
                                                const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = add<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
            
            template <typename T>
            static std::map<std::string, T> subtract(const std::map<std::string, T>& x, 
                                                     const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = subtract<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }                                                     
            
            template <typename T>
            static std::map<std::string, T> subtract_nowarn(const std::map<std::string, T>& x, 
                                                            const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = subtract_nowarn<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }        
                                                     
            #undef for_each_element_in_map                                                     
        };
          
        class MaterialPropertyMaps
        {
          protected:
                                
            MaterialPropertyMap1 Sigma_f;
            MaterialPropertyMap1 nu;
            MaterialPropertyMap1 chi;
            
            MaterialPropertyMap1 Sigma_a;
            MaterialPropertyMap1 nuSigma_f;
            
            std::set<std::string> materials_list;
            std::map<std::string, std::string> region_material_map;
            
            unsigned int G;
            
            bool1 fission_nonzero_structure;
                  
            void extend_to_multigroup(const MaterialPropertyMap0& mrsg_map, MaterialPropertyMap1 *mrmg_map);     
            void extend_to_multiregion(const rank1& srmg_array, MaterialPropertyMap1 *mrmg_map);
            void extend_to_multiregion_multigroup(const rank0& srsg_value, MaterialPropertyMap1 *mrmg_map);
            
            MaterialPropertyMap1 extract_map2_diagonals(const MaterialPropertyMap2& map2) const;
            
            MaterialPropertyMap1 sum_map2_columns(const MaterialPropertyMap2& map2) const;
            MaterialPropertyMap1 sum_map2_rows(const MaterialPropertyMap2& map2) const;
            
            MaterialPropertyMap2 create_map2_by_diagonals(const MaterialPropertyMap1& diags) const;
            
            void fill_with(double c, MaterialPropertyMap1 *mrmg_map);
            void fill_with(double c, MaterialPropertyMap2 *mrmg_map);
                      
            MaterialPropertyMaps(unsigned int G, const std::set<std::string>& mat_list = std::set<std::string>()) 
              : materials_list(mat_list), G(G)  { };
              
            MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map);
                        
            virtual void validate();
            
          public:
            
            std::string get_material(int elem_marker, WeakForm *wf) const;
            std::string get_material(int elem_marker, Mesh *mesh) const;
            
            virtual void set_nu(const MaterialPropertyMap1& nu) {
              this->nu = nu;
            }
            
            virtual void set_nu(const MaterialPropertyMap0& nu) {
              extend_to_multigroup(nu, &this->nu);      
            }
            
            virtual void set_nu(const rank1& nu) {
              extend_to_multiregion(nu, &this->nu);
            }
            
            virtual void set_nu(const rank0& nu) {
              extend_to_multiregion_multigroup(nu, &this->nu);
            }
                  
            virtual void set_chi(const MaterialPropertyMap1& chi) {
              this->chi = chi;
            }
            
            virtual void set_chi(const rank1& chi) {
              extend_to_multiregion(chi, &this->chi);
            }
                        
            virtual void set_Sigma_a(const MaterialPropertyMap1& Sa) {
              this->Sigma_a = Sa;
            }
            
            virtual void set_Sigma_f(const MaterialPropertyMap1& Sf) {
              this->Sigma_f = Sf;
            }
            
            virtual void set_nuSigma_f(const MaterialPropertyMap1 nSf) {
              this->nuSigma_f = nSf;
            }
            
            virtual void set_materials_list(std::set<std::string> mat_list) {
              this->materials_list = mat_list;
            }
            
            const MaterialPropertyMap1& get_Sigma_f() const {
              return this->Sigma_f;
            }
            const MaterialPropertyMap1& get_nu() const {
              return this->nu;
            }
            const MaterialPropertyMap1& get_chi() const {
              return this->chi;
            }
            const bool1& get_fission_nonzero_structure() const {
              return this->fission_nonzero_structure;
            }
            const std::set<std::string>& get_materials_list() const {
              return this->materials_list;
            }
            
            const rank1& get_Sigma_f(const std::string& material) const;
            const rank1& get_nu(const std::string& material) const;
            const rank1& get_chi(const std::string& material) const;
            
            unsigned int get_G() const { return G; } 
            
            friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
        };
      }
      
      namespace Diffusion
      {
        using namespace Definitions;
        using namespace Messages;
        
        class MaterialPropertyMaps : public Common::MaterialPropertyMaps
        {
          protected:
            
            MaterialPropertyMap1 D;
            MaterialPropertyMap1 Sigma_r;
            MaterialPropertyMap2 Sigma_s;
            
            MaterialPropertyMap1 src;
            
            MaterialPropertyMap1 Sigma_t;
            
            bool2 scattering_nonzero_structure;
            
          public:
            
            MaterialPropertyMaps(unsigned int G, const std::set<std::string>& mat_list = std::set<std::string>()) 
              : Common::MaterialPropertyMaps(G, mat_list) { };
            MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map)
              : Common::MaterialPropertyMaps(G, reg_mat_map) { };
            
            // We always need to supply chi, nu, Sigma_f, Sigma_r, Sigma_s and D to our neutronics weak forms. 
            // These parameters are often defined in terms of the other ones, or not specified at all and assumed 
            // to be zero for a particular simplified situation. This method, together with its complement in the
            // parent class, uses the most typical definitions to build the six-parameter set from the given input. 
            // It also checks whether the user did not enter nonsensical values. However, values entered by the 
            // user may sometimes not satisfy the common relations, as some empirical corrections may have been 
            // already included in them.
            virtual void validate();
            
            virtual void set_src(const MaterialPropertyMap1& src) {
              this->src = src;
            }
            
            virtual void set_src(const MaterialPropertyMap0& src) {
              extend_to_multigroup(src, &this->src);            
            }
            
            virtual void set_src(const rank1& src) {
              extend_to_multiregion(src, &this->src);
            }
            
            virtual void set_src(const double& src) {
              extend_to_multiregion_multigroup(src, &this->src);
            }
            
            virtual void set_D(const MaterialPropertyMap1& D) {
              this->D = D;
            }
            
            virtual void set_Sigma_r(const MaterialPropertyMap1& Sr) {
              this->Sigma_r = Sr;
            }
            
            virtual void set_Sigma_t(const MaterialPropertyMap1& St) {
              this->Sigma_t = St;
            }
            
            virtual void set_Sigma_s(const MaterialPropertyMap2& Ss) {
              this->Sigma_s = Ss;
            }
                         
            const MaterialPropertyMap2& get_Sigma_s() const {
              return this->Sigma_s;
            }
            const MaterialPropertyMap1& get_Sigma_r() const {
              return this->Sigma_r;
            }
            const MaterialPropertyMap1& get_D() const {
              return this->D;
            }
            const MaterialPropertyMap1& get_src() const {
              return this->src;
            }
            const bool2& get_scattering_nonzero_structure() const {
              return this->scattering_nonzero_structure;
            }
            
            const rank2& get_Sigma_s(const std::string& material) const;
            const rank1& get_Sigma_r(const std::string& material) const;
            const rank1& get_D(const std::string& material) const;
            const rank1& get_src(const std::string& material) const;
            
            friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
        };
        
        class TransportCorrectedMaterialPropertyMaps : public MaterialPropertyMaps
        {
          protected:
            
            MaterialPropertyMap1 Sigma_s_1_out;
            MaterialPropertyMap1 mu_av;
          
          public:
            
            TransportCorrectedMaterialPropertyMaps(unsigned int G, 
                                                   const MaterialPropertyMap2& Ss_1,
                                                   const RegionMaterialMap& reg_mat_map = RegionMaterialMap());                                            
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   const MaterialPropertyMap1& mu_av,
                                                   const RegionMaterialMap& reg_mat_map = RegionMaterialMap());
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   const MaterialPropertyMap0& mu_av,
                                                   const RegionMaterialMap& reg_mat_map = RegionMaterialMap());                                                                              
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   const RegionMaterialMap& reg_mat_map,
                                                   const rank1& mu_av);
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   const RegionMaterialMap& reg_mat_map,
                                                   const rank0& mu_av);                                       
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   std::set<std::string> mat_list,
                                                   const rank1& mu_av);
            TransportCorrectedMaterialPropertyMaps(unsigned int G,
                                                   std::set<std::string> mat_list,
                                                   const rank0& mu_av);
                                                   
            virtual void set_D(const MaterialPropertyMap1& D) {
              warning("Diffusion coefficient is determined automatically according to the transport correction formula."
                      " Ignoring user setting.");
            }
            
            virtual void validate();
        };
      }  
      
      namespace SPN
      {
        using namespace Definitions;
        using namespace Messages;
        
        extern "C" 
        {
          // LU decomoposition of a general matrix.
          void dgetrf_(int* m, int *n, double* A, int* lda, int* ipiv, int* info);
          
          // Generate inverse of a matrix given its LU decomposition.
          void dgetri_(int* n, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info);
          
          // Compute the norm of A.
          double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work);
          
          // Compute the reciprocal of the condition number of A.
          int dgecon_(const char *norm, const int *n, double *a, const int *lda, const double *anorm, double *rcond, 
                      double *work, int *iwork, int *info);
        }
        
        class MaterialPropertyMaps : public Common::MaterialPropertyMaps
        {
          protected:
            
            MaterialPropertyMap3 Sigma_rn;
            MaterialPropertyMap3 odd_Sigma_rn_inv;
            MaterialPropertyMap1 src0;
            
            MaterialPropertyMap3 Sigma_sn;
            MaterialPropertyMap3 Sigma_tn;
            
            std::map<std::string, bool1> Sigma_rn_is_diagonal;
            
            void extend_to_rank3(const MaterialPropertyMap2& src, MaterialPropertyMap3* dest);
            void extend_to_rank3(const MaterialPropertyMap1& src, MaterialPropertyMap3* dest);
            MaterialPropertyMap3 create_map3_by_diagonals(const MaterialPropertyMap2& diags) const;
            void fill_with(double c, MaterialPropertyMap3 *mmmrmg_map);
            void invert_odd_Sigma_rn();
            
            unsigned int N, N_odd;
            
          public:
            
            MaterialPropertyMaps(unsigned int G, unsigned int N,
                                 std::set<std::string> mat_list = std::set<std::string>()) 
              : Common::MaterialPropertyMaps(G, mat_list), N(N), N_odd((N+1)/2) 
            {
              if ((N % 2) == 0) error(E_EVEN_SPN);
            }
            
            MaterialPropertyMaps(unsigned int G, unsigned int N, 
                                 const RegionMaterialMap& reg_mat_map)
              : Common::MaterialPropertyMaps(G, reg_mat_map), N(N), N_odd((N+1)/2)  
            { 
              if ((N % 2) == 0) error(E_EVEN_SPN);
            }
                                    
            virtual void validate();
            
            virtual void set_src0(const MaterialPropertyMap1& src) {
              this->src0 = src;
            }
            
            virtual void set_src(const MaterialPropertyMap0& src) {
              extend_to_multigroup(src, &this->src0);            
            }
            
            virtual void set_src(const rank1& src) {
              extend_to_multiregion(src, &this->src0);
            }
            
            virtual void set_src(const double& src) {
              extend_to_multiregion_multigroup(src, &this->src0);
            }
                        
            virtual void set_Sigma_rn(const MaterialPropertyMap3& Sr) {
              this->Sigma_rn = Sr;
            }
                        
            virtual void set_Sigma_sn(const MaterialPropertyMap3& Ss) {
              this->Sigma_sn = Ss;
            }
            
            virtual void set_Sigma_tn(const MaterialPropertyMap3& St) {
              this->Sigma_tn = St;
            }
            
            virtual void set_Sigma_tn(const MaterialPropertyMap2& St) {              
              this->Sigma_tn = create_map3_by_diagonals(St);
            }
            
            virtual void set_Sigma_tn(const MaterialPropertyMap1& St) {
              extend_to_rank3(St, &this->Sigma_tn);
            }
                         
            const MaterialPropertyMap3& get_Sigma_sn() const {
              return this->Sigma_sn;
            }
            const MaterialPropertyMap3& get_Sigma_rn() const {
              return this->Sigma_rn;
            }
            const MaterialPropertyMap1& get_src0() const {
              return this->src0;
            }
            
            const bool1  is_Sigma_rn_diagonal() const;
            
            const bool1& is_Sigma_rn_diagonal(const std::string& material) const;
            const rank3& get_Sigma_rn(const std::string& material) const;
            const rank3& get_odd_Sigma_rn_inv(const std::string& material) const;
            const rank1& get_src0(const std::string& material) const;
            
            unsigned int get_N() const { return N; }
            unsigned int get_N_odd() const { return N_odd; }
            
            friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
        };
      }  
      
      template <typename NDArrayType>
      class material_property_map
      {
        private:
          std::map<std::string, NDArrayType> m_map;
        public:
          material_property_map(const std::string& key, const NDArrayType& val) {
            m_map[key] = val;
          }
          
          material_property_map<NDArrayType>& operator()(const std::string& key, const NDArrayType& val) {
            m_map[key] = val;
            return *this;
          }
          
          operator std::map<std::string, NDArrayType>() {
            return m_map;
          }
      };
      
      class region_material_map
      {
        private:
          Definitions::RegionMaterialMap m_map;
        public:
          region_material_map(const std::string& key, const std::string& val) {
            m_map[key] = val;
          }
          
          region_material_map& operator()(const std::string& key, const std::string& val) {
            m_map[key] = val;
            return *this;
          }
          
          operator Definitions::RegionMaterialMap() {
            return m_map;
          }
      };
      
      template <typename NDArrayType>
      class MultiArray
      {
        private:
          std::vector<NDArrayType> m_data;
        public:
          MultiArray(const NDArrayType& val) {
            m_data.push_back(val);
          }
          
          MultiArray<NDArrayType>& operator()(const NDArrayType& val) {
            m_data.push_back(val);
            return *this;
          }
          
          operator std::vector<NDArrayType>() {
            return m_data;
          }
      };
      
      namespace Definitions
      {
        typedef MultiArray<rank0> row;
        typedef MultiArray<rank1> matrix;
        typedef MultiArray<rank2> page;
        typedef MultiArray<bool> bool_row;
        typedef MultiArray< std::vector<bool> > bool_matrix;
        typedef MultiArray< std::vector< std::vector<bool> > > bool_page;
      }
    }
    
    namespace SupportClasses
    {
      namespace Common
      {
        using namespace MaterialProperties::Definitions;
        
        class SourceFilter : public SimpleFilter
        {
          public: 
            SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                         const std::vector<std::string>& source_regions = std::vector<std::string>())
              : SimpleFilter(), matprop(matprop),
                source_regions(source_regions.begin(), source_regions.end())
            {
              pre_init();
            };
            SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                         const std::string& source_region)
              : SimpleFilter(), matprop(matprop)
            { 
              source_regions.insert(source_region);
              pre_init();
            }
            
            virtual void assign_solutions(const Hermes::vector<Solution*>& solutions);
            virtual void assign_solutions(const Hermes::vector<MeshFunction*>& solutions);
            
            double integrate();
                        
          protected:
            const MaterialProperties::Common::MaterialPropertyMaps& matprop;
            std::set<std::string> source_regions;
            std::set<int> markers;
            bool have_solutions;
            
            virtual void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);
            virtual void pre_init();
            virtual void post_init();
        };
      }
      
      namespace SPN
      {
        using namespace MaterialProperties::Definitions;
        
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
          
          class Val : protected Common, public SimpleFilter
          {
            public:       
              Val(unsigned int angular_moment, unsigned int group, unsigned int G, 
                  const Hermes::vector<MeshFunction*>& solutions)
                : Common(angular_moment, group, G), SimpleFilter(solutions, Hermes::vector<int>())
              {};
              Val(unsigned int angular_moment, unsigned int group, unsigned int G,
                  const Hermes::vector<Solution*>& solutions)
                : Common(angular_moment, group, G), SimpleFilter(solutions, Hermes::vector<int>())
              {};
              
            protected:             
              void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);
          };
          
          class ValDxDy : protected Common, public DXDYFilter
          {
            public:
              ValDxDy(unsigned int angular_moment, unsigned int group, unsigned int G, 
                      const Hermes::vector<MeshFunction*>& solutions)
                : Common(angular_moment, group, G), DXDYFilter(solutions)
              {};
              ValDxDy(unsigned int angular_moment, unsigned int group, unsigned int G,
                      const Hermes::vector<Solution*>& solutions)
                : Common(angular_moment, group, G), DXDYFilter(solutions)
              {};
              
            protected:
              void filter_fn(int n, 
                             Hermes::vector<scalar *> values, Hermes::vector<scalar *> dx, Hermes::vector<scalar *> dy, 
                             scalar* rslt, scalar* rslt_dx, scalar* rslt_dy);
          };
          
          static void get_scalar_fluxes(const Hermes::vector<Solution*>& angular_fluxes,
                                        Hermes::vector<MeshFunction*>* scalar_fluxes,
                                        unsigned int G);
          static void get_scalar_fluxes_with_derivatives(const Hermes::vector<Solution*>& angular_fluxes,
                                                         Hermes::vector<MeshFunction*>* scalar_fluxes,
                                                         unsigned int G);
          static void clear_scalar_fluxes(Hermes::vector<MeshFunction*>* scalar_fluxes);
        };
        
        class SourceFilter : public Common::SourceFilter
        {
          public: 
            SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                         const std::vector<std::string>& source_regions = std::vector<std::string>())
              : Common::SourceFilter(matprop, source_regions), G(matprop.get_G()), mg(G)
            {};
            SourceFilter(const MaterialProperties::Common::MaterialPropertyMaps& matprop, const std::string& source_region)
              : Common::SourceFilter(matprop, source_region), G(matprop.get_G()), mg(G) 
            {};
                        
            virtual void assign_solutions(const Hermes::vector<Solution*>& solutions) {
              num = solutions.size();
              Common::SourceFilter::assign_solutions(solutions);
            }
            virtual void assign_solutions(const Hermes::vector<MeshFunction*>& solutions) {
              num = solutions.size();
              Common::SourceFilter::assign_solutions(solutions);
            }
            
          protected:
            unsigned int G;
            MomentGroupFlattener mg;
            
            virtual void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);
        };
      }
    }
                                 
    namespace ElementaryForms
    {             
      namespace Diffusion
      { 
        using namespace MaterialProperties::Diffusion;
        
        class GenericForm
        {
          protected:
            const MaterialPropertyMaps& matprop;
            GeomType geom_type;
            
            GenericForm(const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR)
              : matprop(matprop), geom_type(geom_type) 
            {};
        };
        
        struct VacuumBoundaryCondition
        {
          // TODO: General albedo boundary condition.
          class Jacobian : public WeakForm::MatrixFormSurf
          {
            public:
              Jacobian(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::MatrixFormSurf(g,g,HERMES_ANY), 
                g(g), geom_type(geom_type)
              {};
              
              Jacobian(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::MatrixFormSurf(g,g,area),
                g(g), geom_type(geom_type)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                  Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
                return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
                return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormSurf* clone() {
                return new Jacobian(*this);
              }
                            
            private:
              unsigned int g;
              GeomType geom_type;
          };
          
          class Residual : public WeakForm::VectorFormSurf
          {
            public:
              Residual(unsigned int g, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::VectorFormSurf(g,HERMES_ANY), 
                g(g), geom_type(geom_type)
              {};
              
              Residual(unsigned int g, const std::string& area, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::VectorFormSurf(g,area),
                g(g), geom_type(geom_type)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                 Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormSurf* clone() {
                return new Residual(*this);
              }
                            
            private:
              unsigned int g;
              GeomType geom_type;
          };
        };
        
        struct DiffusionReaction
        {   
          class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
          {
            public:            
              Jacobian(unsigned int g, 
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::MatrixFormVol(g, g, HERMES_ANY, HERMES_SYM),
                  GenericForm(matprop, geom_type),
                  g(g)
              {};
                  
              Jacobian(unsigned int g, const std::string& area,
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : WeakForm::MatrixFormVol(g, g, area, HERMES_SYM),
                  GenericForm(matprop, geom_type),
                  g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;

              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }

              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }

            private:
              
              unsigned int g;
          };
          
          class Residual : public WeakForm::VectorFormVol, protected GenericForm
          {
            public:
              
              Residual(unsigned int g, 
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : WeakForm::VectorFormVol(g, HERMES_ANY),
                  GenericForm(matprop, geom_type),
                  g(g)
              {};
                  
              Residual(unsigned int g, const std::string& area,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : WeakForm::VectorFormVol(g, area),
                  GenericForm(matprop, geom_type), 
                  g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const  {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const  {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int g;
          };
        };
      
        struct FissionYield
        {
          class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
          {
            public:
              
              Jacobian( unsigned int gto, unsigned int gfrom, 
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::MatrixFormVol(gto, gfrom), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {};
              
              Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::MatrixFormVol(gto, gfrom, area), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  -1.0 * matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }
              
            private:
              
              unsigned int gto, gfrom;
          };
      
          class OuterIterationForm : public WeakForm::VectorFormVol, protected GenericForm
          {
            public:
              
              OuterIterationForm( unsigned int g, 
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(g, HERMES_ANY, iterates),
                  GenericForm(matprop, geom_type),
                  g(g), keff(keff)
              {
                if (g >= iterates.size())
                  error(E_INVALID_GROUP_INDEX);
              }
              
              OuterIterationForm( unsigned int g, const std::string& area,
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(g, area, iterates),
                  GenericForm(matprop, geom_type),
                  g(g), keff(keff)
              {
                if (g >= iterates.size())
                  error(E_INVALID_GROUP_INDEX);
              }
              
              OuterIterationForm( unsigned int g, const Hermes::vector<std::string>& areas,
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(g, areas, iterates),
                  GenericForm(matprop, geom_type),
                  g(g), keff(keff)
              {
                if (g >= iterates.size())
                  error(E_INVALID_GROUP_INDEX);
              }
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                  Geom<double> *e, ExtData<scalar> *ext) const {
                return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const  {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }

              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new OuterIterationForm(*this);
              }
              
              void update_keff(double new_keff) { keff = new_keff; }
              
            private:
              
              unsigned int g;
              double keff;
          };
        
          class Residual : public WeakForm::VectorFormVol, protected GenericForm
          {
            public:
              Residual( unsigned int gto, unsigned int gfrom, 
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(gto), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {};
              
              Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(gto, area), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int gto, gfrom;
          };
        };
  
        struct Scattering
        {      
          class Jacobian : public WeakForm::MatrixFormVol, protected GenericForm
          {
            public:
              
              Jacobian( unsigned int gto, unsigned int gfrom, 
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::MatrixFormVol(gto, gfrom), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {
                this->scaling_factor = (gto != gfrom) ? -1 : 0;
              };
              
              Jacobian( unsigned int gto, unsigned int gfrom, const std::string& area,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : WeakForm::MatrixFormVol(gto, gfrom, area), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {
                this->scaling_factor = (gto != gfrom) ? -1 : 0;
              };
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }
              
            private:
              
              unsigned int gto, gfrom;
          };
        
          class Residual : public WeakForm::VectorFormVol, protected GenericForm
          {
            public:
              Residual( unsigned int gto, unsigned int gfrom, 
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(gto), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {
                this->scaling_factor = (gto != gfrom) ? -1 : 0;
              };
              
              Residual( unsigned int gto, unsigned int gfrom, const std::string& area,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : WeakForm::VectorFormVol(gto, area), 
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
              {
                this->scaling_factor = (gto != gfrom) ? -1 : 0;
              };
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int gto, gfrom;
          };
        };
        
        struct ExternalSources
        {
          class LinearForm : public WeakForm::VectorFormVol, protected GenericForm
          {
            public:
              
              LinearForm( unsigned int g, 
                          const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : WeakForm::VectorFormVol(g), 
                  GenericForm(matprop, geom_type),
                  g(g)
              {};
              
              LinearForm( unsigned int g, const std::string& area,
                          const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : WeakForm::VectorFormVol(g, area), 
                  GenericForm(matprop, geom_type),
                  g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
                              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new LinearForm(*this);
              }
            
            private:
              
              unsigned int g;      
          }; 
        };
            
      }                
    
      namespace SPN
      {
        using namespace MaterialProperties::SPN;
        using SupportClasses::SPN::Coeffs;
        using SupportClasses::SPN::MomentGroupFlattener;
        
        //TODO: Make Diffusion::GenericForm only a GenericForm, which takes pointer to
        // MaterialProperties::Common::MaterialPropertyMaps (hence all "matprop." will 
        // have to be changed to "matprop->". The following class will then not be needed.
        class GenericForm
        {
          protected:
            const MaterialPropertyMaps& matprop;
            GeomType geom_type;
            MomentGroupFlattener mg;
            
            GenericForm(const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR)
              : matprop(matprop), geom_type(geom_type), mg(matprop.get_G())
            {};
        };
        
        struct VacuumBoundaryCondition
        {
          // TODO: General albedo boundary condition.
          class Jacobian : protected GenericForm, public WeakForm::MatrixFormSurf
          {
            public:
              Jacobian(unsigned int m, unsigned int n, unsigned int g, 
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormSurf(mg.pos(m,g),mg.pos(n,g),HERMES_ANY),
                  mrow(m), mcol(n), g(g)
              {};
              
              Jacobian(unsigned int m, unsigned int n,  unsigned int g, const std::string& area, 
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormSurf(mg.pos(m,g),mg.pos(n,g),area),
                  mrow(m), mcol(n), g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext ) const {
                return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
                return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormSurf* clone() {
                return new Jacobian(*this);
              }
                            
            private:
              
              unsigned int mrow, mcol;
              unsigned int g;
          };
          
          class Residual : protected GenericForm, public WeakForm::VectorFormSurf
          {
            public:
              Residual(unsigned int m, unsigned int N, unsigned int g,
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormSurf(mg.pos(m,g),HERMES_ANY), 
                  mrow(m), N_odd((N+1)/2), g(g)
              {};
              
              Residual(unsigned int m, unsigned int N, unsigned int g, const std::string& area, 
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormSurf(mg.pos(m,g),area), 
                  mrow(m), N_odd((N+1)/2), g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                 Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormSurf* clone() {
                return new Residual(*this);
              }
                            
            private:
              unsigned int mrow;
              unsigned int N_odd;
              unsigned int g;
          };
        };
        
        struct DiagonalStreamingAndReactions
        {   
          class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
          {
            public:            
              Jacobian(unsigned int m, unsigned int g, 
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,g), mg.pos(m,g), HERMES_ANY, HERMES_SYM),
                  mrow(m), g(g)
              {};
                  
              Jacobian(unsigned int m, unsigned int g, const std::string& area,
                       const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,g), mg.pos(m,g), area, HERMES_SYM),
                  mrow(m), g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;

              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }

              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }

            private:
              
              unsigned int mrow;
              unsigned int g;
          };
          
          class Residual : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              
              Residual(unsigned int m, unsigned int g, 
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR) 
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g), HERMES_ANY),
                  mrow(m), g(g)
              {};
                  
              Residual(unsigned int m, unsigned int g, const std::string& area,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : GenericForm(matprop, geom_type), 
                  WeakForm::VectorFormVol(mg.pos(m,g), area),
                  mrow(m), g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const  {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const  {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int mrow;
              unsigned int g;
          };
        };
      
        struct FissionYield
        {
          class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
          {
            public:
              
              Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
                  mrow(m), mcol(n), gto(gto), gfrom(gfrom)
              {};
              
              Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                        const std::string& area, const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
                  mrow(m), mcol(n), gto(gto), gfrom(gfrom)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }
              
            private:
              
              unsigned int mrow, mcol;
              unsigned int gto, gfrom;
          };
      
          class OuterIterationForm : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              
              OuterIterationForm( unsigned int m, unsigned int g,
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g), HERMES_ANY, iterates),
                  mrow(m), g(g), keff(keff)
              {};
              
              OuterIterationForm( unsigned int m, unsigned int g, const std::string& area,
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g), area, iterates),
                  mrow(m), g(g), keff(keff)
              {};
              
              OuterIterationForm( unsigned int m, unsigned int g, const Hermes::vector<std::string>& areas,
                                  const MaterialPropertyMaps& matprop,
                                  const Hermes::vector<MeshFunction*>& iterates,
                                  double keff = 1.0,
                                  GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g), areas, iterates),
                  mrow(m), g(g), keff(keff)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                  Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const  {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }

              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new OuterIterationForm(*this);
              }
              
              void update_keff(double new_keff) { keff = new_keff; }
              
            private:
              
              unsigned int mrow;
              unsigned int g;
              double keff;
          };
        
          class Residual : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                        const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto)), 
                  mrow(m), N_odd(N_odd), gto(gto)
              {};
              
              Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                        const std::string& area, const MaterialPropertyMaps& matprop, 
                        GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto), area), 
                  mrow(m), N_odd(N_odd), gto(gto)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int mrow;
              unsigned int N_odd;
              unsigned int gto;
          };
        };
        
        struct OffDiagonalStreaming
        {      
          class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
          {
            public:
              
              Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(m,gfrom), HERMES_ANY, sym), 
                  mrow(m), gto(gto), gfrom(gfrom)
              {};
              
              Jacobian( unsigned int m, unsigned int gto, unsigned int gfrom,
                        const std::string& area, const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(m,gfrom), area, sym), 
                  mrow(m), gto(gto), gfrom(gfrom)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }
              
            private:
              
              unsigned int mrow;
              unsigned int gto, gfrom;
          };
        
          class Residual : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              Residual( unsigned int m, unsigned int gto,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto)), 
                  mrow(m), gto(gto)
              {};
              
              Residual( unsigned int m, unsigned int gto,
                        const std::string& area, const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto), area), 
                  mrow(m), gto(gto)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int mrow;
              unsigned int gto;
          };
        };
        
        struct OffDiagonalReactions
        {      
          class Jacobian : protected GenericForm, public WeakForm::MatrixFormVol
          {
            public:
              
              Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), HERMES_ANY, sym), 
                  mrow(m), mcol(n), gto(gto), gfrom(gfrom)
              {};
              
              Jacobian( unsigned int m, unsigned int n, unsigned int gto, unsigned int gfrom,
                        const std::string& area, const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR, SymFlag sym = HERMES_NONSYM )
                : GenericForm(matprop, geom_type),
                  WeakForm::MatrixFormVol(mg.pos(m,gto), mg.pos(n,gfrom), area, sym), 
                  mrow(m), mcol(n), gto(gto), gfrom(gfrom)
              {};
              
              template<typename Real, typename Scalar>
              Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const;
              
              virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const { 
                return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
              }
              
              virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const { 
                return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::MatrixFormVol* clone() {
                return new Jacobian(*this);
              }
              
            private:
              
              unsigned int mrow, mcol;
              unsigned int gto, gfrom;
          };
        
          class Residual : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                        const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto)), 
                  mrow(m), N_odd(N_odd), gto(gto)
              {};
              
              Residual( unsigned int m, unsigned int N_odd, unsigned int gto,
                        const std::string& area, const MaterialPropertyMaps& matprop,
                        GeomType geom_type = HERMES_PLANAR )
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,gto), area), 
                  mrow(m), N_odd(N_odd), gto(gto)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new Residual(*this);
              }
              
            private:
              
              unsigned int mrow;
              unsigned int N_odd;
              unsigned int gto;
          };
        };
        
        struct ExternalSources
        {
          class LinearForm : protected GenericForm, public WeakForm::VectorFormVol
          {
            public:
              
              LinearForm( unsigned int m, unsigned int g,
                          const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g)), 
                  mrow(m), g(g)
              {};
              
              LinearForm( unsigned int m, unsigned int g, const std::string& area,
                          const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                : GenericForm(matprop, geom_type),
                  WeakForm::VectorFormVol(mg.pos(m,g), area), 
                  mrow(m), g(g)
              {};
              
              template<typename Real, typename Scalar>
              Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
              
              virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const {
                return -1.0 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
              }

              virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const {
                return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
              }
                              
              // This is to make the form usable in rk_time_step().
              virtual WeakForm::VectorFormVol* clone() {
                return new LinearForm(*this);
              }
            
            private:
              
              unsigned int mrow;
              unsigned int g;      
          }; 
        };
      }
    }
    
    namespace CompleteWeakForms
    {
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
        using namespace MaterialProperties;
        using namespace ElementaryForms::Diffusion;
               
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
        using namespace MaterialProperties;
        using namespace ElementaryForms::SPN;
        
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
          public:
            DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                       GeomType geom_type = HERMES_PLANAR);
            
            DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                       HermesFunction *minus_f_src,
                                       std::string src_area = HERMES_ANY,
                                       GeomType geom_type = HERMES_PLANAR);
            
            DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                       HermesFunction *minus_f_src,
                                       Hermes::vector<std::string> src_areas,
                                       GeomType geom_type = HERMES_PLANAR);
            
            DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                       const std::vector<HermesFunction*>& minus_f_src,
                                       std::string src_area = HERMES_ANY,
                                       GeomType geom_type = HERMES_PLANAR);
            
            DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                       const std::vector<HermesFunction*>& minus_f_src,
                                       Hermes::vector<std::string> src_areas,
                                       GeomType geom_type = HERMES_PLANAR);
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
    }
    
    namespace SupportClasses
    {
      class SourceIteration
      {
        const Hermes2D& hermes2d;
        DiscreteProblem& dp;
        const std::vector<std::string>& fission_regions;
        
        CompleteWeakForms::Common::WeakFormSourceIteration *wf;
        
        Common::SourceFilter *new_source, *old_source;
        
        public:
          /// \param[in] fission_regions  Strings specifiying the parts of the solution domain where fission occurs.
          /// \param[in]     hermes2d     Class encapsulating global Hermes2D functions.
          /// \param[in]     spaces       Pointers to spaces on which the solutions are defined (one space for each energy group).
          /// \param[in]     wf           Pointer to the weak form of the problem.
          SourceIteration(NeutronicsMethod method, const MaterialProperties::Common::MaterialPropertyMaps& matprop,
                          const std::vector<std::string>& fission_regions, 
                          const Hermes2D& hermes2d, DiscreteProblem& dp);
                          
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
          int eigenvalue_iteration(const Hermes::vector<Solution *>& solutions, 
                                   double tol = 1e-6, MatrixSolverType matrix_solver = SOLVER_UMFPACK);
                              
          double integrate_over_fission_regions(MeshFunction* sln);
      };
    }
  }
}
#endif
