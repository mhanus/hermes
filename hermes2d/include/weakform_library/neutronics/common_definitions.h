#ifndef ___H2D_NEUTRONICS_COMMON_DEFINITIONS_H
#define ___H2D_NEUTRONICS_COMMON_DEFINITIONS_H

// TODO: Distribute includes to the files that really need them.

#include "weakforms_h1.h"
#include "../forms.h"
#include <sstream>

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Neutronics
    {      
      enum NeutronicsMethod { NEUTRONICS_DIFFUSION, NEUTRONICS_SPN };
      enum ReactionType { 
        ABSORPTION, TOTAL, IN_SCATTERING, SELF_SCATTERING, OUT_SCATTERING, FISSION, NU_FISSION 
      };
      
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
          "In addition, external neutron sources were not specified => expect trivial solution.";
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
      
      namespace DataStructures
      {
        template <typename NDArrayType>
        class StdMultiArray
        {
          private:
            std::vector<NDArrayType> m_data;
          public:
            StdMultiArray(const NDArrayType& val) {
              m_data.push_back(val);
            }
            
            StdMultiArray<NDArrayType>& operator()(const NDArrayType& val) {
              m_data.push_back(val);
              return *this;
            }
            
            operator std::vector<NDArrayType>() {
              return m_data;
            }
        };
        
        template <typename NDArrayType>
        class HermesMultiArray
        {
          private:
            Hermes::vector<NDArrayType> m_data;
          public:
            HermesMultiArray(const NDArrayType& val) {
              m_data.push_back(val);
            }
            
            HermesMultiArray<NDArrayType>& operator()(const NDArrayType& val) {
              m_data.push_back(val);
              return *this;
            }
            
            operator Hermes::vector<NDArrayType>() {
              return m_data;
            }
        };
        
        typedef double rank0;
        typedef std::vector<double> rank1;
        typedef std::vector<std::vector<double > > rank2;
        typedef std::vector<std::vector<std::vector<double > > > rank3;
        
        typedef std::vector<bool > bool1;
        typedef std::vector<std::vector<bool > > bool2;
        typedef std::vector<std::vector<std::vector<bool > > > bool3;
        
        typedef StdMultiArray<rank0> row;
        typedef StdMultiArray<rank1> matrix;
        typedef StdMultiArray<rank2> page;
        typedef StdMultiArray<bool> bool_row;
        typedef StdMultiArray< std::vector<bool> > bool_matrix;
        typedef StdMultiArray< std::vector< std::vector<bool> > > bool_page;
        
        class NDArrayMapOp
        {      
          template <typename NDArrayType>
          static rank0 divide(rank0 x, rank0 y) {
            if (x == 0 && y == 0) 
              return 0.0;
            else if (y == 0)
            {
              error(Messages::E_INF_VALUE);
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
              warning(Messages::W_NEG_VALUE);
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
      }
        
    }
  }
}

#endif 