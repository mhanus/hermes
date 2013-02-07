//
//TODO: MaterialRegionMap
//

#ifndef ___H2D_NEUTRONICS_MATERIAL_PROPERTIES_H
#define ___H2D_NEUTRONICS_MATERIAL_PROPERTIES_H

#include "common_definitions.h"
#include "weakform.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{  
  typedef std::map<std::string, rank0> MaterialPropertyMap0;
  typedef std::map<std::string, rank1> MaterialPropertyMap1;
  typedef std::map<std::string, rank2> MaterialPropertyMap2;
  typedef std::map<std::string, rank3> MaterialPropertyMap3;
  
  typedef std::map<std::string, std::string> RegionMaterialMap;
  
  namespace Common { namespace MaterialProperties
  {    
    struct Validation
    {
      struct ensure_trivial { 
        void operator() (MaterialPropertyMap1::value_type x) { 
          MaterialPropertyMap1::mapped_type::iterator it;
          for (it = x.second.begin(); it != x.second.end(); ++it) 
            if (fabs(*it) > 1e-14)
              ErrorHandling::error_function(Messages::E_INVALID_COMBINATION);
        }
      };
      
      struct ensure_size { 
        ensure_size(unsigned int nrows, unsigned int ncols = 0, unsigned int npages = 0) 
        : nrows(nrows), ncols(ncols), npages(npages) {};
        
        void operator() (MaterialPropertyMap1::value_type x) { 
          if (x.second.size() != nrows)
            ErrorHandling::error_function(Messages::E_INVALID_SIZE);
        }
        
        void operator() (MaterialPropertyMap2::value_type x) {
          if (x.second.size() != nrows)
            ErrorHandling::error_function(Messages::E_INVALID_SIZE);
          
          MaterialPropertyMap2::mapped_type::iterator it;
          for (it = x.second.begin(); it != x.second.end(); ++it) 
            if (it->size() != ncols)
              ErrorHandling::error_function(Messages::E_INVALID_SIZE);
        }
        
        void operator() (MaterialPropertyMap3::value_type x) {
          if (x.second.size() != npages)
            ErrorHandling::error_function(Messages::E_MISMATCHED_ORDER_OF_ANISOTROPY, npages);
          
          MaterialPropertyMap3::mapped_type::iterator matrix;
          for (matrix = x.second.begin(); matrix != x.second.end(); ++matrix) 
          {
            if (matrix->size() != nrows)
              ErrorHandling::error_function(Messages::E_INVALID_SIZE);
            
            rank2::iterator row;
            for (row = matrix->begin(); row != matrix->end(); ++row) 
              if (row->size() != nrows)
                ErrorHandling::error_function(Messages::E_INVALID_SIZE);
          }
        }
        
        private:
          unsigned int nrows, ncols, npages;
      };
    };
  
    class MaterialPropertyMaps
    {
      protected:
                            
        MaterialPropertyMap1 Sigma_f;
        MaterialPropertyMap1 nu;
        MaterialPropertyMap1 chi;
        MaterialPropertyMap1 src0;
                    
        MaterialPropertyMap1 Sigma_a;
        MaterialPropertyMap1 nuSigma_f;
        
        std::set<std::string> materials_list;
        std::map<std::string, std::string> region_material_map;
        std::map<std::string, Hermes::vector<std::string> > material_region_map;
        
        Hermes::vector<std::string> fission_materials;
        Hermes::vector<std::string> fission_regions;
        
        unsigned int G;
        
        bool1 fission_nonzero_structure;
        
        rank1 extract_rank2_diagonal(const rank2& x) const {
          rank1 result; result.reserve(G);
          for (unsigned int g = 0; g < G; g++) result.push_back(x[g][g]);
          return result;
        }
              
        void extend_to_multigroup(const MaterialPropertyMap0& mrsg_map, MaterialPropertyMap1 *mrmg_map);     
        void extend_to_multiregion(const rank1& srmg_array, MaterialPropertyMap1 *mrmg_map);
        void extend_to_multiregion_multigroup(const rank0& srsg_value, MaterialPropertyMap1 *mrmg_map);
        
        MaterialPropertyMap1 extract_map2_diagonals(const MaterialPropertyMap2& map2) const;
        
        MaterialPropertyMap1 sum_map2_columns(const MaterialPropertyMap2& map2) const;
        MaterialPropertyMap1 sum_map2_rows(const MaterialPropertyMap2& map2) const;
        
        MaterialPropertyMap2 create_map2_by_diagonals(const MaterialPropertyMap1& diags) const;
        
        void fill_with(double c, MaterialPropertyMap1 *mrmg_map);
        void fill_with(double c, MaterialPropertyMap2 *mrmg_map);
                  
        MaterialPropertyMaps(unsigned int G, const std::set<std::string>& mat_list = std::set<std::string>());          
        MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map);
                    
        virtual void validate();
        
      public:
        
        /// \brief Empty virtual destructor.
        /// Required in order to properly delete derived classes accessed through a pointer to this class.
        virtual ~MaterialPropertyMaps() {}
        
        Hermes::vector<std::string> get_regions(const std::string& material) const;
        std::string get_material(int elem_marker, const Mesh *mesh) const;
        std::string get_material(const std::string& elem_marker) const;
        
        virtual void set_iso_src(const MaterialPropertyMap1& src) {
          this->src0 = src;
        }
        
        virtual void set_iso_src(const MaterialPropertyMap0& src) {
          extend_to_multigroup(src, &this->src0);            
        }
        
        virtual void set_iso_src(const rank1& src) {
          extend_to_multiregion(src, &this->src0);
        }
        
        virtual void set_iso_src(const double& src) {
          extend_to_multiregion_multigroup(src, &this->src0);
        }
        
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
        
        const MaterialPropertyMap1& get_Sigma_a() const {
          return this->Sigma_a;
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
        const MaterialPropertyMap1& get_iso_src() const {
          return this->src0;
        }
        const bool1& get_fission_nonzero_structure() const {
          return this->fission_nonzero_structure;
        }
        const std::set<std::string>& get_materials_list() const {
          return this->materials_list;
        }
        
        virtual rank1 compute_Sigma_a(const std::string& material) const { 
          return get_Sigma_a(material); 
        }
        virtual rank2 compute_Sigma_s(const std::string& material) const = 0;
        virtual rank1 compute_Sigma_t(const std::string& material) const = 0;
        
        const rank1& get_Sigma_a(const std::string& material) const;
        const rank1& get_Sigma_f(const std::string& material) const;
        const rank1& get_nu(const std::string& material) const;
        const rank1& get_chi(const std::string& material) const;
        const rank1& get_iso_src(const std::string& material) const;
        
        void set_fission_materials(const Hermes::vector<std::string>& fission_materials);
        const Hermes::vector<std::string>& get_fission_regions() const { return fission_regions; }
        const Hermes::vector<std::string>& get_fission_materials() const { return fission_materials; }
        
        unsigned int get_G() const { return G; } 
        
        friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
    };
    
  /* MaterialProperties */
  }
  /* Common */
  }
  
  namespace Diffusion { namespace MaterialProperties
  {
    using Common::MaterialProperties::Validation;
    
    class MaterialPropertyMaps : public Common::MaterialProperties::MaterialPropertyMaps
    {
      protected:
        
        MaterialPropertyMap1 D;
        MaterialPropertyMap1 Sigma_r;
        MaterialPropertyMap2 Sigma_s;
        
        MaterialPropertyMap1 Sigma_t;
        
        bool2 scattering_nonzero_structure;
        
      public:
        
        MaterialPropertyMaps(unsigned int G, const std::set<std::string>& mat_list = std::set<std::string>()) 
          : Common::MaterialProperties::MaterialPropertyMaps(G, mat_list) { };
        MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map)
          : Common::MaterialProperties::MaterialPropertyMaps(G, reg_mat_map) { };
          
        /// \brief Empty virtual destructor.
        /// Required in order to properly delete derived classes accessed through a pointer to this class.
        virtual ~MaterialPropertyMaps() {}
        
        // We always need to supply chi, nu, Sigma_f, Sigma_r, Sigma_s and D to our neutronics weak forms. 
        // These parameters are often defined in terms of the other ones, or not specified at all and assumed 
        // to be zero for a particular simplified situation. This method, together with its complement in the
        // parent class, uses the most typical definitions to build the six-parameter set from the given input. 
        // It also checks whether the user did not enter nonsensical values. However, values entered by the 
        // user may sometimes not satisfy the common relations, as some empirical corrections may have been 
        // already included in them.
        virtual void validate();
                    
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
        
        const MaterialPropertyMap1& get_Sigma_t() const {
          return this->Sigma_t;
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
        const bool2& get_scattering_nonzero_structure() const {
          return this->scattering_nonzero_structure;
        }
        
        virtual rank1 compute_Sigma_a(const std::string& material) const;
        virtual rank2 compute_Sigma_s(const std::string& material) const { 
          return get_Sigma_s(material); 
        }
        virtual rank1 compute_Sigma_t(const std::string& material) const;
        
        const rank2& get_Sigma_s(const std::string& material) const;
        const rank1& get_Sigma_r(const std::string& material) const;
        const rank1& get_D(const std::string& material) const;
        
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
          ErrorHandling::warning("Diffusion coefficient is determined automatically according to the transport correction formula."
                  " Ignoring user setting.");
        }
        
        virtual void validate();
    };
    
  /* MaterialProperties */
  }
  /* Diffusion */
  }
  
  namespace SPN { namespace MaterialProperties
  {    
    using Common::MaterialProperties::Validation;
    
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
    
    class MaterialPropertyMaps : public Common::MaterialProperties::MaterialPropertyMaps
    {
      protected:
        
        MaterialPropertyMap3 Sigma_rn;
        MaterialPropertyMap3 odd_Sigma_rn_inv;
        
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
          : Common::MaterialProperties::MaterialPropertyMaps(G, mat_list), N(N), N_odd((N+1)/2) 
        {
          if ((N % 2) == 0) ErrorHandling::error_function(Messages::E_EVEN_SPN);
        }
        
        MaterialPropertyMaps(unsigned int G, unsigned int N, 
                             const RegionMaterialMap& reg_mat_map)
          : Common::MaterialProperties::MaterialPropertyMaps(G, reg_mat_map), N(N), N_odd((N+1)/2)  
        { 
          if ((N % 2) == 0) ErrorHandling::error_function(Messages::E_EVEN_SPN);
        }
                                
        virtual void validate();
                    
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
        
        virtual rank1 compute_Sigma_a(const std::string& material) const;
        virtual rank2 compute_Sigma_s(const std::string& material) const;
        virtual rank1 compute_Sigma_t(const std::string& material) const;
        
        const bool1  is_Sigma_rn_diagonal() const;
        
        const bool1& is_Sigma_rn_diagonal(const std::string& material) const;
        const rank3& get_Sigma_rn(const std::string& material) const;
        const rank3& get_odd_Sigma_rn_inv(const std::string& material) const;
        
        unsigned int get_N() const { return N; }
        unsigned int get_N_odd() const { return N_odd; }
        
        friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
    };
    
  /* MaterialProperties */
  }
  /* SPN */
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
      RegionMaterialMap m_map;
    public:
      region_material_map(const std::string& key, const std::string& val) {
        m_map[key] = val;
      }
      
      region_material_map& operator()(const std::string& key, const std::string& val) {
        m_map[key] = val;
        return *this;
      }
      
      operator RegionMaterialMap() {
        return m_map;
      }
  };
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}

#endif