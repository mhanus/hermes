//
//TODO: MaterialRegionMap
//

#include "neutronics/material_properties.h"
#include <iomanip>

namespace Hermes { namespace Hermes2D { namespace Neutronics
{
  namespace Common { namespace MaterialProperties
  {
    MaterialPropertyMaps::MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map)
      : region_material_map(reg_mat_map), G(G)
    {
      RegionMaterialMap::const_iterator it = reg_mat_map.begin();
      for ( ; it != reg_mat_map.end(); ++it)
      {
        materials_list.insert(it->second);
        material_region_map[it->second].push_back(it->first);
      }
    }
    
    MaterialPropertyMaps::MaterialPropertyMaps(unsigned int G, const std::set<std::string>& mat_list) 
      : materials_list(mat_list), G(G)  
    { 
      std::set<std::string>::const_iterator it = mat_list.begin();
      for ( ; it != mat_list.end(); ++it)
        material_region_map[*it].push_back(*it);
    };

    void MaterialPropertyMaps::extend_to_multigroup(const MaterialPropertyMap0& mrsg_map, 
                                                    MaterialPropertyMap1 *mrmg_map)
    {          
      MaterialPropertyMap0::const_iterator it;
      for (it = mrsg_map.begin(); it != mrsg_map.end(); ++it)
        (*mrmg_map)[it->first].assign(G, it->second);
      
    }
    
    void MaterialPropertyMaps::extend_to_multiregion(const rank1& srmg_array, 
                                                      MaterialPropertyMap1 *mrmg_map)
    {
      if (materials_list.empty())
        error_function(Messages::E_MR_EXTENSION);
      
      std::set<std::string>::const_iterator it;
      for (it = materials_list.begin(); it != materials_list.end(); ++it)
        (*mrmg_map)[*it] = srmg_array;
    }
    
    void MaterialPropertyMaps::extend_to_multiregion_multigroup(const rank0& srsg_value, 
                                                                MaterialPropertyMap1 *mrmg_map)
    {
      if (materials_list.empty())
        error_function(Messages::E_MR_EXTENSION);
      
      std::set<std::string>::const_iterator it;
      for (it = materials_list.begin(); it != materials_list.end(); ++it)
        (*mrmg_map)[*it].assign(G, srsg_value);
    }
    
    void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap1 *mrmg_map)
    {
      if (materials_list.empty())
        error_function(Messages::E_MR_EXTENSION);
      
      std::set<std::string>::const_iterator it;
      for (it = materials_list.begin(); it != materials_list.end(); ++it)
        (*mrmg_map)[*it].assign(G, c);
    }
    
    void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap2 *mrmg_map)
    {
      std::set<std::string>::const_iterator it;
      for (it = materials_list.begin(); it != materials_list.end(); ++it)
        (*mrmg_map)[*it].assign(G, rank1(G, c));
    }
    
    MaterialPropertyMap1 MaterialPropertyMaps::extract_map2_diagonals(const MaterialPropertyMap2& map2) const
    {
      MaterialPropertyMap1 diags;
      
      MaterialPropertyMap2::const_iterator map2_it = map2.begin();
      for ( ; map2_it != map2.end(); ++map2_it)
        diags[map2_it->first] = extract_rank2_diagonal(map2_it->second);
      
      return diags;
    }
    
    MaterialPropertyMap1 MaterialPropertyMaps::sum_map2_columns(const MaterialPropertyMap2& map2) const
    {
      MaterialPropertyMap1 summed;
      
      MaterialPropertyMap2::const_iterator map2_it = map2.begin();
      for ( ; map2_it != map2.end(); ++map2_it)
      {
        summed[map2_it->first].reserve(G);
        for (unsigned int gfrom = 0; gfrom < G; gfrom++)
        {
          double sum = 0.0;
          
          for (unsigned int gto = 0; gto < G; gto++)
            sum += map2_it->second[gto][gfrom];
          
          summed[map2_it->first].push_back(sum);    
        }
      }
      
      return summed;
    }
    
    MaterialPropertyMap1 MaterialPropertyMaps::sum_map2_rows(const MaterialPropertyMap2& map2) const
    {
      MaterialPropertyMap1 summed;
      
      MaterialPropertyMap2::const_iterator map2_it = map2.begin();
      for ( ; map2_it != map2.end(); ++map2_it)
      {
        summed[map2_it->first].reserve(G);
        for (unsigned int gto = 0; gto < G; gto++)
        {
          double sum = 0.0;
          
          for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            sum += map2_it->second[gto][gfrom];
          
          summed[map2_it->first].push_back(sum);    
        }
      }
      
      return summed;
    }
    
    MaterialPropertyMap2 MaterialPropertyMaps::create_map2_by_diagonals(const MaterialPropertyMap1& diags) const 
    {
      MaterialPropertyMap2 map2;
      
      MaterialPropertyMap1::const_iterator diags_it = diags.begin();
      for ( ; diags_it != diags.end(); ++diags_it)
      {
        map2[diags_it->first].resize(G, rank1(G, 0.0));
        
        for (unsigned int g = 0; g < G; g++)
          map2[diags_it->first][g][g] = diags_it->second[g];
      }
      
      return map2;
    }
            
    void MaterialPropertyMaps::validate()
    {             
      fission_nonzero_structure = bool1(G, false);
      
      if (chi.empty())
      {
        fill_with(0.0, &chi);
        MaterialPropertyMap1::iterator it = chi.begin();
        for ( ; it != chi.end(); ++it)
          it->second[0] = 1.0;
        fission_nonzero_structure[0] = true;
      }
      else
      {
        for (unsigned int g = 0; g < G; g++)
        {
          MaterialPropertyMap1::const_iterator it = chi.begin();
          for ( ; it != chi.end(); ++it)
          {
            if (fabs(it->second[g]) > 1e-14)
            {
              fission_nonzero_structure[g] = true;
              break;
            }
          }
        }
      }
      
      if (nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
        nu = NDArrayMapOp::divide<rank1>(nuSigma_f, Sigma_f);
      else if (nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
        nuSigma_f = NDArrayMapOp::multiply<rank1>(nu, Sigma_f);
      else if (Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
        Sigma_f = NDArrayMapOp::divide<rank1>(nuSigma_f, nu);
      else if (!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
      {
        MaterialPropertyMap1 diff = NDArrayMapOp::subtract<rank1>(nuSigma_f, 
                                                                  NDArrayMapOp::multiply<rank1>(nu, Sigma_f) );
        std::for_each(diff.begin(), diff.end(), Validation::ensure_trivial());
      }
      else
      {
        if (src0.empty())
          warning(Messages::W_NO_FISSION);
        fill_with(0.0, &nu);
        fill_with(0.0, &chi);
        fill_with(0.0, &Sigma_f);
        fission_nonzero_structure = bool1(G, false);
      }
      
      if ((nu.size() != Sigma_f.size()) || (nu.size() != chi.size()))
        error_function(Messages::E_NONMATCHING_PROPERTIES);
      
      if (!Sigma_f.empty())
      {
        std::for_each(nu.begin(), nu.end(), Validation::ensure_size(G));
        std::for_each(Sigma_f.begin(), Sigma_f.end(), Validation::ensure_size(G));
        std::for_each(chi.begin(), chi.end(), Validation::ensure_size(G));
      }
      
      if (!src0.empty())
        std::for_each(src0.begin(), src0.end(), Validation::ensure_size(G));
      
      if (!Sigma_a.empty())
      {
        // Warn if \Sigma_a < \Sigma_f for any region (this indicates an unphysical situation, since
        // by definition \Sigma_a = \Sigma_f + \Sigma_c + \Sigma_{n,p} + other possible reactions
        // leading to neutron removal).
        MaterialPropertyMap1::const_iterator ita = Sigma_a.begin();
        MaterialPropertyMap1::const_iterator itf = Sigma_f.begin();
        for ( ; ita != Sigma_a.end(); ++ita, ++itf)
        {
          rank1::const_iterator a = ita->second.begin();
          rank1::const_iterator f = itf->second.begin();
          
          for ( ; a != ita->second.end(); ++a,++f)
            if (*a < *f)
              warning(Messages::W_SA_LT_SF);
        }
      }
    }
    
    Hermes::vector<std::string> MaterialPropertyMaps::get_regions(const std::string& material) const
    {
      std::map<std::string, Hermes::vector<std::string> >::const_iterator it;
      for (it = material_region_map.begin(); it != material_region_map.end(); ++it)
        if (it->first == material)
          return it->second;
    }
        
    std::string MaterialPropertyMaps::get_material(int elem_marker, Mesh *mesh) const 
    { 
      std::string region;
      
      if (elem_marker == HERMES_DUMMY_ELEM_MARKER)
        region = this->nu.begin()->first; 
      else
        region = mesh->get_element_markers_conversion().get_user_marker(elem_marker).marker;
      
      return get_material(region);
    }
    
    std::string MaterialPropertyMaps::get_material(const std::string& elem_marker) const 
    { 
      RegionMaterialMap::const_iterator material = this->region_material_map.find(elem_marker);
      
      if (material != this->region_material_map.end())
        return material->second;
      else
        return elem_marker; // Corresponds to the case when elem_marker <==> material, 
    }
    
    const rank1& MaterialPropertyMaps::get_Sigma_a(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->Sigma_a.find(material);
      if (data != this->Sigma_a.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_Sigma_f(const std::string& material) const
    {
      // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
      // which operator[] is undefined.
      MaterialPropertyMap1::const_iterator data = this->Sigma_f.find(material);
      if (data != this->Sigma_f.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_nu(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->nu.find(material);
      if (data != this->nu.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_chi(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->chi.find(material);
      if (data != this->chi.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_iso_src(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->src0.find(material);
      if (data != this->src0.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
    {
      using namespace std;
      
      os << endl;
      os << setw(12) << "target group" << setw(10) << "chi" << setw(10) << "nu";
      os << setw(10) << "Sigma_f" << setw(14) << "iso. ext. src" << endl; 
      
      MaterialPropertyMap1::const_iterator data_elem = matprop.chi.begin();
      for ( ; data_elem != matprop.chi.end(); ++data_elem)
      {
        std::string mat = data_elem->first;
        
        os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
        os << setw(40) << mat << endl;
        os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
        for (unsigned int gto = 0; gto < matprop.G; gto++)
        {
          os << setw(6) << gto << setw(6) << ' ';
          os << setw(10) << matprop.get_chi(mat)[gto];
          os << setw(10) << matprop.get_nu(mat)[gto];
          os << setw(10) << matprop.get_Sigma_f(mat)[gto];
          os << setw(14);
          if (matprop.src0.empty())
            os << "N/A";
          else
            os << matprop.get_iso_src(mat)[gto];
          
          os << endl;
        }
      }
      
      os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
      os << "All-region fission spectrum: ";
      
      for (unsigned int g = 0; g < matprop.G; g++)
        os << matprop.get_fission_nonzero_structure()[g] << ' ';
      
      return os << endl;
    }
    
  /* MaterialProperties */
  }
  /* Common */
  }
  
  namespace Diffusion { namespace MaterialProperties
  {        
    void MaterialPropertyMaps::validate()
    {
      Common::MaterialProperties::MaterialPropertyMaps::validate();
      
      bool D_given = !D.empty();
      bool Sigma_r_given = !Sigma_r.empty();
      bool Sigma_s_given = !Sigma_s.empty();
      bool Sigma_t_given = !Sigma_t.empty();
      bool Sigma_a_given = !Sigma_a.empty();
      bool Sigma_f_given = !Sigma_f.empty();
      
      if (!Sigma_r_given)
      {
        // If Sigma_r is not given, we can calculate it from Sigma_t and Sigma_s.
        
        if (Sigma_t_given)
        {
          if (!Sigma_s_given)
          {
            if (Sigma_a_given)
            {
              // If Sigma_s is not given, but Sigma_a is, we can calculate Sigma_s from Sigma_t and Sigma_a.
              Sigma_s = create_map2_by_diagonals(NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_a));
            }
            else 
            {
              // If only Sigma_t is given, we assume that all reaction terms are included in Sigma_t; all
              // other x-sections will be set to zero.
              warning(Messages::W_NO_SCATTERING);
              fill_with(0.0, &Sigma_s);
            }
            
            Sigma_s_given = true;
          }
        }
        else
        {
          // If Sigma_t is not given, but Sigma_a and Sigma_s are, we can obtain Sigma_t from the latter two.
          
          if (!Sigma_s_given)
          {
            warning(Messages::W_NO_SCATTERING);
            fill_with(0.0, &Sigma_s);
            Sigma_s_given = true;
          }
          
          if (Sigma_a_given)
            Sigma_t = NDArrayMapOp::add<rank1>(Sigma_a, sum_map2_columns(Sigma_s));
          else 
          {
            // If neither Sigma_r, Sigma_t, Sigma_a are given, we may have a purely fissioning system.
            if (Sigma_f_given)
              Sigma_t = Sigma_f;
            else
              error_function(Messages::E_INSUFFICIENT_DATA);
          }
          
          Sigma_t_given = true;
        }
        
        Sigma_r = NDArrayMapOp::subtract<rank1>(Sigma_t, extract_map2_diagonals(Sigma_s));
        Sigma_r_given = true;
      }
      
      // Now, we surely have Sigma_r ...
      
      scattering_nonzero_structure = bool2(G, std::vector<bool>(G, false));
      
      if (!Sigma_s_given)
      {
        // If Sigma_s is not given, but Sigma_t is, we can obtain the former from the latter and from Sigma_r.
        // Note that the execution will come here only if the user entered Sigma_r himself - otherwise, Sigma_s
        // has been already set in the previous test case.
        
        if (Sigma_t_given)
        {
          Sigma_s = create_map2_by_diagonals(NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_r));
        }
        else
        {
          warning(Messages::W_NO_SCATTERING);
          fill_with(0.0, &Sigma_s);
        }
        
        Sigma_s_given = true;
      }
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        for (unsigned int gfrom = 0; gfrom < G; gfrom++)
        {
          MaterialPropertyMap2::const_iterator Ss_it = Sigma_s.begin();
          for ( ; Ss_it != Sigma_s.end(); ++Ss_it)
          {
            if (fabs(Ss_it->second[gto][gfrom]) > 1e-14)
            {
              scattering_nonzero_structure[gto][gfrom] = true;
              break;
            }
          }
        }
      }
      
      // Now, we surely have Sigma_s and Sigma_r, one parameter to go ...
      
      if (!D_given)
      {
        MaterialPropertyMap1::const_iterator Sr_elem = Sigma_r.begin();
        for ( ; Sr_elem != Sigma_r.end(); ++Sr_elem)
        {
          D[Sr_elem->first].resize(G);
          for (unsigned int g = 0; g < G; g++)
          {
            if (Sigma_t_given)
              D[Sr_elem->first][g] = 1./(3.*Sigma_t[Sr_elem->first][g]);
            else
              D[Sr_elem->first][g] = 1./(3.*Sr_elem->second[g]);
          }
        }
          
        D_given = true;
      }
      
      if ((D.size() != Sigma_r.size()) || (D.size() != Sigma_s.size()))
        error_function(Messages::E_NONMATCHING_PROPERTIES);
      
      std::for_each(Sigma_s.begin(), Sigma_s.end(), Validation::ensure_size(G,G));
      std::for_each(Sigma_r.begin(), Sigma_r.end(), Validation::ensure_size(G));
      std::for_each(D.begin(), D.end(), Validation::ensure_size(G));
    }
    
    rank1 MaterialPropertyMaps::compute_Sigma_a(const std::string& material) const
    {
      if (!Sigma_a.empty())
        return Common::MaterialProperties::MaterialPropertyMaps::get_Sigma_a(material);
      
      MaterialPropertyMap1::const_iterator Sr_mat = this->Sigma_r.find(material);
      MaterialPropertyMap2::const_iterator Ss_mat = this->Sigma_s.find(material);
      if (Sr_mat != this->Sigma_r.end() && Ss_mat != this->Sigma_s.end())
      {
        rank1 Sr = Sr_mat->second;
        rank2 Ss = Ss_mat->second;
        rank1 Sa = Sr;
        
        for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          for (unsigned int gto = 0; gto < G; gto++)
            if (gfrom != gto)
              Sa[gfrom] -= Ss[gto][gfrom];
            
        return Sa;
      }
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    rank1 MaterialPropertyMaps::compute_Sigma_t(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator Sr_mat = this->Sigma_r.find(material);
      MaterialPropertyMap2::const_iterator Ss_mat = this->Sigma_s.find(material);
      if (Sr_mat != this->Sigma_r.end() && Ss_mat != this->Sigma_s.end())
      {
        rank1 Sr = Sr_mat->second;
        rank1 Ssd = extract_rank2_diagonal(Ss_mat->second);
        return NDArrayMapOp::add<rank1>(Sr, Ssd);
      }
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    const rank2& MaterialPropertyMaps::get_Sigma_s(const std::string& material) const
    {
      // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
      // which operator[] is undefined.
      MaterialPropertyMap2::const_iterator data = this->Sigma_s.find(material);
      if (data != this->Sigma_s.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank2()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_Sigma_r(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->Sigma_r.find(material);
      if (data != this->Sigma_r.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    const rank1& MaterialPropertyMaps::get_D(const std::string& material) const
    {
      MaterialPropertyMap1::const_iterator data = this->D.find(material);
      if (data != this->D.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
    {
      using namespace std;
      
      os << static_cast<const Common::MaterialProperties::MaterialPropertyMaps&>(matprop) << endl;
      
      os << setw(12) << "target group" << setw(10) << "D" << setw(10) << "Sigma_r";
      os << setw(22) << "Sigma_s" << endl; 
      
      MaterialPropertyMap1::const_iterator data_elem = matprop.Sigma_r.begin();
      for ( ; data_elem != matprop.Sigma_r.end(); ++data_elem)
      {
        std::string mat = data_elem->first;
        
        os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
        os << setw(40) << mat << endl;
        os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
        for (unsigned int gto = 0; gto < matprop.G; gto++)
        {
          os << setw(6) << gto << setw(6) << ' ';
          os << setw(10) << matprop.get_D(mat)[gto];
          os << setw(10) << matprop.get_Sigma_r(mat)[gto];
          
          for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
            os << setw(8) << matprop.get_Sigma_s(mat)[gto][gfrom];
          
          os << endl;
        }
      }
      
      os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
      os << "All-region scattering spectrum: " << endl;
      for (unsigned int gto = 0; gto < matprop.G; gto++)
      {
        for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
          os << setw(10) << matprop.get_scattering_nonzero_structure()[gto][gfrom];
        os << endl;
      }
      
      return os << endl;
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G, 
                                                                                    const MaterialPropertyMap2& Ss_1,
                                                                                    const RegionMaterialMap& reg_mat_map)
      : MaterialPropertyMaps(G, reg_mat_map)
    {
      std::set<std::string> matlist;
      MaterialPropertyMap2::const_iterator it = Ss_1.begin();
      for( ; it != Ss_1.end(); ++it)
        matlist.insert(it->first);
      
      set_materials_list(matlist);
      
      this->Sigma_s_1_out = sum_map2_columns(Ss_1);
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G, 
                                                                                    const MaterialPropertyMap1& mu_av,
                                                                                    const RegionMaterialMap& reg_mat_map)
      : MaterialPropertyMaps(G, reg_mat_map)
    {
      std::set<std::string> matlist;
      MaterialPropertyMap1::const_iterator it = mu_av.begin();
      for( ; it != mu_av.end(); ++it)
        matlist.insert(it->first);
      
      set_materials_list(matlist);
      
      this->mu_av = mu_av;
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G, 
                                                                                    const MaterialPropertyMap0& mu_av,
                                                                                    const RegionMaterialMap& reg_mat_map)
      : MaterialPropertyMaps(G, reg_mat_map)
    {
      std::set<std::string> matlist;
      MaterialPropertyMap0::const_iterator it = mu_av.begin();
      for( ; it != mu_av.end(); ++it)
        matlist.insert(it->first);
      
      set_materials_list(matlist);
      
      extend_to_multigroup(mu_av, &this->mu_av);
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G,
                                                                                    const RegionMaterialMap& reg_mat_map,
                                                                                    const rank1& mu_av )
      : MaterialPropertyMaps(G, reg_mat_map)
    {
      extend_to_multiregion(mu_av, &this->mu_av);
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G,
                                                                                    const RegionMaterialMap& reg_mat_map,
                                                                                    const rank0& mu_av )
      : MaterialPropertyMaps(G, reg_mat_map)
    {
      extend_to_multiregion_multigroup(mu_av, &this->mu_av);
    }

    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G,
                                                                                    std::set<std::string> mat_list,
                                                                                    const rank1& mu_av )
      : MaterialPropertyMaps(G, mat_list)
    {
      extend_to_multiregion(mu_av, &this->mu_av);
    }
    
    TransportCorrectedMaterialPropertyMaps::TransportCorrectedMaterialPropertyMaps( unsigned int G,
                                                                                    std::set<std::string> mat_list,
                                                                                    const rank0& mu_av )
      : MaterialPropertyMaps(G, mat_list)
    {
      extend_to_multiregion_multigroup(mu_av, &this->mu_av);
    }
    
    void TransportCorrectedMaterialPropertyMaps::validate()
    {          
      bool mu_av_given = !mu_av.empty();
      bool Sigma_s_1_out_given = !Sigma_s_1_out.empty();
      
      if (mu_av_given)
        D = mu_av;
      else if (Sigma_s_1_out_given)
        D = Sigma_s_1_out;
      else
        error_function(Messages::E_INSUFFICIENT_DATA);
      
      MaterialPropertyMaps::validate();
      
      if (Sigma_t.empty())
        error_function(Messages::E_INSUFFICIENT_DATA);
      
      if (!Sigma_s_1_out_given)
        Sigma_s_1_out = NDArrayMapOp::multiply<rank1>(mu_av, sum_map2_columns(Sigma_s));
        
      MaterialPropertyMap1 Sigma_tr = NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_s_1_out);
      MaterialPropertyMap1::const_iterator Str_elem = Sigma_tr.begin();
      for ( ; Str_elem != Sigma_tr.end(); ++Str_elem)
        for (unsigned int g = 0; g < G; g++)
          D[Str_elem->first][g] = 1./(3.*Str_elem->second[g]);
    }
    
  /* MaterialProperties */
  }
  /* Diffusion */
  }
  
  namespace SPN { namespace MaterialProperties
  {
    void MaterialPropertyMaps::extend_to_rank3(const MaterialPropertyMap2& src, MaterialPropertyMap3* dest)
    {
      for (MaterialPropertyMap2::const_iterator src_it = src.begin(); src_it != src.end(); ++src_it)
        for (unsigned int n = 0; n <= N; n++)
          (*dest)[src_it->first].push_back(src_it->second);
    }
    
    void MaterialPropertyMaps::extend_to_rank3(const MaterialPropertyMap1& src, MaterialPropertyMap3* dest)
    {
      MaterialPropertyMap2 matrices = create_map2_by_diagonals(src);
      extend_to_rank3(matrices, dest);
    }
    
    MaterialPropertyMap3 MaterialPropertyMaps::create_map3_by_diagonals(const MaterialPropertyMap2& diags) const 
    {
      MaterialPropertyMap3 map3;
      
      MaterialPropertyMap2::const_iterator diags_it = diags.begin();
      for ( ; diags_it != diags.end(); ++diags_it)
      {
        map3[diags_it->first].resize(N+1, rank2(G, rank1(G, 0.0)));
        
        for (unsigned int n = 0; n <= N; n++)
          for (unsigned int g = 0; g < G; g++)
            map3[diags_it->first][n][g][g] = diags_it->second[n][g];
      }
      
      return map3;
    }
    
    void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap3 *mmmrmg_map)
    {
      std::set<std::string>::const_iterator it;
      for (it = materials_list.begin(); it != materials_list.end(); ++it)
        (*mmmrmg_map)[*it].assign(N+1, rank2(G, rank1(G, c)));
    }
            
    void MaterialPropertyMaps::invert_odd_Sigma_rn()
    {          
      MaterialPropertyMap3::const_iterator Sr_mat = Sigma_rn.begin();
      for ( ; Sr_mat != Sigma_rn.end(); ++Sr_mat)
      {
        std::string mat = Sr_mat->first;
        odd_Sigma_rn_inv[mat].resize(N_odd, rank2(G, rank1(G, 0.0)));
        
        rank3 moment_matrices = Sr_mat->second;
        rank3::const_iterator odd_moment_matrix = moment_matrices.begin()+1;
        rank3::iterator inverted_moment_matrix = odd_Sigma_rn_inv[mat].begin();
        bool1::const_iterator moment_matrix_is_diagonal = Sigma_rn_is_diagonal[mat].begin()+1;
        for (unsigned int num_visited = 0; num_visited < N_odd; num_visited++, 
              odd_moment_matrix += 2, ++inverted_moment_matrix, moment_matrix_is_diagonal += 2)
        {              
          if (*moment_matrix_is_diagonal)
          {
            for (unsigned int g = 0; g < G; g++)
              (*inverted_moment_matrix)[g][g] = 1./(*odd_moment_matrix)[g][g];
          }
          else
          {               
            int lda = G;
            int n = G;
            int sz_dgetri_workspace = G*G;
            int info;
            
            double *A = new double [n*n];
            int *ipiv = new int [n+1];
            double *dgetri_workspace = new double [sz_dgetri_workspace];
            double *dgecon_double_workspace = new double [4*n];
            int *dgecon_int_workspace = new int [n];
            
            rank2::const_iterator row = odd_moment_matrix->begin();  int gto = 0;
            for ( ; row != odd_moment_matrix->end(); ++row, ++gto)
              std::copy(row->begin(), row->end(), &A[gto*G]);
          
            double anorm = dlange_("1", &n, &n, A, &lda, NULL); // array in the last argument is not referenced for "1"-norm
            
            dgetrf_(&n, &n, A, &lda, ipiv, &info);
            if (info != 0) 
              error_function(Messages::E_LAPACK_ERROR, info);
            
            double rcond;
            dgecon_("1", &n, A, &lda, &anorm, &rcond, 
                    dgecon_double_workspace, dgecon_int_workspace, &info);
            if (info != 0) 
              error_function(Messages::E_LAPACK_ERROR, info);
            if (rcond > 1e12)
              warning(Messages::W_SINGULAR_MATRIX);
                    
            dgetri_(&n, A, &lda, ipiv, dgetri_workspace, &sz_dgetri_workspace, &info);
            if (info != 0)
              error_function(Messages::E_LAPACK_ERROR);
                            
            rank2::iterator inv_mtx_row = inverted_moment_matrix->begin(); gto = 0;
            for ( ; inv_mtx_row != inverted_moment_matrix->end(); ++inv_mtx_row, ++gto)
              std::copy(&A[gto*G], &A[(gto+1)*G], inv_mtx_row->begin());
            
            delete [] ipiv;
            delete [] dgetri_workspace;
            delete [] dgecon_double_workspace;
            delete [] dgecon_int_workspace;
            delete [] A;
          }
        }
      }
    }
    
    void MaterialPropertyMaps::validate()
    {
      if (Sigma_tn.empty())
        error_function(Messages::E_SIGMA_T_REQUIRED);
      
      Common::MaterialProperties::MaterialPropertyMaps::validate();
                
      if (Sigma_sn.empty())
      {
        warning(Messages::W_NO_SCATTERING);
        fill_with(0.0, &Sigma_sn);
      }
      
      std::for_each(Sigma_tn.begin(), Sigma_tn.end(), Validation::ensure_size(G,G,N+1));
      
      MaterialPropertyMap3::const_iterator Stn_material = Sigma_tn.begin();
      for ( ; Stn_material != Sigma_tn.end(); ++Stn_material)
      {
        std::string mat = Stn_material->first;
        rank3 Ssn = Sigma_sn[mat];
        rank3 Stn = Sigma_tn[mat];
          
        Sigma_rn_is_diagonal[mat].resize(N+1, true);
        
        if (Ssn.size() > N+1)
        {
          warning(Messages::W_SCATTERING_TRUNCATION, N+1);
          Sigma_sn[mat].erase(Sigma_sn[mat].begin() + N+1, Sigma_sn[mat].end());
          Ssn = Sigma_sn[mat];
        }
        
        if (Ssn.size() == N+1)
        {
          Sigma_rn[mat] = NDArrayMapOp::subtract<rank3>(Stn, Ssn);
          
          rank3 moment_matrices = Sigma_rn[mat];
          rank3::const_iterator moment_matrix = moment_matrices.begin();
          bool1::iterator moment_matrix_is_diagonal = Sigma_rn_is_diagonal[mat].begin();
          for ( ; moment_matrix != moment_matrices.end(); ++moment_matrix, ++moment_matrix_is_diagonal)
          { 
            if (moment_matrix->size() != G)
              error_function(Messages::E_INVALID_SIZE);
            
            for (unsigned int gto = 0; gto < G && *moment_matrix_is_diagonal; gto++)
            {
              if (moment_matrix->at(gto).size() != G)
                error_function(Messages::E_INVALID_SIZE);
              
              for (unsigned int gfrom = 0; gfrom < G && *moment_matrix_is_diagonal; gfrom++)
                if (gfrom != gto && fabs((*moment_matrix)[gto][gfrom]) > 1e-12)
                  *moment_matrix_is_diagonal = false;
            }
          }
        }
        else if (Ssn.size() < N+1)
        {
          Sigma_rn[mat] = Stn;
          
          // Totally absorbing material in case of Ssn.size() == 0.
          for (unsigned int n = 0; n < Ssn.size(); n++)
          {
            Sigma_rn[mat][n] = NDArrayMapOp::subtract<rank2>(Stn[n], Ssn[n]);
            
            bool isdiag = true;
            for (unsigned int gto = 0; gto < G && isdiag; gto++)
            {
              if (Ssn[n].size() != G)
                error_function(Messages::E_INVALID_SIZE);
              
              for (unsigned int gfrom = 0; gfrom < G && isdiag; gfrom++)
              {
                if (Ssn[n][gto].size() != G)
                  error_function(Messages::E_INVALID_SIZE);
                
                if (gfrom != gto && fabs(Ssn[n][gto][gfrom]) > 1e-12)
                  isdiag = false;
              }
            }
            
            Sigma_rn_is_diagonal[mat][n] = isdiag;
          }
        }
      }
      
      invert_odd_Sigma_rn();
    }
    
    rank1 MaterialPropertyMaps::compute_Sigma_a(const std::string& material) const
    {
      if (!Sigma_a.empty())
        return Common::MaterialProperties::MaterialPropertyMaps::get_Sigma_a(material);
      
      MaterialPropertyMap3::const_iterator Sr_mat = this->Sigma_rn.find(material);
      MaterialPropertyMap3::const_iterator Ss_mat = this->Sigma_sn.find(material);
      if (Sr_mat != this->Sigma_rn.end() && Ss_mat != this->Sigma_sn.end())
      {
        rank1 Sr = extract_rank2_diagonal(Sr_mat->second.front());
        rank2 Ss = Ss_mat->second.front();
        rank1 Sa = Sr;
        
        for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          for (unsigned int gto = 0; gto < G; gto++)
            if (gfrom != gto)
              Sa[gfrom] -= Ss[gto][gfrom];
            
          return Sa;
      }
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    rank1 MaterialPropertyMaps::compute_Sigma_t(const std::string& material) const
    {
      MaterialPropertyMap3::const_iterator Stn_mat = this->Sigma_tn.find(material);
      if (Stn_mat != this->Sigma_tn.end())
        return extract_rank2_diagonal(Stn_mat->second.front());
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    rank2 MaterialPropertyMaps::compute_Sigma_s(const std::string& material) const
    {
      MaterialPropertyMap3::const_iterator Ss_mat = this->Sigma_sn.find(material);
      if (Ss_mat != this->Sigma_sn.end())
        return Ss_mat->second.front();
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank2()); // To avoid MSVC problems; execution should never come to this point.
      }
    }

    const rank3& MaterialPropertyMaps::get_Sigma_rn(const std::string& material) const
    {
      MaterialPropertyMap3::const_iterator data = this->Sigma_rn.find(material);
      if (data != this->Sigma_rn.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank3()); // To avoid MSVC problems; execution should never come to this point.
      }
    }
    
    const rank3& MaterialPropertyMaps::get_odd_Sigma_rn_inv(const std::string& material) const
    {
      MaterialPropertyMap3::const_iterator data = this->odd_Sigma_rn_inv.find(material);
      if (data != this->odd_Sigma_rn_inv.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new rank3()); // To avoid MSVC problems; execution should never come to this point.
      }
    }

    const bool1& MaterialPropertyMaps::is_Sigma_rn_diagonal(const std::string& material) const
    {
      std::map<std::string, bool1>::const_iterator data = this->Sigma_rn_is_diagonal.find(material);
      if (data != this->Sigma_rn_is_diagonal.end())
        return data->second;
      else
      {
        error_function(Messages::E_INVALID_MARKER);
        return *(new bool1()); // To avoid MSVC problems; execution should never come to this point.
      }
    }

    const bool1 MaterialPropertyMaps::is_Sigma_rn_diagonal() const
    {
      bool1 ret(N+1, true);
      
      for (unsigned int k = 0; k <= N; k++)
      {
        std::map<std::string, bool1>::const_iterator data_it = Sigma_rn_is_diagonal.begin();
        for ( ; data_it != Sigma_rn_is_diagonal.end(); ++data_it)
        {
          if (data_it->second[k] == false)
          {
            ret[k] = false;
            break;
          }
        }
      }
        
      return ret;    
    }
    
    std::ostream& operator<<(std::ostream& os, const MaterialPropertyMaps& matprop)
    {
      using namespace std;
      
      os << static_cast<const Common::MaterialProperties::MaterialPropertyMaps&>(matprop) << endl;
      
      int gto_width = 12;
      int elem_width = 14;
      int total_width = 2*gto_width + 2*elem_width*matprop.G;          
                          
      MaterialPropertyMap3::const_iterator Srn_elem = matprop.Sigma_rn.begin();
      MaterialPropertyMap3::const_iterator Srn_inv_elem = matprop.odd_Sigma_rn_inv.begin();
      for ( ; Srn_elem != matprop.Sigma_rn.end(); ++Srn_elem, ++Srn_inv_elem)
      {
        std::string mat = Srn_elem->first;
        rank3 Srn_moments = Srn_elem->second;
        rank3 Srn_inv_moments = Srn_inv_elem->second;
        
        os << setw(total_width) << setfill('_') << ' ' << endl << setfill(' ');
        os << setw(total_width/2+mat.length()/2) << mat << endl << endl;
               
        os << setw(gto_width)  << "trgt group";
        os << setw(matprop.G*elem_width) << "Sigma_rn";
        os << setw(matprop.G*elem_width) << "odd_Srn_inv";
        
        rank3::const_iterator Srn_moment = Srn_moments.begin();
        rank3::const_iterator Srn_inv_moment = Srn_inv_moments.begin();
        int moment = 0;
        for ( ; Srn_moment != Srn_moments.end(); ++Srn_moment, ++moment)
        {
          os << endl << setw(total_width/2+2) << "moment " << moment << endl;

          for (unsigned int gto = 0; gto < matprop.G; gto++)
          {
            os << setw(gto_width) << gto;
            for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
              os << setw(elem_width) << (*Srn_moment)[gto][gfrom];
            
            if (moment % 2)
              for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
                os << setw(elem_width) << (*Srn_inv_moment)[gto][gfrom];
                          
            os << endl;
          }
          
          if (moment % 2)
            ++Srn_inv_moment;
        }
      }
      
      return os << endl;
    }
    
  /* MaterialProperties */
  }
  /* SPN */
  }    
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}