#include "../hermes2d.h"

#include <algorithm>
#include <iomanip>

namespace WeakFormsNeutronics
{
  namespace Monoenergetic
  {    
    namespace Diffusion 
    {
      DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( Hermes::vector<std::string> regions, 
                                                              Hermes::vector<double> D_map, 
                                                              Hermes::vector<double> Sigma_a_map, 
                                                              Hermes::vector<double> Q_map ) : WeakForm(1) 
      {
        using namespace WeakFormsH1;
        
        for (unsigned int i = 0; i < regions.size(); i++)
        {
          /* Jacobian */
          // Diffusion.
          add_matrix_form(new DefaultJacobianDiffusion(0, 0, regions[i], new HermesFunction(D_map[i]), 
                                                       HERMES_SYM));
          // Absorption.
          add_matrix_form(new DefaultMatrixFormVol(0, 0, regions[i], new HermesFunction(Sigma_a_map[i]), 
                                                   HERMES_SYM));
          
          /* Residual */
          // Diffusion.
          add_vector_form(new DefaultResidualDiffusion(0, regions[i], new HermesFunction(D_map[i])));
          // Absorption.
          add_vector_form(new DefaultResidualVol(0, regions[i], new HermesFunction(Sigma_a_map[i])));
          // Sources.
          add_vector_form(new DefaultVectorFormVol(0, regions[i], new HermesFunction(-Q_map[i])));
        }
      }
    }
  }
      
  namespace Multigroup
  { 
    namespace MaterialProperties
    {
      namespace Common
      {
        MaterialPropertyMaps::MaterialPropertyMaps(unsigned int G, const RegionMaterialMap& reg_mat_map)
          : region_material_map(reg_mat_map), G(G)
        {
          RegionMaterialMap::const_iterator it = reg_mat_map.begin();
          for ( ; it != reg_mat_map.end(); ++it)
            materials_list.insert(it->second);
        }

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
            error(E_MR_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            (*mrmg_map)[*it] = srmg_array;
        }
        
        void MaterialPropertyMaps::extend_to_multiregion_multigroup(const rank0& srsg_value, 
                                                                    MaterialPropertyMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            (*mrmg_map)[*it].assign(G, srsg_value);
        }
        
        void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          
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
          using namespace ValidationFunctors;
          
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
            std::for_each(diff.begin(), diff.end(), ensure_trivial());
          }
          else
          {
            if (src0.empty())
              warning(W_NO_FISSION);
            fill_with(0.0, &nu);
            fill_with(0.0, &chi);
            fill_with(0.0, &Sigma_f);
            fission_nonzero_structure = bool1(G, false);
          }
          
          if ((nu.size() != Sigma_f.size()) || (nu.size() != chi.size()))
            error(E_NONMATCHING_PROPERTIES);
          
          if (!Sigma_f.empty())
          {
            std::for_each(nu.begin(), nu.end(), ensure_size(G));
            std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size(G));
            std::for_each(chi.begin(), chi.end(), ensure_size(G));
          }
          
          if (!src0.empty())
            std::for_each(src0.begin(), src0.end(), ensure_size(G));
          
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
                  warning(W_SA_LT_SF);
            }
          }
        }
        
        std::string MaterialPropertyMaps::get_material(int elem_marker, WeakForm *wf) const 
        { 
          std::string region;
          
          if (elem_marker == HERMES_DUMMY_ELEM_MARKER)
            region = this->nu.begin()->first; 
          else
            region = wf->get_element_markers_conversion()->get_user_marker(elem_marker);
          
          RegionMaterialMap::const_iterator material = this->region_material_map.find(region);
          
          if (material != this->region_material_map.end())
            return material->second;
          else
            return region; // Corresponds to the case when region <==> material, 
        }
        
        std::string MaterialPropertyMaps::get_material(int elem_marker, Mesh *mesh) const 
        { 
          std::string region;
          
          if (elem_marker == HERMES_DUMMY_ELEM_MARKER)
            region = this->nu.begin()->first; 
          else
            region = mesh->get_element_markers_conversion().get_user_marker(elem_marker);
          
          RegionMaterialMap::const_iterator material = this->region_material_map.find(region);
          
          if (material != this->region_material_map.end())
            return material->second;
          else
            return region; // Corresponds to the case when region <==> material, 
        }
        
        const rank1& MaterialPropertyMaps::get_Sigma_a(const std::string& material) const
        {
          MaterialPropertyMap1::const_iterator data = this->Sigma_a.find(material);
          if (data != this->Sigma_a.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            string mat = data_elem->first;
            
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
      }
      
      namespace Diffusion
      {        
        void MaterialPropertyMaps::validate()
        {
          Common::MaterialPropertyMaps::validate();
          
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
                  Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_a));
                }
                else 
                {
                  // If only Sigma_t is given, we assume that all reaction terms are included in Sigma_t; all
                  // other x-sections will be set to zero.
                  warning(W_NO_SCATTERING);
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
                warning(W_NO_SCATTERING);
                fill_with(0.0, &Sigma_s);
                Sigma_s_given = true;
              }
              
              if (Sigma_a_given)
                Sigma_t = Common::NDArrayMapOp::add<rank1>(Sigma_a, sum_map2_columns(Sigma_s));
              else 
              {
                // If neither Sigma_r, Sigma_t, Sigma_a are given, we may have a purely fissioning system.
                if (Sigma_f_given)
                  Sigma_t = Sigma_f;
                else
                  error(E_INSUFFICIENT_DATA);
              }
              
              Sigma_t_given = true;
            }
            
            Sigma_r = Common::NDArrayMapOp::subtract<rank1>(Sigma_t, extract_map2_diagonals(Sigma_s));
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
              Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_r));
            }
            else
            {
              warning(W_NO_SCATTERING);
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
            error(E_NONMATCHING_PROPERTIES);
          
          using ValidationFunctors::ensure_size;
          std::for_each(Sigma_s.begin(), Sigma_s.end(), ensure_size(G,G));
          std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size(G));
          std::for_each(D.begin(), D.end(), ensure_size(G));
        }
        
        rank1 MaterialPropertyMaps::compute_Sigma_a(const std::string& material) const
        {
          if (!Sigma_a.empty())
            return Common::MaterialPropertyMaps::get_Sigma_a(material);
          
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
            error(E_INVALID_MARKER);
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
            return Common::NDArrayMapOp::add<rank1>(Sr, Ssd);
          }
          else
          {
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        
        std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
        {
          using namespace std;
          
          os << static_cast<const Common::MaterialPropertyMaps&>(matprop) << endl;
          
          os << setw(12) << "target group" << setw(10) << "D" << setw(10) << "Sigma_r";
          os << setw(22) << "Sigma_s" << endl; 
          
          MaterialPropertyMap1::const_iterator data_elem = matprop.Sigma_r.begin();
          for ( ; data_elem != matprop.Sigma_r.end(); ++data_elem)
          {
            string mat = data_elem->first;
            
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
            error(E_INSUFFICIENT_DATA);
          
          MaterialPropertyMaps::validate();
         
          if (Sigma_t.empty())
            error(E_INSUFFICIENT_DATA);
          
          if (!Sigma_s_1_out_given)
            Sigma_s_1_out = Common::NDArrayMapOp::multiply<rank1>(mu_av, sum_map2_columns(Sigma_s));
            
          MaterialPropertyMap1 Sigma_tr = Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_s_1_out);
          MaterialPropertyMap1::const_iterator Str_elem = Sigma_tr.begin();
          for ( ; Str_elem != Sigma_tr.end(); ++Str_elem)
            for (unsigned int g = 0; g < G; g++)
              D[Str_elem->first][g] = 1./(3.*Str_elem->second[g]);
        }
      }
      
      namespace SPN
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
                  error(E_LAPACK_ERROR, info);
                
                double rcond;
                dgecon_("1", &n, A, &lda, &anorm, &rcond, 
                        dgecon_double_workspace, dgecon_int_workspace, &info);
                if (info != 0) 
                  error(E_LAPACK_ERROR, info);
                if (rcond > 1e12)
                  warning(W_SINGULAR_MATRIX);
                        
                dgetri_(&n, A, &lda, ipiv, dgetri_workspace, &sz_dgetri_workspace, &info);
                if (info != 0)
                  error(E_LAPACK_ERROR);
                                
                rank2::iterator inv_mtx_row = inverted_moment_matrix->begin(); gto = 0;
                for ( ; inv_mtx_row != inverted_moment_matrix->end(); ++inv_mtx_row, ++gto)
                  std::copy(&A[gto*G], &A[gto*(G+1)], inv_mtx_row->begin());
                
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
            error(E_SIGMA_T_REQUIRED);
          
          Common::MaterialPropertyMaps::validate();
                    
          if (Sigma_sn.empty())
          {
            warning(W_NO_SCATTERING);
            fill_with(0.0, &Sigma_sn);
          }
          
          std::for_each(Sigma_tn.begin(), Sigma_tn.end(), ValidationFunctors::ensure_size(G,G,N+1));
          
          MaterialPropertyMap3::const_iterator Stn_material = Sigma_tn.begin();
          for ( ; Stn_material != Sigma_tn.end(); ++Stn_material)
          {
            std::string mat = Stn_material->first;
            rank3 Ssn = Sigma_sn[mat];
            rank3 Stn = Sigma_tn[mat];
              
            Sigma_rn_is_diagonal[mat].resize(N+1, true);
            
            if (Ssn.size() > N+1)
            {
              warning(W_SCATTERING_TRUNCATION, N+1);
              Sigma_sn[mat].erase(Sigma_sn[mat].begin() + N+1, Sigma_sn[mat].end());
              Ssn = Sigma_sn[mat];
            }
            
            if (Ssn.size() == N+1)
            {
              Sigma_rn[mat] = Common::NDArrayMapOp::subtract<rank3>(Stn, Ssn);
              
              rank3 moment_matrices = Stn_material->second;
              rank3::const_iterator moment_matrix = moment_matrices.begin();
              bool1::iterator moment_matrix_is_diagonal = Sigma_rn_is_diagonal[mat].begin();
              for ( ; moment_matrix != moment_matrices.end(); ++moment_matrix, ++moment_matrix_is_diagonal)
              { 
                if (moment_matrix->size() != G)
                  error(E_INVALID_SIZE);
                
                for (unsigned int gto = 0; gto < G && *moment_matrix_is_diagonal; gto++)
                {
                  if (moment_matrix->at(gto).size() != G)
                    error(E_INVALID_SIZE);
                  
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
                Sigma_rn[mat][n] = Common::NDArrayMapOp::subtract<rank2>(Stn[n], Ssn[n]);
                
                bool isdiag = true;
                for (unsigned int gto = 0; gto < G && isdiag; gto++)
                {
                  if (Ssn[n].size() != G)
                    error(E_INVALID_SIZE);
                  
                  for (unsigned int gfrom = 0; gfrom < G && isdiag; gfrom++)
                  {
                    if (Ssn[n][gto].size() != G)
                      error(E_INVALID_SIZE);
                    
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
            return Common::MaterialPropertyMaps::get_Sigma_a(material);
          
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
            error(E_INVALID_MARKER);
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
          
          os << static_cast<const Common::MaterialPropertyMaps&>(matprop) << endl;
          
          int gto_width = 12;
          int elem_width = 14;
          int total_width = 2*gto_width + 2*elem_width*matprop.G;          
                              
          MaterialPropertyMap3::const_iterator Srn_elem = matprop.Sigma_rn.begin();
          MaterialPropertyMap3::const_iterator Srn_inv_elem = matprop.odd_Sigma_rn_inv.begin();
          for ( ; Srn_elem != matprop.Sigma_rn.end(); ++Srn_elem, ++Srn_inv_elem)
          {
            string mat = Srn_elem->first;
            rank3 Srn_moments = Srn_elem->second;
            rank3 Srn_inv_moments = Srn_inv_elem->second;
            
            os << setw(total_width) << setfill('_') << ' ' << endl << setfill(' ');
            os << setw(total_width/2+mat.length()/2) << mat << endl << endl;
            
            os << setw(gto_width)  << "trgt group";
            os << setw(elem_width) << "Sigma_rn";
            os << setw(elem_width) << "odd_Srn_inv";
            
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
                {
                  for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
                    os << setw(elem_width) << (*Srn_inv_moment)[gto][gfrom];
                  
                  ++Srn_inv_moment;
                }
                
                os << endl;
              }
            }
          }
          
          return os << endl;
        }
      }
    }
    
    namespace ElementaryForms
    {  
      namespace SPN
      {        
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = int_u_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          
          return Coeffs::D_grad_F(mrow, mcol) * result;
        }
        template
        scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;                                                               
        
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result = 0;
          
          for (unsigned int mcol = 0; mcol < N_odd; mcol++)
          {
            double coeff = Coeffs::D_grad_F(mrow, mcol);

            unsigned int i = mg.pos(mcol,g);
            
            if (geom_type == HERMES_PLANAR) 
              result += coeff * int_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v);
            else if (geom_type == HERMES_AXISYM_X) 
              result += coeff * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
            else 
              result += coeff * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
          }
          
          return result;
        }
        template
        scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<scalar> *u_ext[],
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Ord> *u_ext[],
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
                                                           
        template<typename Real, typename Scalar>
        Scalar DiagonalStreamingAndReactions::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
        {
          Scalar result;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);     
          
          double Sigma_r_elem = 0.;
          for (unsigned int k = 0; k <= mrow; k++)
            Sigma_r_elem += Coeffs::system_matrix(mrow, mrow, k) * matprop.get_Sigma_rn(mat)[2*k][g][g];
          
          double D_elem = -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][g][g];
          
          // cout << "DiagonalStreamingAndReactions::Jacobian (mom. #" << mrow << ") | " << mat << " | Sigma_r = " << Sigma_r_elem << " | D = " << D_elem << endl;

          if (geom_type == HERMES_PLANAR) 
          {
            result = D_elem * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
                     Sigma_r_elem * int_u_v<Real, Scalar>(n, wt, u, v);
          }
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
            {
              result = D_elem * int_y_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) + 
                       Sigma_r_elem * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            }
            else 
            {
              result = D_elem * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) + 
                       Sigma_r_elem * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
            }
          }
          return result;
        }
        
        template<typename Real, typename Scalar>
        Scalar DiagonalStreamingAndReactions::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
        { 
          Scalar result;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);     
          
          double Sigma_r_elem = 0.;
          for (unsigned int k = 0; k <= mrow; k++)
            Sigma_r_elem += Coeffs::system_matrix(mrow, mrow, k) * matprop.get_Sigma_rn(mat)[2*k][g][g];
          
          double D_elem = -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][g][g];
          
          // cout << "DiagonalStreamingAndReactions::Residual (mom. #" << mrow << ") | " << mat << " | Sigma_r = " << Sigma_r_elem << " | D = " << D_elem << endl;          
          
          unsigned int i = mg.pos(mrow,g);
          
          if (geom_type == HERMES_PLANAR) 
            result = D_elem * int_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[i], v) +
                     Sigma_r_elem * int_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = D_elem * int_y_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[i], v, e) + 
                     Sigma_r_elem * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
          else 
            result = D_elem * int_x_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[i], v, e) + 
                     Sigma_r_elem * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
          
          return result;
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
        {
          if (!matprop.get_fission_nonzero_structure()[gto])
            return 0.0;
          
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = int_u_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          return result * (-Coeffs::system_matrix(mrow, mcol, 0)) * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        {  
          if (!matprop.get_fission_nonzero_structure()[g])
            return 0.0;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
                    
          Scalar result = 0;
          for (int i = 0; i < n; i++) 
          {
            Scalar local_res = 0;
            for (int gfrom = 0; gfrom < ext->nf; gfrom++)
              local_res += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext->fn[gfrom]->val[i]; // scalar flux in group 'gfrom'
            
            // cout << "FissionYield::OuterIterationForm (mom. #" << mrow << " (x, y) = (" << e->x[i] << ", " << e->y[i] << "), " << mat << ") : ";
            // cout << (local_res * Coeffs::even_moment(0, mrow) * chi_elem[g] / keff) << endl;
            
            local_res = local_res * wt[i] * v->val[i];
            
            if (geom_type == HERMES_AXISYM_X)
              local_res = local_res * e->y[i];
            else if (geom_type == HERMES_AXISYM_Y)
              local_res = local_res * e->x[i];
            
            result += local_res;
          }
          
          return result * Coeffs::even_moment(0, mrow) * chi_elem[g] / keff;
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          if (!matprop.get_fission_nonzero_structure()[gto])
            return 0.0;
                   
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          Scalar result = 0;
          for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
          {
            double nSf = nu_elem[gfrom] * Sigma_f_elem[gfrom];
            
            for (unsigned int mcol = 0; mcol <= N_odd; mcol++)
            {
              if (geom_type == HERMES_PLANAR) 
                result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_u_ext_v<Real, Scalar>(n, wt, u_ext[mg.pos(mcol,gfrom)], v);
              else if (geom_type == HERMES_AXISYM_X) 
                result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[mg.pos(mcol,gfrom)], v, e);
              else 
                result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[mg.pos(mcol,gfrom)], v, e);
            }
          }
          
          return result * chi_elem[gto];
        }
        
        template<typename Real, typename Scalar>
        Scalar OffDiagonalStreaming::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const
        {
          if (gfrom == gto)
            return 0;
          
          Scalar result = 0;
          
          if (geom_type == HERMES_PLANAR) 
            result = int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = int_y_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
          
          std::string mat = matprop.get_material(e->elem_marker, wf);     
          
          // cout << "OffDiagonalStreaming::Jacobian (mom. #" << mrow << ") | " << mat << " | D = " << -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto][gfrom] << endl;          
          
          return -result * Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto][gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar OffDiagonalStreaming::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          Scalar result = 0;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 D_elem = matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto];
          
          for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
          { 
            if (gfrom != gto)
            {
              unsigned int i = mg.pos(mrow, gfrom);
              
              if (geom_type == HERMES_PLANAR) 
                result += D_elem[gfrom] * int_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v);
              else if (geom_type == HERMES_AXISYM_X) 
                result += D_elem[gfrom] * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
              else 
                result += D_elem[gfrom] * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[i], v, e);
            }
          }
          
          // cout << "OffDiagonalStreaming::Residual (mom. #" << mrow << ") | " << mat << " | D = " << -Coeffs::D(mrow) * D_elem[0] << endl;          
          
          return -result * Coeffs::D(mrow);
        }
        
        template<typename Real, typename Scalar>
        Scalar OffDiagonalReactions::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const
        {
          if (mrow == mcol)
            return 0;
          
          Scalar result = 0;
                                        
          if (geom_type == HERMES_PLANAR)
            result = int_u_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          
          double Sigma_rn_elem = 0.;
          for (unsigned int k = 0; k <= mrow; k++)
            Sigma_rn_elem += Coeffs::system_matrix(mrow, mcol, k) * matprop.get_Sigma_rn(mat)[2*k][gto][gfrom];
          
          // cout << "OffDiagonalReactions::Jacobian (mom. #(" << mrow << "," << mcol << ") | " << mat << " | Sigma_r = " << Sigma_rn_elem << endl;
          
          return result * Sigma_rn_elem;
        }
        
        template<typename Real, typename Scalar>
        Scalar OffDiagonalReactions::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank3 Sigma_rn_elem = matprop.get_Sigma_rn(mat);
          
          Scalar result = 0;
          unsigned int i = mg.pos(mrow, gto);
          for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
          {            
            for (unsigned int mcol = 0; mcol < N_odd; mcol++)
            {
              unsigned int j = mg.pos(mcol, gfrom);
              
              if (i != j)
              {
                double coeff = 0.;
                for (unsigned int k = 0; k <= std::min(mrow, mcol); k++)
                  coeff += Sigma_rn_elem[2*k][gto][gfrom] * Coeffs::system_matrix(mrow, mcol, k);
                
                // cout << "OffDiagonalReactions::Residual (mom. #(" << mrow << "," << mcol << ") | " << mat << " | coeff = " << coeff << endl;
                
                if (geom_type == HERMES_PLANAR) 
                  result += coeff * int_u_ext_v<Real, Scalar>(n, wt, u_ext[j], v);
                else if (geom_type == HERMES_AXISYM_X) 
                  result += coeff * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[j], v, e);
                else 
                  result += coeff * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[j], v, e);
              }
            }
          }
          
          return result;
        }
        
        template<typename Real, typename Scalar>
        Scalar ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          std::string mat = matprop.get_material(e->elem_marker, wf);
          
          if (geom_type == HERMES_PLANAR) 
            return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_v<Real>(n, wt, v);
          else if (geom_type == HERMES_AXISYM_X) 
            return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_y_v<Real>(n, wt, v, e);
          else 
            return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_x_v<Real>(n, wt, v, e);
        }
      }
      
      namespace Diffusion
      { 
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = 0.5 * int_u_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = 0.5 * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = 0.5 * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          
          return result;
        }
        template
        scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;                                                               
        
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = 0.5 * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = 0.5 * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
          else 
            result = 0.5 * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
          
          return result;
        }
        template
        scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<scalar> *u_ext[],
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Ord> *u_ext[],
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;        
       
        template<typename Real, typename Scalar>
        Scalar DiffusionReaction::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        {
          Scalar result;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);     
          rank1 D_elem = matprop.get_D(mat);
          rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
          
          if (geom_type == HERMES_PLANAR) 
          {
            result = D_elem[g] * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
                      Sigma_r_elem[g] * int_u_v<Real, Scalar>(n, wt, u, v);
          }
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
            {
              result = D_elem[g] * int_y_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) + 
                        Sigma_r_elem[g] * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            }
            else 
            {
              result = D_elem[g] * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) + 
                        Sigma_r_elem[g] * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
            }
          }
          return result;
        }
        
        template<typename Real, typename Scalar>
        Scalar DiffusionReaction::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          std::string mat = matprop.get_material(e->elem_marker, wf);        
          rank1 D_elem = matprop.get_D(mat);
          rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
          
          if (geom_type == HERMES_PLANAR) 
          {
            result = D_elem[g] * int_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v) +
                     Sigma_r_elem[g] * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
          }
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
            {
              result = D_elem[g] * int_y_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) + 
                       Sigma_r_elem[g] * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
            }
            else 
            {
              result = D_elem[g] * int_x_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) + 
                       Sigma_r_elem[g] * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
            }
          }
          return result;
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
        {
          if (!matprop.get_fission_nonzero_structure()[gto])
            return 0.0;
          
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          }
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          if (!matprop.get_fission_nonzero_structure()[g])
            return 0.0;
            
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          if ((unsigned)ext->nf != nu_elem.size() || (unsigned)ext->nf != Sigma_f_elem.size())
            error(E_INVALID_GROUP_INDEX);
          
          Scalar result = 0;
          for (int i = 0; i < n; i++) 
          {
            Scalar local_res = 0;
            for (int gfrom = 0; gfrom < ext->nf; gfrom++)
              local_res += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext->fn[gfrom]->val[i];
                      
            local_res = local_res * wt[i] * v->val[i];
            
            if (geom_type == HERMES_AXISYM_X)
              local_res = local_res * e->y[i];
            else if (geom_type == HERMES_AXISYM_Y)
              local_res = local_res * e->x[i];
            
            result += local_res;
          }
          
          return result * chi_elem[g] / keff;
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          if (!matprop.get_fission_nonzero_structure()[gto])
            return 0.0;
          
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
            else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
          }
          
          std::string mat = matprop.get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar Scattering::Jacobian::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const  
        {
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          }
          
          return result * matprop.get_Sigma_s(matprop.get_material(e->elem_marker, wf))[gto][gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar Scattering::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
            else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
          }
          
          return result * matprop.get_Sigma_s(matprop.get_material(e->elem_marker, wf))[gto][gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          std::string mat = matprop.get_material(e->elem_marker, wf);
          
          if (geom_type == HERMES_PLANAR) 
            return matprop.get_iso_src(mat)[g] * int_v<Real>(n, wt, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
              return matprop.get_iso_src(mat)[g] * int_y_v<Real>(n, wt, v, e);
            else 
              return matprop.get_iso_src(mat)[g] * int_x_v<Real>(n, wt, v, e);
          }
        }
      }
    }
    
    namespace CompleteWeakForms
    {             
      namespace Diffusion
      {   
        void DefaultWeakFormFixedSource::lhs_init(unsigned int G, 
                                                  const MaterialPropertyMaps& matprop, GeomType geom_type)
        {
          bool2 Ss_nnz = matprop.get_scattering_nonzero_structure();
          bool1 chi_nnz = matprop.get_fission_nonzero_structure();
          
          for (unsigned int gto = 0; gto < G; gto++)
          {
            add_matrix_form(new DiffusionReaction::Jacobian(gto, matprop, geom_type));
            add_vector_form(new DiffusionReaction::Residual(gto, matprop, geom_type));
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              if (Ss_nnz[gto][gfrom] && gto != gfrom)
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, matprop, geom_type));
              }
              
              if (chi_nnz[gto])
              {
                add_matrix_form(new FissionYield::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new FissionYield::Residual(gto, gfrom, matprop, geom_type));
              }
            }
          }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                                               GeomType geom_type) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new ExternalSources::LinearForm(gto, matprop, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                HermesFunction *minus_f_src, const std::string& src_area,
                                                                GeomType geom_type  ) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, minus_f_src, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                HermesFunction *minus_f_src,
                                                                const Hermes::vector<std::string>& src_areas,
                                                                GeomType geom_type  ) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, minus_f_src, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                const std::vector<HermesFunction*>& minus_f_src,
                                                                const std::string& src_area, 
                                                                GeomType geom_type ) : WeakForm(matprop.get_G())
        {
          if (minus_f_src.size() != matprop.get_G())
            error(E_INVALID_SIZE);
          
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, minus_f_src[gto], geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                const std::vector<HermesFunction*>& minus_f_src,
                                                                const Hermes::vector<std::string>& src_areas,
                                                                GeomType geom_type ) : WeakForm(matprop.get_G())
        {
          if (minus_f_src.size() != matprop.get_G())
            error(E_INVALID_SIZE);
          
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, minus_f_src[gto], geom_type));
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                                                        const Hermes::vector<MeshFunction*>& iterates,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type ) 
          : WeakForm(matprop.get_G()), Common::WeakFormSourceIteration(initial_keff_guess)
        {      
          init(matprop, iterates, geom_type);
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                                                        const Hermes::vector<Solution*>& iterates,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type ) 
          : WeakForm(matprop.get_G()), Common::WeakFormSourceIteration(initial_keff_guess)
        {      
          Hermes::vector<MeshFunction *> iterates_mf;
          for (unsigned int i = 0; i < iterates.size(); i++)
            iterates_mf.push_back(static_cast<MeshFunction*>(iterates[i]));
          
          init(matprop, iterates_mf, geom_type);
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                                                        const Hermes::vector<MeshFunction*>& iterates, 
                                                                        const Hermes::vector<std::string>& src_areas,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type )
          : WeakForm(matprop.get_G()), Common::WeakFormSourceIteration(initial_keff_guess)
        {
          init(matprop, iterates, geom_type, src_areas);
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                                                        const Hermes::vector<Solution*>& iterates, 
                                                                        const Hermes::vector<std::string>& src_areas,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type ) 
          : WeakForm(matprop.get_G()), Common::WeakFormSourceIteration(initial_keff_guess)
        {
          Hermes::vector<MeshFunction *> iterates_mf;
          for (unsigned int i = 0; i < iterates.size(); i++)
            iterates_mf.push_back(static_cast<MeshFunction*>(iterates[i]));
          
          init(matprop, iterates_mf, geom_type, src_areas);
        }
        
        void DefaultWeakFormSourceIteration::init(const MaterialPropertyMaps& matprop,
                                                  const Hermes::vector<MeshFunction*>& iterates, 
                                                  GeomType geom_type, const Hermes::vector<std::string>& areas)
        {
          bool2 Ss_nnz = matprop.get_scattering_nonzero_structure();
          
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
          {
            add_matrix_form(new DiffusionReaction::Jacobian(gto, matprop, geom_type));
            add_vector_form(new DiffusionReaction::Residual(gto, matprop, geom_type));
            
            for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
            {
              if (Ss_nnz[gto][gfrom] && gto != gfrom)
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, matprop, geom_type));
              }
            }
            
            FissionYield::OuterIterationForm* keff_iteration_form;
            
            if (areas.size() > 0) 
              keff_iteration_form = new FissionYield::OuterIterationForm( gto, areas, matprop, iterates, keff, geom_type );
            else
              keff_iteration_form = new FissionYield::OuterIterationForm( gto, matprop, iterates, keff, geom_type );
            
            keff_iteration_forms.push_back(keff_iteration_form);
            add_vector_form(keff_iteration_form);
          }
        }
        
        void DefaultWeakFormSourceIteration::update_keff(double new_keff) 
        { 
          keff = new_keff;
          
          std::vector<FissionYield::OuterIterationForm*>::iterator it = keff_iteration_forms.begin();
          for ( ; it != keff_iteration_forms.end(); ++it)
            (*it)->update_keff(new_keff); 
        }
      }
    
      namespace SPN
      {
        WeakFormHomogeneous::WeakFormHomogeneous(unsigned int N, const MaterialPropertyMaps& matprop,
                                                 GeomType geom_type, bool include_fission) 
          : WeakForm(matprop.get_G()*(N+1)/2), G(matprop.get_G()), N_odd((N+1)/2)
        {
          mg.set_G(G);
         
          bool1 diagonal_moments = matprop.is_Sigma_rn_diagonal();
          bool2 present(N_odd * G, bool1(N_odd * G, false));
          bool2 sym(N_odd * G, bool1(N_odd * G, false));
          
          for (unsigned int m = 0; m < N_odd; m++)
          {
            for (unsigned int gto = 0; gto < G; gto++)
            {
              unsigned int i = mg.pos(m, gto);
              
              for (unsigned int n = m; n < N_odd; n++)
              {
                for (unsigned int gfrom = 0; gfrom < G; gfrom++)
                {
                  unsigned int j = mg.pos(n, gfrom);
                  
                  if ( (j-i)%G )
                  {
                    for (unsigned int k = 0; k < 2*m; k+=2)
                    {
                      if (!diagonal_moments[k])
                      {
                        present[j][i] = present[i][j] = true;
                        break;
                      }
                    }
                    
                    if (!present[j][i])
                    {
                      if (j < (n+1)*G && diagonal_moments[2*m+1])
                      {
                        present[j][i] = present[i][j] = true;
                      }
                    }
                  }
                  else
                  {
                    present[j][i] = true;
                    sym[j][i] = sym[i][j] = true;
                  }
                }
              }
            }
          }
          
          bool1 chi_nnz = matprop.get_fission_nonzero_structure();
          
          for (unsigned int gto = 0; gto < G; gto++)
          {
            for (unsigned int m = 0; m < N_odd; m++)
            {
              unsigned int i = mg.pos(m, gto);
              
              add_matrix_form(new DiagonalStreamingAndReactions::Jacobian(m, gto, matprop, geom_type));
              add_vector_form(new DiagonalStreamingAndReactions::Residual(m, gto, matprop, geom_type));
              
              if (include_fission && chi_nnz[gto])
                add_vector_form(new FissionYield::Residual(m, N_odd, gto, matprop, geom_type));
              
              add_vector_form(new OffDiagonalReactions::Residual(m, N_odd, gto, matprop, geom_type));
              
              if (G > 1)
                add_vector_form(new OffDiagonalStreaming::Residual(m, gto, matprop, geom_type));
              
              for (unsigned int gfrom = 0; gfrom < G; gfrom++)
              {
                if (gfrom != gto)
                  add_matrix_form(new OffDiagonalStreaming::Jacobian(m, gto, gfrom, matprop, geom_type));
                
                for (unsigned int n = 0; n < N_odd; n++)
                {
                  unsigned int j = mg.pos(n, gfrom);
                  
                  if (include_fission && chi_nnz[gto])
                    add_matrix_form( new FissionYield::Jacobian(m, n, gto, gfrom, matprop, geom_type) );
                  
                  //// cout << "(" << i << "," << j << ") : P" << present[i][j] << " S" << sym[i][j] << endl;
                  
                  if (i != j)
                  {
                    if (present[i][j])
                      add_matrix_form( new OffDiagonalReactions::Jacobian(m, n, gto, gfrom, matprop, geom_type, 
                                                                          sym[i][j] ? HERMES_SYM : HERMES_NONSYM) );
                  }
                }
              }
            }
          }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, unsigned int N,
                                                               GeomType geom_type) 
          : WeakFormHomogeneous(N, matprop, geom_type, true)
        {
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
              add_vector_form(new ExternalSources::LinearForm(m, gto, matprop, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, unsigned int N, 
                                                                HermesFunction *minus_isotropic_source, std::string src_area,
                                                                GeomType geom_type  )
          : WeakFormHomogeneous(N, matprop, geom_type, true)
        {
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
            {
              VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_area, minus_isotropic_source, geom_type);
              src->scaling_factor = Coeffs::even_moment(0, m);
              source_terms.push_back(src);
              add_vector_form(src);
            }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, unsigned int N, 
                                                                HermesFunction *minus_isotropic_source,
                                                                Hermes::vector<std::string> src_areas,
                                                                GeomType geom_type  )
          : WeakFormHomogeneous(N, matprop, geom_type, true)
        {
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
            {
              VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_areas, minus_isotropic_source, geom_type);
              src->scaling_factor = Coeffs::even_moment(0, m);
              source_terms.push_back(src);
              add_vector_form(src);
            }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, unsigned int N, 
                                                                const std::vector<HermesFunction*>& minus_isotropic_sources,
                                                                std::string src_area, 
                                                                GeomType geom_type )
          : WeakFormHomogeneous(N, matprop, geom_type, true)
        {
          if (minus_isotropic_sources.size() != G)
            error(E_INVALID_SIZE);
         
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
            {
              VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_area, minus_isotropic_sources[gto], geom_type);
              src->scaling_factor = Coeffs::even_moment(0, m);
              source_terms.push_back(src);
              add_vector_form(src);
            }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, unsigned int N, 
                                                                const std::vector<HermesFunction*>& minus_isotropic_sources,
                                                                Hermes::vector<std::string> src_areas,
                                                                GeomType geom_type )
          : WeakFormHomogeneous(N, matprop, geom_type, true)
        {
          if (minus_isotropic_sources.size() != G)
            error(E_INVALID_SIZE);
          
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
            {
              VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_areas, minus_isotropic_sources[gto], geom_type);
              src->scaling_factor = Coeffs::even_moment(0, m);
              source_terms.push_back(src);
              add_vector_form(src);
            }
        }
        
        DefaultWeakFormFixedSource::~DefaultWeakFormFixedSource()
        {
          std::vector<VectorFormVol*>::const_iterator it = source_terms.begin();
          for ( ; it != source_terms.end(); ++it)
            delete *it;
          source_terms.clear();
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop, unsigned int N,
                                                                        const Hermes::vector<Solution*>& iterates, 
                                                                        const Hermes::vector<std::string>& src_areas,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type )
          : WeakFormHomogeneous(N, matprop, geom_type, false), Common::WeakFormSourceIteration(initial_keff_guess)
        {      
          SupportClasses::SPN::MomentFilter::get_scalar_fluxes_with_derivatives(iterates, &scalar_flux_iterates, G);
          
          for (unsigned int m = 0; m < N_odd; m++)
          {
            for (unsigned int gto = 0; gto < G; gto++)
            {            
              FissionYield::OuterIterationForm* keff_iteration_form;
              
              if (src_areas.size() > 0) 
                keff_iteration_form = new FissionYield::OuterIterationForm( m, gto, src_areas, matprop, scalar_flux_iterates, 
                                                                            initial_keff_guess, geom_type );
              else
                keff_iteration_form = new FissionYield::OuterIterationForm( m, gto, matprop, scalar_flux_iterates, 
                                                                            initial_keff_guess, geom_type );
              
              keff_iteration_forms.push_back(keff_iteration_form);
              add_vector_form(keff_iteration_form);
            }
          }
        }
        
        DefaultWeakFormSourceIteration::~DefaultWeakFormSourceIteration()
        {
          std::vector<FissionYield::OuterIterationForm*>::const_iterator it = keff_iteration_forms.begin();
          for ( ; it != keff_iteration_forms.end(); ++it)
            delete *it;
          keff_iteration_forms.clear();
          
          SupportClasses::SPN::MomentFilter::clear_scalar_fluxes(&scalar_flux_iterates);
        }
        
        void DefaultWeakFormSourceIteration::update_keff(double new_keff) 
        { 
          keff = new_keff;
          
          std::vector<FissionYield::OuterIterationForm*>::iterator it = keff_iteration_forms.begin();
          for ( ; it != keff_iteration_forms.end(); ++it)
            (*it)->update_keff(new_keff); 
        }
      }
    }
    
    namespace SupportClasses
    {
      namespace Common
      {
        void SourceFilter::filter_fn(int n, Hermes::vector<scalar*> values, scalar* result)
        {
          int marker = this->get_active_element()->marker;
          std::string material = matprop.get_material(marker, mesh);
          
          memset(result, 0, n*sizeof(scalar));
          
          if (markers.empty() || markers.find(marker) != markers.end())
          {           
            rank1 Sigma_f = matprop.get_Sigma_f(material);
            rank1 nu = matprop.get_nu(material);
          
            for (int i = 0; i < n; i++) 
              for (unsigned int j = 0; j < values.size(); j++)
                result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
          }
        }
        
        void SourceFilter::pre_init()
        {
          num = matprop.get_G();
          if(num > 10)
            error("Unable to create an instance of SourceFilter: Hermes is currently able to handle"
                  "only 10 functions in filters.");
          
          for (int i = 0; i < 10; i++)
          {
            item[i] = H2D_FN_VAL & H2D_FN_COMPONENT_0;
            tables[i] = NULL;
            sln[i] = NULL;
          }
          
          num_components = 1;
          have_solutions = false;
        }

        void SourceFilter::assign_solutions(const Hermes::vector< Solution* >& solutions)
        {
          if (solutions.size() != (unsigned) num)
            error("SourceFilter: Number of solutions does not match the size of data.");
         
          free();
          for (int i = 0; i < num; i++)
            sln[i] = solutions[i];
          init();
          post_init();
        }

        void SourceFilter::assign_solutions(const Hermes::vector< MeshFunction* >& solutions)
        {
          if (solutions.size() != (unsigned) num)
            error("SourceFilter: Number of solutions does not match the size of data.");
          
          free();
          for (int i = 0; i < num; i++)
            sln[i] = solutions[i];
          init();
          post_init();
        }
        
        void SourceFilter::post_init()
        {
          set_quad_2d(&g_quad_2d_std);
          have_solutions = true;
          
          std::set<std::string>::const_iterator it = source_regions.begin();
          for ( ; it != source_regions.end(); ++it)
            markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it));
        }
        
        void SourceFilter::set_active_element(Element* e)
        {
          SimpleFilter::set_active_element(e);
          
          order = sln[0]->get_fn_order();
          for (int i = 1; i < num; i++)
            if (sln[i]->get_fn_order() > order)
              order = sln[i]->get_fn_order();
        }

        
        double SourceFilter::integrate(GeomType geom_type)
        {
          if (!have_solutions)
            return 0.0;
                              
          Quad2D* quad = get_quad_2d(); // Needed for h1_integrate_expression.
          double integral = 0.0;
          Element* e;
                    
          for_all_active_elements(e, mesh)
          {
            if (markers.empty() || markers.find(e->marker) != markers.end())
            {
              update_limit_table(e->get_mode());
              this->set_active_element(e);
              RefMap* ru = this->get_refmap();
              int o = this->get_fn_order() + ru->get_inv_ref_order();
              if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
                o++;
              limit_order(o);
              this->set_quad_order(o, H2D_FN_VAL);
              scalar *uval = this->get_fn_values();
              double result = 0.0;
              
              if (geom_type == HERMES_PLANAR)
              {
                h1_integrate_expression(uval[i]);
              }
              else if (geom_type == HERMES_AXISYM_X)
              {
                double* y = ru->get_phys_y(o);
                h1_integrate_expression(y[i] * uval[i]);
              }
              else
              {
                double* x = ru->get_phys_x(o);
                h1_integrate_expression(x[i] * uval[i]);
              }
              
              integral += result;
            }
          }
          
          if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
            integral *= 2*M_PI;
          
          return integral;
        }

      }
      
      namespace SPN
      {
        const double Coeffs::SYSTEM_MATRIX[N_MAX][N_MAX][N_MAX] =            
        {
          {
            {-1., 0, 0, 0, 0}, 
            {2./3., 0, 0, 0, 0}, 
            {-8./15., 0, 0, 0, 0}, 
            {16./35., 0, 0, 0, 0}, 
            {-128./315., 0, 0, 0, 0}
          },
          {
            {-4./9., -5./9., 0, 0, 0}, 
            {16./45., 4./9., 0, 0, 0}, 
            {-32./105., -8./21., 0, 0, 0}, 
            {256./945., 64./189., 0, 0, 0},
            {0, 0, 0, 0, 0}
          },
          {
            {-64./225., -16./45., -9./25., 0, 0}, 
            {128./525., 32./105., 54./175., 0, 0}, 
            {-1024./4725., -256./945., -48./175., 0, 0},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0}
          },
          {
            {-256./1225., -64./245., -324./1225., -13./49., 0}, 
            {2048./11025., 512./2205., 288./1225., 104./441., 0},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0}
          },
          {
            {-16384./99225., -4096./19845., -256./1225., -832./3969., -17./81.},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0}
          }
        };
          
        const double Coeffs::D_GRAD_F[N_MAX][N_MAX] = 
        {
          {-1./2., 1./8., -1./16., 5./128., -7./256.},
          {-7./24, 41./384., -1./16., 131./3072., 0},
          {-407./1920., 233./2560., -1777./30720., 0, 0},
          {3023./17920., 90989./1146880., 0, 0, 0},
          {-1456787./10321920., 0, 0, 0, 0}
        };
        
        const double Coeffs::EVEN_MOMENTS[N_MAX][N_MAX] = 
        {
          {1., -2./3., 8./15., -16./35., 128./315.},
          {1./3., -4./15., 8./35., -64./315., 0},
          {1./5., -6./35., 16./105., 0, 0},
          {1./7., -8./63., 0, 0, 0},
          {1./9., 0, 0, 0, 0}
        };
            
        double Coeffs::system_matrix(unsigned int m, unsigned int n, unsigned int k_of_Sigma_t2k)
        {
          if (m >= N_MAX || n >= N_MAX)
            error("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
                  "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
          
          if (m > n) std::swap(m,n);
          
          if (k_of_Sigma_t2k > n)
            error("In the m-th SPN equation, m = %d, the coefficients at the unknown generalized fluxes involve"
                  "Sigma_tk with k up to %d. Entered k = %d.", n, 2*n, 2*k_of_Sigma_t2k);
          
          return SYSTEM_MATRIX[m][n-m][k_of_Sigma_t2k];
        }
        
        double Coeffs::D_grad_F(unsigned int m, unsigned int n)
        {
          if (m >= N_MAX || n >= N_MAX)
            error("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
                  "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
          
          if (m > n) std::swap(m,n); 
          return D_GRAD_F[m][n-m];
        }
        
        double Coeffs::even_moment(unsigned int m, unsigned int n)
        {
          if (m >= N_MAX)
            error("For the maximum implemented SPN order (N = %d), there are %d even moments."
                  "Tried to access %d. moment.", 2*N_MAX-1, N_MAX, m+1);
          if (n >= N_MAX-m)
            error("The %d. even moment may be expressed in terms of %d odd moments."
                  "Tried to use %d. odd moment.", m, N_MAX-m, n+1);
          if (m > n)
            error("The %d. even moment may be expressed in terms of odd moments starting with %d."
                  "Tried to use moment %d.", m, m, n+1);
                  
          return EVEN_MOMENTS[m][n-m];
        }
                
        MomentFilter::Common::Common(unsigned int angular_moment, unsigned int group, unsigned int G) 
          : odd_req_mom(false), g(group), mg(G)
        {
          if (angular_moment % 2)
          { 
            warning("Using MomentFilter to access the odd moments of flux is inefficient. "
                    "The vector of solutions itself contains these moments (for each energy group).");
            req_mom_idx = (angular_moment-1)/2;
            odd_req_mom = true;
          }
          else
            req_mom_idx = angular_moment/2;
        }
        
        void MomentFilter::Val::filter_fn(int n, Hermes::vector< scalar* > values, scalar* result)
        {
          if (odd_req_mom)
            memcpy(result, values.at(req_mom_idx), n*sizeof(scalar));
          else
          {            
            for (int i = 0; i < n; i++) 
            {
              result[i] = 0;
              unsigned int exp_mom_idx = req_mom_idx;
              for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
                result[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * values.at(sol_idx)[i];
            }
          }
        }
        
        void MomentFilter::Val::set_active_element(Element* e)
        {
          SimpleFilter::set_active_element(e);

          order = -1;
          
          unsigned int exp_mom_idx = req_mom_idx;
          for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
            if (sln[sol_idx]->get_fn_order() > order)
              order = sln[sol_idx]->get_fn_order();
        }
        
        void MomentFilter::ValDxDy::filter_fn(int n, Hermes::vector<scalar *> values, 
                                              Hermes::vector<scalar *> dx, Hermes::vector<scalar *> dy, 
                                              scalar* rslt, scalar* rslt_dx, scalar* rslt_dy)
        {
          if (odd_req_mom)
          {
            memcpy(rslt, values.at(req_mom_idx), n*sizeof(scalar));
            memcpy(rslt_dx, dx.at(req_mom_idx), n*sizeof(scalar));
            memcpy(rslt_dy, dy.at(req_mom_idx), n*sizeof(scalar));
          }
          else
          {            
            for (int i = 0; i < n; i++) 
            {
              rslt[i] = rslt_dx[i] = rslt_dy[i] = 0;
              unsigned int exp_mom_idx = req_mom_idx;
              for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
              {
                rslt[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * values.at(sol_idx)[i];
                rslt_dx[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * dx.at(sol_idx)[i];
                rslt_dy[i] += Coeffs::even_moment(req_mom_idx, exp_mom_idx) * dy.at(sol_idx)[i];
              }
            }
          }
        }
        
        void MomentFilter::ValDxDy::set_active_element(Element* e)
        {
          DXDYFilter::set_active_element(e);
          
          order = -1;
          
          unsigned int exp_mom_idx = req_mom_idx;
          for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
            if (sln[sol_idx]->get_fn_order() > order)
              order = sln[sol_idx]->get_fn_order();
        }
        
        // TODO: Templatize.
        void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution* >& angular_fluxes, 
                                             Hermes::vector< MeshFunction* >* scalar_fluxes,
                                             unsigned int G)
        {          
          scalar_fluxes->reserve(G);
          for (unsigned int g = 0; g < G; g++)
            scalar_fluxes->push_back(new MomentFilter::Val(0, g, G, angular_fluxes));
        }
        
        void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution* >& angular_fluxes, 
                                                              Hermes::vector< MeshFunction* >* scalar_fluxes,
                                                              unsigned int G)
        {          
          scalar_fluxes->reserve(G);
          for (unsigned int g = 0; g < G; g++)
            scalar_fluxes->push_back(new MomentFilter::ValDxDy(0, g, G, angular_fluxes));
        }
        
        void MomentFilter::clear_scalar_fluxes(Hermes::vector< MeshFunction* >* scalar_fluxes)
        {
          Hermes::vector< MeshFunction* >::const_iterator it = scalar_fluxes->begin();
          for( ; it != scalar_fluxes->end(); ++it)
            delete *it;
          scalar_fluxes->clear();
        }
        
        void SourceFilter::filter_fn(int n, Hermes::vector< scalar* > values, scalar* result)
        {
          int marker = this->get_active_element()->marker;
          std::string material = matprop.get_material(marker, mesh);
          
          memset(result, 0, n*sizeof(scalar));
          
          if (markers.empty() || markers.find(marker) != markers.end())
          {           
            rank1 Sigma_f = matprop.get_Sigma_f(material);
            rank1 nu = matprop.get_nu(material);
            
            for (int i = 0; i < n; i++) 
            {
              for (unsigned int g = 0; g < G; g++)
              {
                scalar group_scalar_flux = 0;
                
                unsigned int exp_mom_idx = 0;
                for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
                  group_scalar_flux += Coeffs::even_moment(0, exp_mom_idx) * values.at(sol_idx)[i];
                
                result[i] += nu[g] * Sigma_f[g] * group_scalar_flux;
              }
            }
          }
        }
      }
      
      double PostProcessor::integrate(MeshFunction* solution, const Hermes::vector<std::string>& areas) const
      {
        Quad2D* quad = &g_quad_2d_std;
        solution->set_quad_2d(quad);
        Mesh* mesh = solution->get_mesh();
        
        std::set<int> markers;
        Hermes::vector<std::string>::const_iterator it = areas.begin();
        for ( ; it != areas.end(); ++it)
          markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it));
        
        double integral = 0.0;
        Element* e;
        for_all_active_elements(e, mesh)
        {
          if (markers.empty() || markers.find(e->marker) != markers.end())
          {
            update_limit_table(e->get_mode());
            solution->set_active_element(e);
            RefMap* ru = solution->get_refmap();
            int o = solution->get_fn_order() + ru->get_inv_ref_order();
            limit_order(o);
            solution->set_quad_order(o, H2D_FN_VAL);
            scalar *uval = solution->get_fn_values();
            double result = 0.0;
            
            if (geom_type == HERMES_PLANAR)
            {
              h1_integrate_expression(uval[i]);
            }
            else if (geom_type == HERMES_AXISYM_X)
            {
              double* y = ru->get_phys_y(o);
              h1_integrate_expression(y[i] * uval[i]);
            }
            else
            {
              double* x = ru->get_phys_x(o);
              h1_integrate_expression(x[i] * uval[i]);
            }
            
            integral += result;
          }
        }
        
        if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
          integral *= 2.0*M_PI;
        
        return integral;
      }

      void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution* >* solutions, double integrated_fission_source) const
      {
        if (integrated_fission_source < 1e-12)
          error("PostProcessor::normalize_to_unit_fission_source : Invalid fission source.");
        
        Hermes::vector< Solution* >::iterator sln = solutions->begin();
        for ( ; sln != solutions->end(); ++sln)
          (*sln)->multiply(1./integrated_fission_source);
      }

      void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution* >* solutions, 
                                                           const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                           const Hermes::vector< std::string >& src_areas) const
      {
        Common::SourceFilter *sf;
        
        if (method == NEUTRONICS_DIFFUSION)
          sf = new Common::SourceFilter(*solutions, matprop, src_areas);
        else if (method == NEUTRONICS_SPN)
          sf = new SPN::SourceFilter(*solutions, matprop, src_areas);
        
        normalize_to_unit_fission_source(solutions, sf->integrate(geom_type));
        
        delete sf;
      }

      void PostProcessor::normalize_to_unit_power(Hermes::vector< Solution* >* solutions, 
                                                  const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                  double power_per_fission, const Hermes::vector< std::string >& src_areas) const
      {
        // TODO
      }
      
      double PostProcessor::get_integrated_group_reaction_rates_internal( ReactionType reaction, MeshFunction* solution, 
                                                                          const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                                          const Hermes::vector< string >& regions,
                                                                          unsigned int this_group, int other_group) const
      {
        using namespace MaterialProperties::Messages;
        
        if (this_group > matprop.get_G())
          error(E_INVALID_GROUP_INDEX);
        
        Quad2D* quad = &g_quad_2d_std;
        solution->set_quad_2d(quad);
        Mesh* mesh = solution->get_mesh();
        
        std::set<int> markers;
        Hermes::vector<std::string>::const_iterator it = regions.begin();
        for ( ; it != regions.end(); ++it)
          markers.insert(mesh->get_element_markers_conversion().get_internal_marker(*it));
                
        double integral = 0.0;
        Element* e;
        for_all_active_elements(e, mesh)
        {
          if (markers.empty() || markers.find(e->marker) != markers.end())
          {
            update_limit_table(e->get_mode());
            solution->set_active_element(e);
            RefMap* ru = solution->get_refmap();
            int o = solution->get_fn_order() + ru->get_inv_ref_order();
            limit_order(o);
            solution->set_quad_order(o, H2D_FN_VAL);
            scalar *uval = solution->get_fn_values();
            double result = 0.0;
            
            if (geom_type == HERMES_PLANAR)
            {
              h1_integrate_expression(uval[i]);
            }
            else if (geom_type == HERMES_AXISYM_X)
            {
              double* y = ru->get_phys_y(o);
              h1_integrate_expression(y[i] * uval[i]);
            }
            else
            {
              double* x = ru->get_phys_x(o);
              h1_integrate_expression(x[i] * uval[i]);
            }
            
            std::string mat = matprop.get_material(e->marker, mesh);
            double xsec;
            
            switch (reaction)
            {
              case ABSORPTION:
                xsec = matprop.compute_Sigma_a(mat)[this_group];
                break;
              case TOTAL:
                xsec = matprop.compute_Sigma_t(mat)[this_group];
                break;
              case IN_SCATTERING:
                if (other_group > (int) matprop.get_G() || other_group < 0)
                  error(E_INVALID_GROUP_INDEX);
                
                xsec = matprop.compute_Sigma_s(mat)[this_group][other_group];
                break;
              case SELF_SCATTERING:
                xsec = matprop.compute_Sigma_s(mat)[this_group][this_group];
                break;
              case OUT_SCATTERING:
                if (other_group > (int) matprop.get_G() || other_group < 0)
                  error(E_INVALID_GROUP_INDEX);
                
                xsec = matprop.compute_Sigma_s(mat)[other_group][this_group];
                break;
              case FISSION:
                xsec = matprop.get_Sigma_f(mat)[this_group];
                break;
              case NU_FISSION:
                xsec = matprop.get_Sigma_f(mat)[this_group] * matprop.get_nu(mat)[this_group];
                break;
            };
            
            integral += xsec * result;
          }
        }
        
        if (geom_type == HERMES_AXISYM_X || geom_type == HERMES_AXISYM_Y)
          integral *= 2.0*M_PI;
        
        return integral;
      }

      
      void PostProcessor::get_integrated_group_reaction_rates(ReactionType reaction, 
                                                              const Hermes::vector< Solution* >& solutions, Hermes::vector< double >* results, 
                                                              const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                              unsigned int group, const Hermes::vector< std::string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        Hermes::vector<std::string>::const_iterator region = regions.begin(); 

        for ( ; region != regions.end(); ++region)
        {
          double result = 0.0;
         
          if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
          {
            for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
            {
              if (reaction == IN_SCATTERING && group != g_other)
                result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, *region, group, g_other);
              else if (reaction == OUT_SCATTERING && group != g_other)
                result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group, g_other);
            }
          }
          else
            result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group);
          
          results->push_back(result);
        }
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
      }
      
      double PostProcessor::get_integrated_group_reaction_rates(ReactionType reaction, const Hermes::vector< Solution* >& solutions, 
                                                                const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                                unsigned int group, const Hermes::vector< string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        double result = 0.0;
          
        if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
        {
          for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
          {
            if (reaction == IN_SCATTERING && group != g_other)
              result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, regions, group, g_other);
            else if (reaction == OUT_SCATTERING && group != g_other)
              result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, regions, group, g_other);
          }
        }
        else
          result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, regions, group);
          
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
        
        return result;
      }
      
      void PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution* >& solutions, 
                                                        Hermes::vector< double >* results, 
                                                        const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                        const Hermes::vector< string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        Hermes::vector<std::string>::const_iterator region = regions.begin(); 
        
        for ( ; region != regions.end(); ++region)
        {
          double result = 0.0;
          for (unsigned int group = 0; group < matprop.get_G(); group++)
          {
            if (reaction == IN_SCATTERING || reaction == OUT_SCATTERING)
            {
              for (unsigned int g_other = 0; g_other < matprop.get_G(); g_other++)
              {
                if (reaction == IN_SCATTERING && group != g_other)
                  result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[g_other], matprop, *region, group, g_other);
                else if (reaction == OUT_SCATTERING && group != g_other)
                  result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group, g_other);
              }
            }
            else
              result += get_integrated_group_reaction_rates_internal(reaction, scalar_fluxes[group], matprop, *region, group);
          }
            
          results->push_back(result);
        }

        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
      }
      
      double PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution* >& solutions, 
                                                          const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                          const Hermes::vector< string >& regions) const
      {
        double result = 0.0;
        for (unsigned int group = 0; group < matprop.get_G(); group++)
          result += get_integrated_group_reaction_rates(reaction, solutions, matprop, group, regions);
        return result;
      }

      void PostProcessor::get_integrated_group_scalar_fluxes( const Hermes::vector< Solution* >& solutions, 
                                                              Hermes::vector< double >* results,
                                                              unsigned int group, unsigned int G,
                                                              const Hermes::vector< std::string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        Hermes::vector<std::string>::const_iterator region = regions.begin(); 
        for ( ; region != regions.end(); ++region)
          results->push_back(integrate(scalar_fluxes[group], *region));
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
      }
      
      double PostProcessor::get_integrated_group_scalar_fluxes(const Hermes::vector< Solution* >& solutions, 
                                                               unsigned int group, unsigned int G,
                                                               const Hermes::vector< std::string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        double result = integrate(scalar_fluxes[group], regions);
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
        
        return result;
      }
      
      void PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution* >& solutions, Hermes::vector< double >* results, 
                                                       unsigned int G, const Hermes::vector< std::string >& regions) const
      {
        Hermes::vector<MeshFunction*> scalar_fluxes;
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
        else if (method == NEUTRONICS_DIFFUSION)
          for (Hermes::vector<Solution*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
            scalar_fluxes.push_back(*it);
        
        Hermes::vector<std::string>::const_iterator region = regions.begin(); 
        for ( ; region != regions.end(); ++region)
        {
          double result = 0.0;
          for (unsigned int group = 0; group < G; group++)
            result += integrate(scalar_fluxes[group], *region);
          results->push_back(result);
        }
        
        if (method == NEUTRONICS_SPN)
          SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
      }

      double PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution* >& solutions, unsigned int G,
                                                         const Hermes::vector< string >& regions) const
      {
        double result = 0.0;
        for (unsigned int group = 0; group < G; group++)
          result += get_integrated_group_scalar_fluxes(solutions, group, G, regions);
        return result;
      }
      
      void PostProcessor::get_areas(Mesh *mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const
      {
        Solution unity(mesh, 1);
        
        Hermes::vector<std::string>::const_iterator region = regions.begin(); 
        for ( ; region != regions.end(); ++region)
          results->push_back(integrate(&unity, *region));
      }

      double PostProcessor::get_area(Mesh* mesh, const Hermes::vector< std::string >& regions) const
      {
        Solution unity(mesh, 1);
        return integrate(&unity, regions);
      }




      SourceIteration::SourceIteration( NeutronicsMethod method, const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                        const Hermes2D& hermes2d, const std::vector< string >& fission_regions, GeomType geom_type )
        : hermes2d(hermes2d), method(method), fission_regions(fission_regions), geom_type(geom_type)
      {
        if (method == NEUTRONICS_DIFFUSION)
        {
          new_source = new Common::SourceFilter(matprop, fission_regions);
          old_source = new Common::SourceFilter(matprop, fission_regions);
        }
        else if (method == NEUTRONICS_SPN)
        {
          new_source = new SPN::SourceFilter(matprop, fission_regions);
          old_source = new SPN::SourceFilter(matprop, fission_regions);
        }
      }
      
      int SourceIteration::eigenvalue_iteration(const Hermes::vector<Solution *>& solutions, DiscreteProblem& dp,
                                                double tol_keff, double tol_flux, MatrixSolverType matrix_solver)
      {
        // Sanity checks.
        if (dp.get_spaces().size() != solutions.size()) 
          error("Spaces and solutions supplied to power_iteration do not match.");
        
        CompleteWeakForms::Common::WeakFormSourceIteration *wf;
        
        if (method == NEUTRONICS_DIFFUSION)
          wf = static_cast<CompleteWeakForms::Diffusion::DefaultWeakFormSourceIteration*>(dp.get_weak_formulation());
        else if (method == NEUTRONICS_SPN)
          wf = static_cast<CompleteWeakForms::SPN::DefaultWeakFormSourceIteration*>(dp.get_weak_formulation());
                                  
        // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
        // updating the eigenvalue. 
        bool meshes_changed;
        Hermes::vector<Solution*> new_solutions;
        for (unsigned int i = 0; i < solutions.size(); i++) 
        {
          Mesh *m = dp.get_space(i)->get_mesh();
          if (solutions[i]->get_mesh()->get_seq() != m->get_seq())
             meshes_changed = true;
          new_solutions.push_back(new Solution(m));
        }
        
        // Assign the solution vectors to fission source calculators.
        new_source->assign_solutions(new_solutions);
        old_source->assign_solutions(solutions);
        
        // Create a matrix solver according to the user's preference. Note that the matrix doesn't change within the 
        // power iteration loop, so the first computed LU factorization may be completely reused in following iterations.
        SparseMatrix* mat = create_matrix(matrix_solver);
        Vector* rhs = create_vector(matrix_solver);
        Solver* solver = create_linear_solver(matrix_solver, mat, rhs);
        solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
        
        // Initial coefficient vector for the Newton's method.
        int ndof = Space::get_num_dofs(dp.get_spaces());
        scalar* coeff_vec = new scalar[ndof];
        
        // Force the Jacobian assembling in the first iteration.
        bool Jacobian_changed = true;
        
        // NOTE:
        //  According to the 'physical' definitions of keff as a fraction of neutron sources and sinks, keff should be 
        //  determined in an iteration like this: 
        //    k_new = k_old * src_new/src_old.  (1)
        //  This iteration can be expanded into the following form:
        //    k_new = k_0 * src_new/src_0,      (2)
        //  so if k_0 is set to src_0 at the beginning, we could potentially save one integration (of src_old) in each 
        //  iteration. Although the results are the same for both approaches (1) and (2) (up to a different eigenvector
        //  scaling), the calculation is slower surprisingly in case (2), even if compared to (1) where both old and new 
        //  sources are integrated in each iteration (this is actually not neccessary as you can see in the implementation
        //  below).
        //        
        double src_new = old_source->integrate(geom_type);
        bool eigen_done = false; int it = 0;
        do 
        {
          memset(coeff_vec, 0.0, ndof*sizeof(scalar));

          if (!hermes2d.solve_newton(coeff_vec, &dp, solver, mat, rhs, Jacobian_changed, 1e-8, 3, true)) 
            error("Newton's iteration failed.");
              
          // The matrix doesn't change within the power iteration loop, so it does not need to be reassembled again.
          Jacobian_changed = false;
          
          // Convert coefficients vector into a set of Solution pointers.
          Solution::vector_to_solutions(solver->get_solution(), dp.get_spaces(), new_solutions);

          // Compute the eigenvalue for current iteration.
          double src_old = src_new;
          src_new = new_source->integrate(geom_type);
          double k_new = wf->get_keff() * src_new / src_old;
          
          double diff_keff = fabs(wf->get_keff() - k_new) / k_new;
          info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, diff_keff);
          
          // Stopping criteria.
          if (diff_keff < tol_keff) 
            eigen_done = true;

          // cout << "Iteration: " << it << ", flux diffces: " << endl; 
          
          if (tol_flux > 0)
          {
            for (unsigned int i = 0; i < solutions.size(); i++)
            {
              PostProcessor pp(method, geom_type);
              
              // Normalize both flux iterates with the same criterion (unit integrated fission source).
              solutions[i]->multiply(1./src_old);
              new_solutions[i]->multiply(1./src_new);
              
              //cout << i << " = " << fabs(pp.integrate(new_solutions[i]) - pp.integrate(solutions[i])) << ", ";
              
              // Compare the two solutions.
              if (fabs(pp.integrate(new_solutions[i]) - pp.integrate(solutions[i])) >= tol_flux * pp.integrate(new_solutions[i]))
                eigen_done = false;
              
              // Go back to unnormalized fluxes.
              solutions[i]->multiply(src_old);
              new_solutions[i]->multiply(src_new);
              
              if (eigen_done == false)
                break;
            }
          }
          
          //cout << endl;

          // Update the final eigenvalue.
          wf->update_keff(k_new);
    
          // Store the new eigenvector approximation in the result.
          for (unsigned int i = 0; i < solutions.size(); i++)  
            solutions[i]->copy(new_solutions[i]);
          
          if (meshes_changed)
          {
            old_source->assign_solutions(solutions);
            // meshes_changed = false; should be here in principle, since meshes do not change
            // during the iteration, but it fails after some iteration
          }
          
          it++;
        }
        while (!eigen_done);
        
        // Free memory.
        delete mat;
        delete rhs;
        delete solver;
        for (unsigned int i = 0; i < solutions.size(); i++) 
          delete new_solutions[i];
        
        return it;
      }
      
      namespace Common
      {
        const std::string Views::base_title_flux = "Neutron flux: group ";
        const std::string Views::base_title_order = "Polynomial orders: group ";
        const std::string Views::base_title_mesh = "Core mesh for group ";
        
        void Views::init(unsigned int nu, unsigned int ne, unsigned int ng) 
        {
          n_unknowns = nu; n_equations = ne; n_groups = ng;
          
          if (nu > 0 && ne > 0 && ng > 0)
          {
            sviews = new ScalarView* [n_unknowns];
            oviews = new OrderView* [n_equations];
            if (display_meshes)
              mviews = new MeshView* [n_equations];
          }
          else
          {
            sviews = NULL;
            oviews = NULL;
            mviews = NULL;
          }
        }
        
        Views::~Views()
        {
          if (sviews != NULL)
          {
            for (unsigned int i = 0; i < n_unknowns; i++)
              delete sviews[i];
            delete [] sviews;
          }
          
          if (oviews != NULL)
          {
            for (unsigned int i = 0; i < n_equations; i++)
              delete oviews[i];
            delete [] oviews;
          }
          
          if (mviews != NULL)
          {
            for (unsigned int i = 0; i < n_equations; i++)
              delete mviews[i];
            delete [] mviews;
          }
        }
        
        void Views::inspect_meshes(Hermes::vector< Mesh* > meshes)
        {
          if (display_meshes)
          {
            show_meshes(meshes);
            View::wait();
            
            for (unsigned int i = 0; i < n_equations; i++)
              delete mviews[i];
            delete [] mviews;
            
            mviews = NULL;
          }
        }
        
        void Views::inspect_solutions(Hermes::vector< Solution* > solutions)
        {
          show_solutions(solutions);
          View::wait();
          
          for (unsigned int i = 0; i < n_unknowns; i++)
            delete sviews[i];
          delete [] sviews;
          
          sviews = NULL;
        }
        
        void Views::inspect_orders(Hermes::vector< Space* > spaces)
        {
          show_orders(spaces);
          View::wait();
          
          for (unsigned int i = 0; i < n_equations; i++)
            delete oviews[i];
          delete [] oviews;
          
          oviews = NULL;
        }
      }
      
      namespace Diffusion
      {
        Views::Views(unsigned int G, bool display_meshes) : Common::Views(G, G, G, display_meshes)
        {
          for (unsigned int g = 0; g < n_groups; g++)
          {
            std::string title_flux = base_title_flux + itos(g);
            std::string title_order = base_title_order + itos(g);
            
            sviews[g] = new ScalarView(title_flux.c_str(), new WinGeom(0, g*452, 450, 450));
            sviews[g]->show_mesh(false);
            sviews[g]->set_3d_mode(true);
            oviews[g] = new OrderView(title_order.c_str(), new WinGeom(0, n_groups*452 + g*452, 450, 450));
          }
          
          if (display_meshes)
            for (unsigned int g = 0; g < n_groups; g++)
            {
              std::string title = base_title_mesh + itos(g);
              mviews[g] = new MeshView(title.c_str(), new WinGeom(0, g*352, 350, 350));
            }
        }
        
        void Views::show_meshes(Hermes::vector< Mesh* > meshes)
        {
          if (display_meshes)
            for (unsigned int g = 0; g < n_groups; g++)
              mviews[g]->show(meshes[g]);
        }
        
        void Views::show_solutions(Hermes::vector< Solution* > solutions)
        {
          for (unsigned int g = 0; g < n_groups; g++)
            sviews[g]->show(solutions[g]);
        }
        
        void Views::show_orders(Hermes::vector< Space* > spaces)
        {
          for (unsigned int g = 0; g < n_groups; g++)
            oviews[g]->show(spaces[g]);
        }
        
        void Views::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                       Hermes::vector< Solution* > solutions, bool mode_3D)
        {
          Linearizer lin;
          for (unsigned int g = 0; g < n_groups; g++)
          {
            std::string appendix = std::string("_group_") + itos(g);
            std::string file = base_filename + appendix + std::string(".vtk");
            std::string var = base_varname + appendix;
            lin.save_solution_vtk(solutions[g], file.c_str(), var.c_str(), mode_3D);
            info("Scalar flux in group %d saved in VTK format to file %s.", g, file.c_str());
          }
        }
        
        void Views::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space* > spaces)
        {
          Orderizer ord;
          for (unsigned int g = 0; g < n_groups; g++)
          {
            std::string file = base_filename + std::string("_group_") + itos(g) + std::string(".vtk");
            ord.save_orders_vtk(spaces[g], file.c_str());
            info("Information about approximation space for group %d saved in VTK format to file %s.", g, file.c_str());
          }
        }
      }
      
      namespace SPN
      {
        Views::Views(unsigned int spn_order, unsigned int G, bool display_meshes) : Common::Views(display_meshes), mg(G)
        {
          n_moments = spn_order+1;
          n_odd_moments = (n_moments+1)/2;
          Common::Views::init(G * n_moments, G * n_odd_moments, G);
          
          for (unsigned int g = 0; g < n_groups; g++)
          {
            std::string title_flux = base_title_flux + itos(g) + std::string(", moment ");
            std::string title_order = base_title_order + itos(g) + std::string(", moment ");
            for (unsigned int m = 0; m < n_moments; m++)
            {
              unsigned int i = mg.pos(m,g);
              
              sviews[i] = new ScalarView((title_flux + itos(m)).c_str(), new WinGeom(m*452, g*452, 450, 450));
              sviews[i]->show_mesh(false);
              sviews[i]->set_3d_mode(true);
              
              if (m%2) 
                oviews[mg.pos((m-1)/2,g)] = new OrderView((title_order + itos(m)).c_str(), new WinGeom(m*452, n_groups*452 + g*452, 450, 450));
            }
          }
          
          if (display_meshes)
          {
            for (unsigned int g = 0; g < n_groups; g++)
            {
              std::string title = base_title_mesh + itos(g) + std::string(", moment ");
              for (unsigned int m = 0; m < n_odd_moments; m++)
                mviews[mg.pos(m,g)] = new MeshView((title + itos(m)).c_str(), new WinGeom(m*352, g*352, 350, 350));
            }
          }
          
          for (unsigned int i = 0; i < MAX_SOLUTIONS_SETS; i++)
          {
            moment_filters[i] = new MomentFilter::Val** [n_odd_moments];
            for (unsigned int m = 0; m < n_odd_moments; m++)
            {
              moment_filters[i][m] = new MomentFilter::Val* [n_groups];
              
              for (unsigned int g = 0; g < n_groups; g++)
                moment_filters[i][m][g] = NULL;
            }
          }
        }
        
        Views::~Views()
        { 
          for (unsigned int i = 0; i < MAX_SOLUTIONS_SETS; i++)
          {
            for (unsigned int m = 0; m < n_odd_moments; m++)
            {
              for (unsigned int g = 0; g < n_groups; g++)
                if (moment_filters[i][m][g])
                  delete moment_filters[i][m][g];
                
              delete [] moment_filters[i][m];
            }
            delete [] moment_filters[i];
          }
        }
        
        void Views::show_meshes(Hermes::vector< Mesh* > meshes)
        {
          if (display_meshes)
            for (unsigned int g = 0; g < n_groups; g++)
              for (unsigned int m = 0; m < n_odd_moments; m++)
                mviews[mg.pos(m,g)]->show(meshes[mg.pos(m,g)]);
        }
        
        void Views::show_solutions_internal(Hermes::vector< Solution* > solutions, unsigned int solutions_set)
        {
          if (solutions_set > MAX_SOLUTIONS_SETS)
            error("Change Views::MAX_SOLUTIONS_SETS and rebuild to allow visualizing more than 10 solutions sets.");
          
          for (unsigned int g = 0; g < n_groups; g++)
          {
            for (unsigned int m = 0; m < n_odd_moments; m++)
            {
              unsigned int i = mg.pos(m,g);
              unsigned int j = mg.pos(2*m,g);
              unsigned int k = mg.pos(2*m+1,g);
              
              if (moment_filters[solutions_set][m][g] == NULL)
                moment_filters[solutions_set][m][g] = new MomentFilter::Val(2*m, g, n_groups, solutions);
              else
                moment_filters[solutions_set][m][g]->reinit();
              
              sviews[j]->show(moment_filters[solutions_set][m][g]);
              sviews[k]->show(solutions[i]);
            }
          }
        }
        
        void Views::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                       Hermes::vector< Solution* > solutions, unsigned int solutions_set, bool mode_3D)
        {
          if (solutions_set > MAX_SOLUTIONS_SETS)
            error("Change Views::MAX_SOLUTIONS_SETS and rebuild to allow visualizing more than 10 solutions sets.");
          
          Linearizer lin;
          for (unsigned int g = 0; g < n_groups; g++)
          {
            std::string appendix = std::string("_group_") + itos(g);
            for (unsigned int m = 0; m < n_odd_moments; m++)
            {
              if (moment_filters[solutions_set][m][g] == NULL)
                moment_filters[solutions_set][m][g] = new MomentFilter::Val(2*m, g, n_groups, solutions);
              
              std::string file = base_filename + std::string("_moment_") + itos(2*m) + appendix + std::string(".vtk");
              std::string var = base_varname + std::string("_moment_") + itos(2*m) + appendix;
              lin.save_solution_vtk(moment_filters[solutions_set][m][g], file.c_str(), var.c_str(), mode_3D);
              info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m, g, file.c_str());
              
              file = base_filename + std::string("_moment_") + itos(2*m+1) + appendix + std::string(".vtk");
              var = base_varname + std::string("_moment_") + itos(2*m+1) + appendix;
              lin.save_solution_vtk(solutions[mg.pos(m,g)], file.c_str(), var.c_str(), mode_3D);
              info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m+1, g, file.c_str());
            }
          }
        }
        
        void Views::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space* > spaces)
        {
          Orderizer ord;
          for (unsigned int g = 0; g < n_groups; g++)
            for (unsigned int m = 0; m < n_odd_moments; m++)
            {
              std::string file = base_filename + std::string("_moment_") + itos(2*m+1) + std::string("_group_") + itos(g) + std::string(".vtk");
              ord.save_orders_vtk(spaces[mg.pos(m,g)], file.c_str());
              info("Information about approximation space for moment %d, group %d saved in VTK format to file %s.", 2*m+1, g, file.c_str());
            }
        }
        
        void Views::inspect_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set)
        {
          show_solutions(solutions, solutions_set);
          View::wait();
          
          for (unsigned int i = 0; i < n_unknowns; i++)
            delete sviews[i];
          delete [] sviews;
          
          sviews = NULL;
        }

        void Views::show_orders(Hermes::vector< Space* > spaces)
        {
          for (unsigned int g = 0; g < n_groups; g++)
            for (unsigned int m = 0; m < n_odd_moments; m++)
              oviews[mg.pos(m,g)]->show(spaces[mg.pos(m,g)]);
        }
      }
    }
  }
}
