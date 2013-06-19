//
// TODO: -  Unify the way MaterialPropertyMaps enter the diffusion and SPN forms (is the static casting in diffusion really
//          neccessary?).
//       -  Unify the way SPN coefficients from class Coeffs are used in weakforms and weakform_parts (some are here, some 
//          in the other file).
//
#include "neutronics/weakforms.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics
{
  namespace SimpleMonoenergeticDiffusionWeakForms
  {    
    FixedSourceProblem::FixedSourceProblem(Hermes::vector<std::string> regions, 
                                           Hermes::vector<double> D_map, 
                                           Hermes::vector<double> Sigma_a_map, 
                                           Hermes::vector<double> Q_map ) : WeakForm<double>(1) 
    {
      using namespace WeakFormsH1;
      
      for (unsigned int i = 0; i < regions.size(); i++)
      {
        /* Jacobian */
        // Diffusion.
        add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, regions[i], new Hermes1DFunction<double>(D_map[i]), 
                                                      HERMES_SYM));
        // Absorption.
        add_matrix_form(new DefaultMatrixFormVol<double>(0, 0, regions[i], new Hermes2DFunction<double>(Sigma_a_map[i]), 
                                                  HERMES_SYM));
        
        /* Residual */
        // Diffusion.
        add_vector_form(new DefaultResidualDiffusion<double>(0, regions[i], new Hermes1DFunction<double>(D_map[i])));
        // Absorption.
        add_vector_form(new DefaultResidualVol<double>(0, regions[i], new Hermes2DFunction<double>(Sigma_a_map[i])));
        // Sources.
        add_vector_form(new DefaultVectorFormVol<double>(0, regions[i], new Hermes2DFunction<double>(-Q_map[i])));
      }
    }
  }
  
  namespace Common { namespace WeakForms
  {
    void NeutronicsProblem::get_source_part(WeakForm<double>* wf)
    {
      Hermes::vector<VectorFormVol<double> *> vfv = this->get_vfvol();
      Hermes::vector<VectorFormDG<double> *> vfDG = this->get_vfDG();
      Hermes::vector<VectorFormSurf<double> *> vfsurf = this->get_vfsurf();
      
      Hermes::vector<VectorFormVol<double> *>::const_iterator it1 = vfv.begin();
      for (; it1 != vfv.end(); ++it1)
        wf->add_vector_form((*it1)->clone());
      
      Hermes::vector<VectorFormDG<double> *>::const_iterator it2 = vfDG.begin();
      for (; it2 != vfDG.end(); ++it2)
        wf->add_vector_form_DG((*it2)->clone());
      
      Hermes::vector<VectorFormSurf<double> *>::const_iterator it3 = vfsurf.begin();
      for (; it3 != vfsurf.end(); ++it3)
        wf->add_vector_form_surf((*it3)->clone());
    }
  }
  }
   
  namespace Diffusion { namespace WeakForms
  {   
    void DiffusionWeakForm::add_forms_nonlinear()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool2 Ss_nnz = mp->get_scattering_nonzero_structure();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        bool do_include_fission = (include_fission != NONE);
        
        const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
        if (do_include_fission && !fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            do_include_fission = false;
                  
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 D = mp->get_D(*material);
        rank1 Sigma_r = mp->get_Sigma_r(*material);
        rank2 Sigma_s = mp->get_Sigma_s(*material);
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
        {
          add_matrix_form(new DiffusionReaction::Jacobian(regions, gto, D[gto], Sigma_r[gto], geom_type));
          add_vector_form(new DiffusionReaction::Residual(regions, gto, D[gto], Sigma_r[gto], geom_type));
          
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
          {
            if (Ss_nnz[gto][gfrom] && gto != gfrom)
            {
              add_matrix_form(new Scattering::Jacobian(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
              add_vector_form(new Scattering::Residual(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
            }
            
            if (do_include_fission && chi_nnz[gto])
            {
              if (include_fission == IMPLICIT)
                add_matrix_form( new FissionYield::Jacobian(regions, gto, gfrom, 
                                                                chi[gto], nu[gfrom], Sigma_f[gfrom], 
                                                                geom_type) );
              add_vector_form( new FissionYield::Residual(regions, gto, gfrom, 
                                                              chi[gto], nu[gfrom], Sigma_f[gfrom],
                                                              geom_type) );
            }
          }
        }
      }
    }
    
    void DiffusionWeakForm::add_forms()
    {    
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool2 Ss_nnz = mp->get_scattering_nonzero_structure();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        bool do_include_fission = (include_fission != NONE);
        
        const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
        if (do_include_fission && !fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            do_include_fission = false;
                  
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 D = mp->get_D(*material);
        rank1 Sigma_r = mp->get_Sigma_r(*material);
        rank2 Sigma_s = mp->get_Sigma_s(*material);
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
        {
          add_matrix_form(new DiffusionReaction::Jacobian(regions, gto, D[gto], Sigma_r[gto], geom_type));
          
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
          {
            if (Ss_nnz[gto][gfrom] && gto != gfrom)
              add_matrix_form(new Scattering::Jacobian(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
                       
            if (do_include_fission && chi_nnz[gto])
            {
              if (include_fission == IMPLICIT)
                add_matrix_form( new FissionYield::Jacobian(regions, gto, gfrom, 
                                                                chi[gto], nu[gfrom], Sigma_f[gfrom], 
                                                                geom_type) );              
              else if (include_fission == EXPLICIT)
                add_vector_form( new FissionYield::Residual(regions, gto, gfrom, 
                                                                chi[gto], -nu[gfrom], Sigma_f[gfrom],
                                                                geom_type) );
            }
          }
        }
      }
    }
    
    void DiffusionWeakForm::add_fission_sparse_structure()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
        if (!fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            continue;
                  
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
            if (chi_nnz[gto])
              add_matrix_form( new FissionYield::Jacobian(regions, gto, gfrom, 0.0, 0.0, 0.0, geom_type) );              
      }
    }
    
    void DiffusionWeakForm::get_fission_yield_part(WeakForm<double>* wf)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {        
        const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
        if (!fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            continue;
                  
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
            if (chi_nnz[gto])
              wf->add_matrix_form( new FissionYield::Jacobian(regions, gto, gfrom, 
                                                              chi[gto], -nu[gfrom], Sigma_f[gfrom], 
                                                              geom_type) );
      }
    }
    
    void DiffusionWeakForm::get_diffusion_reaction_part(WeakForm<double>* wf, 
                                                        const Hermes::vector<std::string>& vacuum_boundaries)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {          
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 D = mp->get_D(*material);
        rank1 Sigma_r = mp->get_Sigma_r(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
          wf->add_matrix_form(new DiffusionReaction::Jacobian(regions, gto, D[gto], Sigma_r[gto], geom_type));
      }
      
      if (!vacuum_boundaries.empty())
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
          wf->add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(vacuum_boundaries, gto));
    }
    
    void DiffusionWeakForm::get_scattering_part(WeakForm<double>* wf)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool2 Ss_nnz = mp->get_scattering_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {           
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank2 Sigma_s = mp->get_Sigma_s(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
            if (Ss_nnz[gto][gfrom] && gto != gfrom)
              wf->add_matrix_form(new Scattering::Jacobian(regions, gto, gfrom, -Sigma_s[gto][gfrom], geom_type));
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           GeomType geom_type, bool solve_by_newton) 
      : DiffusionWeakForm(&matprop, geom_type)
    { 
      bool explicit_sources = !matprop.get_iso_src().empty();
      
      if (explicit_sources)
      {
        include_fission = IMPLICIT;
        add_forms_from_homogeneous_part(solve_by_newton);
      }
      else
      {
        include_fission = EXPLICIT;
        add_forms_from_homogeneous_part(solve_by_newton);
      }
      
      if (explicit_sources)
      {
        std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
        for ( ; material != matprop.get_materials_list().end(); ++material)
        {
          Hermes::vector<std::string> regions = matprop.get_regions(*material);
          
          rank1 src_data = matprop.get_iso_src(*material);
          for (unsigned int gto = 0; gto < G; gto++)
            add_vector_form(new ExternalSources::LinearForm(regions, gto, solve_by_newton ? -src_data[gto] : src_data[gto], geom_type));
        }
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           Hermes2DFunction<double> *src, const std::string& src_area,
                                           GeomType geom_type, bool solve_by_newton  ) 
      : DiffusionWeakForm(&matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           Hermes2DFunction<double> *src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type, bool solve_by_newton  ) 
      : DiffusionWeakForm(&matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<Hermes2DFunction<double>*>& src,
                                           const std::string& src_area, 
                                           GeomType geom_type, bool solve_by_newton ) 
      : DiffusionWeakForm(&matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      if (src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, src[gto], geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<Hermes2DFunction<double>*>& src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type, bool solve_by_newton ) 
      : DiffusionWeakForm(&matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      if (src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, src[gto], geom_type));
    }
    
  /* WeakForms */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakForms
  {
    void SPNWeakForm::determine_symmetries()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
            
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
      
      sym.assign(N_odd * G, bool1(N_odd * G, false));
      present.assign(N_odd * G, bool1(N_odd * G, false));
      
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
              
              if ( (j-i)%G ) // inter-group coupling
              { 
                for (unsigned int k = 0; k <= 2*m; k+=2)
                {
                  if (!diagonal_moments[k])
                  {
                    // if any of the even Sigma_r moments appearing in the (m,n) equation is
                    // non-symmetric due to inter-group coupling...
                    present[j][i] = present[i][j] = true;
                    break;
                  }
                }
              }
              else  // intra-group
              {
                present[j][i] = true;
                sym[j][i] = sym[i][j] = true;
              }
            }
          }
        }
      }
    }
    
    void SPNWeakForm::add_forms_nonlinear()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
            
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();        
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
/* DEBUG      
      std::cout << std::endl;
      for (unsigned int gto = 0; gto < G; gto++)
      {
        for (unsigned int m = 0; m < N_odd; m++)
        {
          unsigned int i = mg.pos(m, gto);
          
          for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          {
            for (unsigned int n = 0; n < N_odd; n++)
            {
              unsigned int j = mg.pos(n, gfrom);
              
              std::cout << "(" << m << "," << n << " ; " << gto << "," << gfrom << ") --- i = " << i << ", j = " << j;
              std::cout << "\t p" << present[i][j] << " s" << sym[i][j] << std::endl;
            }
          }
        }
      }
      std::cout << std::endl;
*/    
      //int dssrJ = 0, dssrR = 0, fyJ = 0, fyR = 0, odsJ = 0, odsR = 0, odrJ = 0, odrR = 0; // DEBUG
      
      
      const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        bool do_include_fission = (include_fission != NONE);
        
        if (do_include_fission && !fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            do_include_fission = false;
          
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank3 Sigma_rn = mp->get_Sigma_rn(*material);
        rank3 odd_Sigma_rn_inv = mp->get_odd_Sigma_rn_inv(*material);
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < G; gto++)
        {
          for (unsigned int m = 0; m < N_odd; m++)
          {
            unsigned int i = mg.pos(m, gto);
            
            double Sigma_r = 0.;
            for (unsigned int k = 0; k <= m; k++)
              Sigma_r += Coeffs::system_matrix(m, m, k) * Sigma_rn[2*k][gto][gto];

            double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gto];
      
            add_matrix_form(new DiagonalStreamingAndReactions::Jacobian(regions, m, gto, G, D, Sigma_r, geom_type));// dssrJ++;
            add_vector_form(new DiagonalStreamingAndReactions::Residual(regions, m, gto, G, D, Sigma_r, geom_type));// dssrR++;
            
            if (do_include_fission && chi_nnz[gto]) {
              add_vector_form(new FissionYield::Residual(regions, m, N_odd, gto, G, chi[gto], nu, Sigma_f, geom_type));// fyR++;
            }
            
            add_vector_form(new OffDiagonalReactions::Residual(regions, m, N_odd, gto, G, Sigma_rn, geom_type));// odrR++;
            
            if (G > 1) {
              add_vector_form(new OffDiagonalStreaming::Residual(regions, m, gto, G, 
                                                                     odd_Sigma_rn_inv[m][gto], geom_type));// odsR++;
            }
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              if (gfrom != gto) {
                double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gfrom];
                add_matrix_form(new OffDiagonalStreaming::Jacobian(regions, m, gto, gfrom, G, D, geom_type));// odsJ++;
              }
              
              for (unsigned int n = 0; n < N_odd; n++)
              {
                unsigned int j = mg.pos(n, gfrom);
                
                if (do_include_fission && chi_nnz[gto] && include_fission == IMPLICIT) {
                  add_matrix_form( new FissionYield::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                  chi[gto], nu[gfrom], Sigma_f[gfrom], geom_type) );// fyJ++;
                }
                
                //// cout << "(" << i << "," << j << ") : P" << present[i][j] << " S" << sym[i][j] << endl;
                
                if (i != j)
                {
                  if (present[i][j]) {
                    double Sigma_rn_local = 0.;
                    for (unsigned int k = 0; k < m; k++)
                      Sigma_rn_local += Coeffs::system_matrix(m, n, k) * Sigma_rn[2*k][gto][gfrom];
      
                    add_matrix_form( new OffDiagonalReactions::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                            Sigma_rn_local, geom_type, 
                                                                            sym[i][j] ? HERMES_SYM : HERMES_NONSYM) );// odrJ++;
                  }
                }
              }
            }
          }
        }
      }
/* DEBUG       
      std::cout << "DiagonalStreamingAndReactions::Jacobian: " << dssrJ << std::endl;
      std::cout << "DiagonalStreamingAndReactions::Residual: " << dssrR << std::endl;
      std::cout << "FissionYield::Jacobian: "  << fyJ << std::endl;
      std::cout << "FissionYield::Residual: " << fyR << std::endl;
      std::cout << "OffDiagonalStreaming::Jacobian: " << odsJ << std::endl;
      std::cout << "OffDiagonalStreaming::Residual: " << odsR << std::endl;
      std::cout << "OffDiagonalReactions::Jacobian: " << odrJ << std::endl;
      std::cout << "OffDiagonalReactions::Residual: " << odrR << std::endl;
*/      
    }
    
    void SPNWeakForm::add_forms()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
/* DEBUG     
      std::cout << std::endl;
      for (unsigned int gto = 0; gto < G; gto++)
      {
        for (unsigned int m = 0; m < N_odd; m++)
        {
          unsigned int i = mg.pos(m, gto);
          
          for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          {
            for (unsigned int n = 0; n < N_odd; n++)
            {
              unsigned int j = mg.pos(n, gfrom);
              
              std::cout << "(" << m << "," << n << " ; " << gto << "," << gfrom << ") --- i = " << i << ", j = " << j;
              std::cout << "\t p" << present[i][j] << " s" << sym[i][j] << std::endl;
            }
          }
        }
      }
      std::cout << std::endl;

      int dssrJ = 0, dssrR = 0, fyJ = 0, fyR = 0, odsJ = 0, odsR = 0, odrJ = 0, odrR = 0; // DEBUG
*/      
      
      const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        // cout << endl << "          " << *material << endl;
        bool do_include_fission = (include_fission != NONE);
        
        if (do_include_fission && !fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            do_include_fission = false;
          
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank3 Sigma_rn = mp->get_Sigma_rn(*material);
        rank3 odd_Sigma_rn_inv = mp->get_odd_Sigma_rn_inv(*material);
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < G; gto++)
        {
          // cout << "gto=" << gto << endl;
          for (unsigned int m = 0; m < N_odd; m++)
          {
            unsigned int i = mg.pos(m, gto);
            
            double Sigma_r = 0.;
            for (unsigned int k = 0; k <= m; k++)
              Sigma_r += Coeffs::system_matrix(m, m, k) * Sigma_rn[2*k][gto][gto];

            double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gto];
            // cout << "  D: " << D << endl;
            // cout << "  Sr: " << Sigma_r << endl;
      
            add_matrix_form(new DiagonalStreamingAndReactions::Jacobian(regions, m, gto, G, D, Sigma_r, geom_type));// dssrJ++;
            
            if (do_include_fission && chi_nnz[gto] && include_fission == EXPLICIT) {
              // cout << "  nSf: " << nu[gto]*Sigma_f[gto] << endl;
              add_vector_form(new FissionYield::Residual(regions, m, N_odd, gto, G, -chi[gto], nu, Sigma_f, geom_type));// fyR++;
            }
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              // cout << "  gfrom=" << gfrom << endl;
              if (gfrom != gto && !diagonal_moments[2*m+1]) {
                double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gfrom];
                // cout << "    D: " << D << endl;
                add_matrix_form(new OffDiagonalStreaming::Jacobian(regions, m, gto, gfrom, G, D, geom_type));// odsJ++;
              }
              
              for (unsigned int n = 0; n < N_odd; n++)
              {
                unsigned int j = mg.pos(n, gfrom);
                
                if ((G == 1 && n >= m) || G > 1)
                  if (do_include_fission && chi_nnz[gto] && include_fission == IMPLICIT) {
                    add_matrix_form( new FissionYield::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                    chi[gto], nu[gfrom], Sigma_f[gfrom], geom_type) );// fyJ++;
                  }
                
                // cout << "(" << i << "," << j << ") : P" << present[i][j] << " S" << sym[i][j] << endl;
                
                if (i != j)
                {
                  if (present[i][j]) {
                    double Sigma_rn_local = 0.;
                    for (unsigned int k = 0; k <= m; k++)
                      Sigma_rn_local += Coeffs::system_matrix(m, n, k) * Sigma_rn[2*k][gto][gfrom];
      
                    // cout << "    Srnl: " << Sigma_rn_local << endl;
                    add_matrix_form( new OffDiagonalReactions::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                            Sigma_rn_local, geom_type, 
                                                                            sym[i][j] ? HERMES_SYM : HERMES_NONSYM) );// odrJ++;
                  }
                }
              }
            }
          }
        }
      }
/* DEBUG      
      std::cout << std::endl;
      std::cout << "DiagonalStreamingAndReactions::Jacobian: " << dssrJ << std::endl;
      std::cout << "DiagonalStreamingAndReactions::Residual: " << dssrR << std::endl;
      std::cout << "FissionYield::Jacobian: "  << fyJ << std::endl;
      std::cout << "FissionYield::Residual: " << fyR << std::endl;
      std::cout << "OffDiagonalStreaming::Jacobian: " << odsJ << std::endl;
      std::cout << "OffDiagonalStreaming::Residual: " << odsR << std::endl;
      std::cout << "OffDiagonalReactions::Jacobian: " << odrJ << std::endl;
      std::cout << "OffDiagonalReactions::Residual: " << odrR << std::endl;
*/       
    }
    
    void SPNWeakForm::add_fission_sparse_structure()
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {       
        if (!fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            continue;
          
        Hermes::vector<std::string> regions = mp->get_regions(*material);

        for (unsigned int gto = 0; gto < G; gto++)
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
              for (unsigned int n = (G==1) ? m : 0; n < N_odd; n++)
                if (chi_nnz[gto])
                  add_matrix_form( new FissionYield::Jacobian(regions, m, n, gto, gfrom, G, 0.0, 0.0, 0.0, geom_type) );
      }
    }
    
    void SPNWeakForm::get_fission_yield_part(WeakForm<double>* wf)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 chi_nnz = mp->get_fission_nonzero_structure();

      const Hermes::vector<std::string>& fission_materials = mp->get_fission_materials();
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {        
        if (!fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            continue;
          
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < G; gto++)
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
              for (unsigned int n = (G==1) ? m : 0; n < N_odd; n++)
                if (chi_nnz[gto])
                  wf->add_matrix_form( new FissionYield::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                  chi[gto], -nu[gfrom], Sigma_f[gfrom], geom_type) );
      }
    }
    
    void SPNWeakForm::get_symmetric_diffusion_reaction_part(WeakForm<double>* wf, 
                                                            const Hermes::vector<std::string>& vacuum_boundaries)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {    
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank3 Sigma_rn = mp->get_Sigma_rn(*material);
        rank3 odd_Sigma_rn_inv = mp->get_odd_Sigma_rn_inv(*material);
        
        for (unsigned int gto = 0; gto < G; gto++)
        {
          for (unsigned int m = 0; m < N_odd; m++)
          {
            unsigned int i = mg.pos(m, gto);
            
            double Sigma_r = 0.;
            for (unsigned int k = 0; k <= m; k++)
              Sigma_r += Coeffs::system_matrix(m, m, k) * Sigma_rn[2*k][gto][gto];

            double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gto];
                  
            wf->add_matrix_form(new DiagonalStreamingAndReactions::Jacobian(regions, m, gto, G, D, Sigma_r, geom_type));
                  
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {              
              for (unsigned int n = 0; n < N_odd; n++)
              {
                unsigned int j = mg.pos(n, gfrom);
                
                if (i != j)
                {
                  if (present[i][j] && sym[i][j]) 
                  {
                    double Sigma_rn_local = 0.;
                    for (unsigned int k = 0; k <= m; k++)
                      Sigma_rn_local += Coeffs::system_matrix(m, n, k) * Sigma_rn[2*k][gto][gfrom];
      
                    wf->add_matrix_form( new OffDiagonalReactions::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                            Sigma_rn_local, geom_type, 
                                                                            HERMES_SYM) );
                  }
                }
              }
            }
          }
        }
      }
      
      if (!vacuum_boundaries.empty())
        for (unsigned int m = 0; m < N_odd; m++)
          for (unsigned int n = 0; n < N_odd; n++)
            for (unsigned int gto = 0; gto < G; gto++)
              wf->add_matrix_form_surf(new WeakFormParts::VacuumBoundaryCondition::Jacobian(vacuum_boundaries, m, n, gto, G));
    }
    
    void SPNWeakForm::get_nonsymmetric_diffusion_reaction_part(WeakForm<double>* wf)
    {
      const MaterialPropertyMaps* mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {    
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank3 Sigma_rn = mp->get_Sigma_rn(*material);
        rank3 odd_Sigma_rn_inv = mp->get_odd_Sigma_rn_inv(*material);
        
        for (unsigned int gto = 0; gto < G; gto++)
        {
          for (unsigned int m = 0; m < N_odd; m++)
          {
            unsigned int i = mg.pos(m, gto);
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              if (gfrom != gto && !diagonal_moments[2*m+1]) {
                double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gfrom];
                wf->add_matrix_form(new OffDiagonalStreaming::Jacobian(regions, m, gto, gfrom, G, D, geom_type));
              }
              
              for (unsigned int n = 0; n < N_odd; n++)
              {
                unsigned int j = mg.pos(n, gfrom);
                
                if (i != j)
                {
                  if (present[i][j] && !sym[i][j]) 
                  {
                    double Sigma_rn_local = 0.;
                    for (unsigned int k = 0; k <= m; k++)
                      Sigma_rn_local += Coeffs::system_matrix(m, n, k) * Sigma_rn[2*k][gto][gfrom];
      
                    wf->add_matrix_form( new OffDiagonalReactions::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                            Sigma_rn_local, geom_type, 
                                                                            HERMES_NONSYM) );
                  }
                }
              }
            }
          }
        }
      }
    }
      
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                           GeomType geom_type, bool solve_by_newton) 
      : SPNWeakForm(N, &matprop, geom_type)
    { 
      bool explicit_sources = !matprop.get_iso_src().empty();
      
      if (explicit_sources)
      {
        include_fission = IMPLICIT;
        add_forms_from_homogeneous_part(solve_by_newton);
        
        std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
        for ( ; material != matprop.get_materials_list().end(); ++material)
        {
          Hermes::vector<std::string> regions = matprop.get_regions(*material);
          
          rank1 src_data = matprop.get_iso_src(*material);
          
          for (unsigned int m = 0; m < N_odd; m++)
            for (unsigned int gto = 0; gto < G; gto++)
              add_vector_form(new ExternalSources::LinearForm(regions, m, gto, G, solve_by_newton ? -src_data[gto] : src_data[gto], geom_type));
        }
      }
      else
      {
        include_fission = EXPLICIT;
        add_forms_from_homogeneous_part(solve_by_newton);
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           Hermes2DFunction<double> *iso_src, std::string src_area,
                                           GeomType geom_type, bool solve_by_newton  )
      : SPNWeakForm(N, &matprop, geom_type)
    { 
      add_forms_from_homogeneous_part(solve_by_newton);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, iso_src, geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           Hermes2DFunction<double> *iso_src,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type, bool solve_by_newton  )
      : SPNWeakForm(N, &matprop, geom_type)
    {  
      add_forms_from_homogeneous_part(solve_by_newton);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, iso_src, geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<Hermes2DFunction<double>*>& iso_src,
                                           std::string src_area, 
                                           GeomType geom_type, bool solve_by_newton )
      : SPNWeakForm(N, &matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      if (iso_src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
            
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, iso_src[gto], geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }                                                                                   
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<Hermes2DFunction<double>*>& iso_src,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type, bool solve_by_newton )
      : SPNWeakForm(N, &matprop, geom_type)
    {
      add_forms_from_homogeneous_part(solve_by_newton);
      
      if (iso_src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, iso_src[gto], geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }            
  /* WeakForms */
  }
  /* Diffusion */
  }
  
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 