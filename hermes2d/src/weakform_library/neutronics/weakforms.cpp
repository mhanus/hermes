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
    
  namespace Diffusion { namespace WeakForms
  {   
    void DiffusionWeakForm::add_forms_from_homogeneous_part(Common::WeakForms::NeutronicsProblem* wf, 
                                                            const MaterialPropertyMaps *mp, GeomType geom_type, bool include_fission)
    {      
      bool2 Ss_nnz = mp->get_scattering_nonzero_structure();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
        Hermes::vector<std::string> regions = mp->get_regions(*material);
        
        rank1 D = mp->get_D(*material);
        rank1 Sigma_r = mp->get_Sigma_r(*material);
        rank2 Sigma_s = mp->get_Sigma_s(*material);
        rank1 Sigma_f = mp->get_Sigma_f(*material);
        rank1 chi = mp->get_chi(*material);
        rank1 nu = mp->get_nu(*material);
        
        for (unsigned int gto = 0; gto < mp->get_G(); gto++)
        {
          wf->add_matrix_form(new DiffusionReaction::Jacobian(regions, gto, D[gto], Sigma_r[gto], geom_type));
          wf->add_vector_form(new DiffusionReaction::Residual(regions, gto, D[gto], Sigma_r[gto], geom_type));
          
          for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
          {
            if (Ss_nnz[gto][gfrom] && gto != gfrom)
            {
              wf->add_matrix_form(new Scattering::Jacobian(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
              wf->add_vector_form(new Scattering::Residual(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
            }
            
            if (include_fission && chi_nnz[gto])
            {
              wf->add_matrix_form( new FissionYield::Jacobian(regions, gto, gfrom, 
                                                              chi[gto], nu[gfrom], Sigma_f[gfrom], 
                                                              geom_type) );
              wf->add_vector_form( new FissionYield::Residual(regions, gto, gfrom, 
                                                              chi[gto], nu[gfrom], Sigma_f[gfrom],
                                                              geom_type) );
            }
          }
        }
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           GeomType geom_type) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {            
      std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
      for ( ; material != matprop.get_materials_list().end(); ++material)
      {
        Hermes::vector<std::string> regions = matprop.get_regions(*material);
        
        rank1 src_data = matprop.get_iso_src(*material);
        for (unsigned int gto = 0; gto < G; gto++)
          add_vector_form(new ExternalSources::LinearForm(regions, gto, src_data[gto], geom_type));
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           Hermes2DFunction<double> *minus_f_src, const std::string& src_area,
                                           GeomType geom_type  ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, minus_f_src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           Hermes2DFunction<double> *minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type  ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, minus_f_src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                                           const std::string& src_area, 
                                           GeomType geom_type ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      if (minus_f_src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, minus_f_src[gto], geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<Hermes2DFunction<double>*>& minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      if (minus_f_src.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int gto = 0; gto < G; gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, minus_f_src[gto], geom_type));
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                const Hermes::vector<Solution<double>*>& iterates,
                                                double initial_keff_guess, 
                                                GeomType geom_type ) 
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess)
    {            
      init_rhs(iterates);
    }

    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                const Hermes::vector<Solution<double>*>& iterates, 
                                                const Hermes::vector<std::string>& fission_materials,
                                                double initial_keff_guess, 
                                                GeomType geom_type ) 
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess, fission_materials)
    {
      init_rhs(iterates);
    }
    
    void KeffEigenvalueProblem::init_rhs(const Hermes::vector<Solution<double>*>& iterates)
    {
      const Diffusion::MaterialProperties::MaterialPropertyMaps *mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      stored_flux_solutions.reserve(G);
      scalar_flux_iterates.reserve(G);
      for (unsigned int gto = 0; gto < G; gto++)
      { 
        stored_flux_solutions.push_back(iterates[gto]);
        scalar_flux_iterates.push_back(static_cast<MeshFunction<double>*>(stored_flux_solutions.back()));
      }
      
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
          add_vector_form(new FissionYield::OuterIterationForm( regions, gto, 
                                                                chi[gto], nu, Sigma_f,
                                                                scalar_flux_iterates, keff, 
                                                                geom_type ) );
      }
    }
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
    { 
      keff = new_keff;
      
      Hermes::vector<Form<double> *>::const_iterator it = get_forms().begin();
      for ( ; it != get_forms().end(); ++it)
      { 
        FissionYield::OuterIterationForm* keff_iteration_form = dynamic_cast<FissionYield::OuterIterationForm*>(*it);
        if (keff_iteration_form != NULL)
          keff_iteration_form->update_keff(new_keff); 
      }
    }
    
    void KeffEigenvalueProblem::update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed)
    {
      for (unsigned int gto = 0; gto < G; gto++)
        stored_flux_solutions[gto]->copy(new_solutions[gto]);
    }
    
  /* WeakForms */
  }
  /* Diffusion */
  }

  namespace SPN { namespace WeakForms
  {
    void SPNWeakForm::add_forms_from_homogeneous_part(Common::WeakForms::NeutronicsProblem* wf, 
                                                      const MaterialPropertyMaps *mp, GeomType geom_type, bool include_fission)
    {
      unsigned int G = mp->get_G();
      
      bool1 diagonal_moments = mp->is_Sigma_rn_diagonal();
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
                for (unsigned int k = 0; k <= 2*m; k+=2)
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
      
      std::set<std::string>::const_iterator material = mp->get_materials_list().begin();
      for ( ; material != mp->get_materials_list().end(); ++material)
      {
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
      
            wf->add_matrix_form(new DiagonalStreamingAndReactions::Jacobian(regions, m, gto, G, D, Sigma_r, geom_type));// dssrJ++;
            wf->add_vector_form(new DiagonalStreamingAndReactions::Residual(regions, m, gto, G, D, Sigma_r, geom_type));// dssrR++;
            
            if (include_fission && chi_nnz[gto]) {
              wf->add_vector_form(new FissionYield::Residual(regions, m, N_odd, gto, G, chi[gto], nu, Sigma_f, geom_type));// fyR++;
            }
            
            wf->add_vector_form(new OffDiagonalReactions::Residual(regions, m, N_odd, gto, G, Sigma_rn, geom_type));// odrR++;
            
            if (G > 1) {
              wf->add_vector_form(new OffDiagonalStreaming::Residual(regions, m, gto, G, 
                                                                     odd_Sigma_rn_inv[m][gto], geom_type));// odsR++;
            }
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              if (gfrom != gto) {
                double D = -Coeffs::D(m) * odd_Sigma_rn_inv[m][gto][gfrom];
                wf->add_matrix_form(new OffDiagonalStreaming::Jacobian(regions, m, gto, gfrom, G, D, geom_type));// odsJ++;
              }
              
              for (unsigned int n = 0; n < N_odd; n++)
              {
                unsigned int j = mg.pos(n, gfrom);
                
                if (include_fission && chi_nnz[gto]) {
                  wf->add_matrix_form( new FissionYield::Jacobian(regions, m, n, gto, gfrom, G, 
                                                                  chi[gto], nu[gfrom], Sigma_f[gfrom], geom_type) );// fyJ++;
                }
                
                //// cout << "(" << i << "," << j << ") : P" << present[i][j] << " S" << sym[i][j] << endl;
                
                if (i != j)
                {
                  if (present[i][j]) {
                    double Sigma_rn_local = 0.;
                    for (unsigned int k = 0; k < m; k++)
                      Sigma_rn_local += Coeffs::system_matrix(m, n, k) * Sigma_rn[2*k][gto][gfrom];
      
                    wf->add_matrix_form( new OffDiagonalReactions::Jacobian(regions, m, n, gto, gfrom, G, 
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
   
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                           GeomType geom_type) 
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {      
      std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
      for ( ; material != matprop.get_materials_list().end(); ++material)
      {
        Hermes::vector<std::string> regions = matprop.get_regions(*material);
        
        rank1 src_data = matprop.get_iso_src(*material);
        
        for (unsigned int m = 0; m < N_odd; m++)
          for (unsigned int gto = 0; gto < G; gto++)
            add_vector_form(new ExternalSources::LinearForm(regions, m, gto, G, src_data[gto], geom_type));
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           Hermes2DFunction<double> *minus_isotropic_source, std::string src_area,
                                           GeomType geom_type  )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, minus_isotropic_source, geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           Hermes2DFunction<double> *minus_isotropic_source,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type  )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {            
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, minus_isotropic_source, geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                                           std::string src_area, 
                                           GeomType geom_type )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {
      if (minus_isotropic_sources.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
            
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, minus_isotropic_sources[gto], geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }                                                                                   
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<Hermes2DFunction<double>*>& minus_isotropic_sources,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {
      if (minus_isotropic_sources.size() != G)
        ErrorHandling::error_function(Messages::E_INVALID_SIZE);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, minus_isotropic_sources[gto], geom_type);
          src->setScalingFactor(Coeffs::even_moment(0, m));
          add_vector_form(src);
        }
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                                 const Hermes::vector<Solution<double>*>& iterates, 
                                                 const Hermes::vector<std::string>& fission_materials,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type )
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type, 
                                                 initial_keff_guess, fission_materials),
        SPNWeakForm(N, matprop.get_G())
    { 
      stored_flux_solutions.reserve(iterates.size());
      for (Hermes::vector<Solution<double>*>::const_iterator it = iterates.begin(); it != iterates.end(); ++it)
        stored_flux_solutions.push_back(*it);
      
      SupportClasses::MomentFilter::get_scalar_fluxes_with_derivatives(stored_flux_solutions, &scalar_flux_iterates, G);
      
      std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
      for ( ; material != matprop.get_materials_list().end(); ++material)
      {
        if (!fission_materials.empty()) 
          if (std::find(fission_materials.begin(), fission_materials.end(), *material) == fission_materials.end())
            continue;
          
        Hermes::vector<std::string> regions = matprop.get_regions(*material);
        
        rank1 Sigma_f = matprop.get_Sigma_f(*material);
        rank1 chi = matprop.get_chi(*material);
        rank1 nu = matprop.get_nu(*material);
        
        for (unsigned int m = 0; m < N_odd; m++)
          for (unsigned int gto = 0; gto < G; gto++)
            add_vector_form(new FissionYield::OuterIterationForm( regions, 
                                                                  m, gto, G, 
                                                                  chi[gto], nu, Sigma_f, 
                                                                  scalar_flux_iterates, initial_keff_guess, 
                                                                  geom_type ));
      }
    }
    
    KeffEigenvalueProblem::~KeffEigenvalueProblem()
    {  
      SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_flux_iterates);
    }
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
    { 
      keff = new_keff;
      
      Hermes::vector<Form<double> *>::const_iterator it = get_forms().begin();
      for ( ; it != get_forms().end(); ++it)
      { 
        FissionYield::OuterIterationForm* keff_iteration_form = dynamic_cast<FissionYield::OuterIterationForm*>(*it);
        if (keff_iteration_form != NULL)
          keff_iteration_form->update_keff(new_keff); 
      }
    }
    
    void KeffEigenvalueProblem::update_fluxes(const Hermes::vector<Solution<double>*>& new_solutions, bool meshes_changed)
    {
      Hermes::vector<Solution<double>*>::const_iterator new_solution = new_solutions.begin();
      Hermes::vector<Solution<double>*>::const_iterator stored_flux_solution = stored_flux_solutions.begin();
      
      for ( ; new_solution != new_solutions.end(); ++new_solution, ++stored_flux_solution)
      {
        // FIXME: It seems that Space::construct_refined_spaces prevents automatic determination of meshes_changed.
        //
        //if ((*new_solution)->get_mesh()->get_seq() != (*stored_flux_solution)->get_mesh()->get_seq())
        //  meshes_changed = true;
        
        (*stored_flux_solution)->copy(*new_solution);
      }
      
      if (meshes_changed)
      {
        Hermes::vector<MeshFunction<double>*>::const_iterator scalar_flux_iterate = scalar_flux_iterates.begin();
        for ( ; scalar_flux_iterate != scalar_flux_iterates.end(); ++scalar_flux_iterate)
          (*scalar_flux_iterate)->reinit();
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