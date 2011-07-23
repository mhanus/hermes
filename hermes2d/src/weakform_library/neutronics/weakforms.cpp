#include "neutronics/weakforms.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace WeakForms 
{
  namespace Simple
  {    
    FixedSourceProblem::FixedSourceProblem(Hermes::vector<std::string> regions, 
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
  
  namespace Diffusion
  {   
    void FixedSourceProblem::lhs_init(unsigned int G, 
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
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           GeomType geom_type) : WeakForm(matprop.get_G())
    {
      lhs_init(matprop.get_G(), matprop, geom_type);
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
        add_vector_form(new ExternalSources::LinearForm(gto, matprop, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           HermesFunction *minus_f_src, const std::string& src_area,
                                           GeomType geom_type  ) : WeakForm(matprop.get_G())
    {
      lhs_init(matprop.get_G(), matprop, geom_type);
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, minus_f_src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           HermesFunction *minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type  ) : WeakForm(matprop.get_G())
    {
      lhs_init(matprop.get_G(), matprop, geom_type);
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, minus_f_src, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<HermesFunction*>& minus_f_src,
                                           const std::string& src_area, 
                                           GeomType geom_type ) : WeakForm(matprop.get_G())
    {
      if (minus_f_src.size() != matprop.get_G())
        error(Messages::E_INVALID_SIZE);
      
      lhs_init(matprop.get_G(), matprop, geom_type);
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, minus_f_src[gto], geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<HermesFunction*>& minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type ) : WeakForm(matprop.get_G())
    {
      if (minus_f_src.size() != matprop.get_G())
        error(Messages::E_INVALID_SIZE);
      
      lhs_init(matprop.get_G(), matprop, geom_type);
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, minus_f_src[gto], geom_type));
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                 const Hermes::vector<MeshFunction*>& iterates,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type ) 
      : WeakForm(matprop.get_G()), Common::KeffEigenvalueProblem(initial_keff_guess)
    {      
      init(matprop, iterates, geom_type);
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                 const Hermes::vector<Solution*>& iterates,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type ) 
      : WeakForm(matprop.get_G()), Common::KeffEigenvalueProblem(initial_keff_guess)
    {      
      Hermes::vector<MeshFunction *> iterates_mf;
      for (unsigned int i = 0; i < iterates.size(); i++)
        iterates_mf.push_back(static_cast<MeshFunction*>(iterates[i]));
      
      init(matprop, iterates_mf, geom_type);
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                 const Hermes::vector<MeshFunction*>& iterates, 
                                                 const Hermes::vector<std::string>& src_areas,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type )
      : WeakForm(matprop.get_G()), Common::KeffEigenvalueProblem(initial_keff_guess)
    {
      init(matprop, iterates, geom_type, src_areas);
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                 const Hermes::vector<Solution*>& iterates, 
                                                 const Hermes::vector<std::string>& src_areas,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type ) 
      : WeakForm(matprop.get_G()), Common::KeffEigenvalueProblem(initial_keff_guess)
    {
      Hermes::vector<MeshFunction *> iterates_mf;
      for (unsigned int i = 0; i < iterates.size(); i++)
        iterates_mf.push_back(static_cast<MeshFunction*>(iterates[i]));
      
      init(matprop, iterates_mf, geom_type, src_areas);
    }
    
    void KeffEigenvalueProblem::init(const MaterialPropertyMaps& matprop,
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
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
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
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                           GeomType geom_type) 
      : WeakFormHomogeneous(N, matprop, geom_type, true)
    {
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
          add_vector_form(new ExternalSources::LinearForm(m, gto, matprop, geom_type));
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
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
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
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
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<HermesFunction*>& minus_isotropic_sources,
                                           std::string src_area, 
                                           GeomType geom_type )
      : WeakFormHomogeneous(N, matprop, geom_type, true)
    {
      if (minus_isotropic_sources.size() != G)
        error(Messages::E_INVALID_SIZE);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_area, minus_isotropic_sources[gto], geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<HermesFunction*>& minus_isotropic_sources,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type )
      : WeakFormHomogeneous(N, matprop, geom_type, true)
    {
      if (minus_isotropic_sources.size() != G)
        error(Messages::E_INVALID_SIZE);
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol *src = new WeakFormsH1::DefaultVectorFormVol(mg.pos(m,gto), src_areas, minus_isotropic_sources[gto], geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::~FixedSourceProblem()
    {
      std::vector<VectorFormVol*>::const_iterator it = source_terms.begin();
      for ( ; it != source_terms.end(); ++it)
        delete *it;
      source_terms.clear();
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                                 const Hermes::vector<Solution*>& iterates, 
                                                 const Hermes::vector<std::string>& src_areas,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type )
      : WeakFormHomogeneous(N, matprop, geom_type, false), Common::KeffEigenvalueProblem(initial_keff_guess)
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
    
    KeffEigenvalueProblem::~KeffEigenvalueProblem()
    {
      std::vector<FissionYield::OuterIterationForm*>::const_iterator it = keff_iteration_forms.begin();
      for ( ; it != keff_iteration_forms.end(); ++it)
        delete *it;
      keff_iteration_forms.clear();
      
      SupportClasses::SPN::MomentFilter::clear_scalar_fluxes(&scalar_flux_iterates);
    }
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
    { 
      keff = new_keff;
      
      std::vector<FissionYield::OuterIterationForm*>::iterator it = keff_iteration_forms.begin();
      for ( ; it != keff_iteration_forms.end(); ++it)
        (*it)->update_keff(new_keff); 
    }
  }
    
/* WeakForms */
}
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 