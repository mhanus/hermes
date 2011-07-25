#include "neutronics/weakforms.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace WeakForms 
{
  namespace Simple
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
        add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, regions[i], new HermesFunction<double>(D_map[i]), 
                                                      HERMES_SYM));
        // Absorption.
        add_matrix_form(new DefaultMatrixFormVol<double>(0, 0, regions[i], new HermesFunction<double>(Sigma_a_map[i]), 
                                                  HERMES_SYM));
        
        /* Residual */
        // Diffusion.
        add_vector_form(new DefaultResidualDiffusion<double>(0, regions[i], new HermesFunction<double>(D_map[i])));
        // Absorption.
        add_vector_form(new DefaultResidualVol<double>(0, regions[i], new HermesFunction<double>(Sigma_a_map[i])));
        // Sources.
        add_vector_form(new DefaultVectorFormVol<double>(0, regions[i], new HermesFunction<double>(-Q_map[i])));
      }
    }
  }
  
  namespace Common
  {    
    HomogeneousPart::~HomogeneousPart()
    {
      std::vector<MatrixFormVol<double>*>::const_iterator matrices_it = matrix_forms.begin();
      for ( ; matrices_it != matrix_forms.end(); ++matrices_it)
        delete *matrices_it;
      matrix_forms.clear();
      
      std::vector<VectorFormVol<double>*>::const_iterator vectors_it = vector_forms.begin();
      for ( ; vectors_it != vector_forms.end(); ++vectors_it)
        delete *vectors_it;
      vector_forms.clear();
    }
    
    void NeutronicsProblem::add_forms_from_homogeneous_part()
    {
      std::vector<MatrixFormVol<double>*>::const_iterator matrices_it = homogeneous_part->matrix_forms.begin();
      for ( ; matrices_it != homogeneous_part->matrix_forms.end(); ++matrices_it)
        add_matrix_form(*matrices_it);
      
      std::vector<VectorFormVol<double>*>::const_iterator vectors_it = homogeneous_part->vector_forms.begin();
      for ( ; vectors_it != homogeneous_part->vector_forms.end(); ++vectors_it)
        add_vector_form(*vectors_it);
    }
    
    KeffEigenvalueProblem::~KeffEigenvalueProblem()
    {
      std::vector<SupportClasses::Common::SourceFilter*>::const_iterator it = source_filters.begin();
      for ( ; it != source_filters.end(); ++it)
        delete *it;
      source_filters.clear();
    }
  }
  
  namespace Diffusion
  {   
    HomogeneousPart::HomogeneousPart(const MaterialProperties::Common::MaterialPropertyMaps* matprop,
                                     GeomType geom_type, bool include_fission)
    {
      const MaterialPropertyMaps *mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      bool2 Ss_nnz = mp->get_scattering_nonzero_structure();
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      for (unsigned int gto = 0; gto < mp->get_G(); gto++)
      {
        matrix_forms.push_back(new DiffusionReaction::Jacobian(gto, *mp, geom_type));
        vector_forms.push_back(new DiffusionReaction::Residual(gto, *mp, geom_type));
        
        for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
        {
          if (Ss_nnz[gto][gfrom] && gto != gfrom)
          {
            matrix_forms.push_back(new Scattering::Jacobian(gto, gfrom, *mp, geom_type));
            vector_forms.push_back(new Scattering::Residual(gto, gfrom, *mp, geom_type));
          }
          
          if (include_fission && chi_nnz[gto])
          {
            matrix_forms.push_back(new FissionYield::Jacobian(gto, gfrom, *mp, geom_type));
            vector_forms.push_back(new FissionYield::Residual(gto, gfrom, *mp, geom_type));
          }
        }
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           GeomType geom_type) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        VectorFormVol<double> *src = new ExternalSources::LinearForm(gto, matprop, geom_type);
        source_terms.push_back(src);
        add_vector_form(src);
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           HermesFunction<double> *minus_f_src, const std::string& src_area,
                                           GeomType geom_type  ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, minus_f_src, geom_type);
        source_terms.push_back(src);
        add_vector_form(src);
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           HermesFunction<double> *minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type  ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, minus_f_src, geom_type);
        source_terms.push_back(src);
        add_vector_form(src);
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<HermesFunction<double>*>& minus_f_src,
                                           const std::string& src_area, 
                                           GeomType geom_type ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      if (minus_f_src.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_area, minus_f_src[gto], geom_type);
        source_terms.push_back(src);
        add_vector_form(src);
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, 
                                           const std::vector<HermesFunction<double>*>& minus_f_src,
                                           const Hermes::vector<std::string>& src_areas,
                                           GeomType geom_type ) 
      : NeutronicsProblem(matprop.get_G(), &matprop, geom_type)
    {
      if (minus_f_src.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(gto, src_areas, minus_f_src[gto], geom_type);
        source_terms.push_back(src);
        add_vector_form(src);
      }
    }
    
    FixedSourceProblem::~FixedSourceProblem()
    {
      std::vector<VectorFormVol<double>*>::const_iterator it = source_terms.begin();
      for ( ; it != source_terms.end(); ++it)
        delete *it;
      source_terms.clear();
    }
    
   KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                             const Hermes::vector<MeshFunction<double>*>& iterates,
                                             double initial_keff_guess, 
                                             GeomType geom_type ) 
      : Common::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      init_rhs(iterates);
    }

    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                const Hermes::vector<Solution<double>*>& iterates,
                                                double initial_keff_guess, 
                                                GeomType geom_type ) 
      : Common::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess)
    {      
      Hermes::vector<MeshFunction<double> *> iterates_mf;
      for (unsigned int i = 0; i < iterates.size(); i++)
        iterates_mf.push_back(static_cast<MeshFunction<double>*>(iterates[i]));
      
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      init_rhs(iterates_mf);
    }

    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                const Hermes::vector<MeshFunction<double>*>& iterates, 
                                                const Hermes::vector<std::string>& fission_regions,
                                                double initial_keff_guess, 
                                                GeomType geom_type )
      : Common::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess, fission_regions)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      init_rhs(iterates);
    }

    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop,
                                                const Hermes::vector<Solution<double>*>& iterates, 
                                                const Hermes::vector<std::string>& fission_regions,
                                                double initial_keff_guess, 
                                                GeomType geom_type ) 
      : Common::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess, fission_regions)
    {
      Hermes::vector<MeshFunction<double> *> iterates_mf;
      for (unsigned int i = 0; i < iterates.size(); i++)
        iterates_mf.push_back(static_cast<MeshFunction<double>*>(iterates[i]));
      
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      init_rhs(iterates_mf);
    }
    
    void KeffEigenvalueProblem::init_rhs(const Hermes::vector<MeshFunction<double>*>& iterates)
    {
      const MaterialPropertyMaps *mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      for (unsigned int gto = 0; gto < G; gto++)
      { 
        FissionYield::OuterIterationForm* keff_iteration_form;
        
        if (fission_regions.size() > 0) 
          keff_iteration_form = new FissionYield::OuterIterationForm( gto, fission_regions, *mp, iterates, keff, geom_type );
        else
          keff_iteration_form = new FissionYield::OuterIterationForm( gto, *mp, iterates, keff, geom_type );
        
        keff_iteration_forms.push_back(keff_iteration_form);
        add_vector_form(keff_iteration_form);
      }
    }
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
    { 
      keff = new_keff;
      
      std::vector<FissionYield::OuterIterationForm*>::const_iterator it = keff_iteration_forms.begin();
      for ( ; it != keff_iteration_forms.end(); ++it)
        (*it)->update_keff(new_keff); 
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter()
    {
      source_filters.push_back(new SupportClasses::Common::SourceFilter(*matprop, fission_regions, geom_type));
      return source_filters.back();
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter(const Hermes::vector<Solution<double>*>& solutions)
    {
      source_filters.push_back(new SupportClasses::Common::SourceFilter(solutions, *matprop, fission_regions, geom_type));
      return source_filters.back();
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions)
    {
      source_filters.push_back(new SupportClasses::Common::SourceFilter(solutions, *matprop, fission_regions, geom_type));
      return source_filters.back();
    }
  }

  namespace SPN
  {
    HomogeneousPart::HomogeneousPart(unsigned int N_odd,
                                     const MaterialProperties::Common::MaterialPropertyMaps* matprop,
                                     GeomType geom_type, bool include_fission)
    {
      unsigned int G = matprop->get_G();
      
      const MaterialPropertyMaps *mp = static_cast<const MaterialPropertyMaps*>(matprop);
      MomentGroupFlattener mg(G);
      
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
      
      bool1 chi_nnz = mp->get_fission_nonzero_structure();
      
      for (unsigned int gto = 0; gto < G; gto++)
      {
        for (unsigned int m = 0; m < N_odd; m++)
        {
          unsigned int i = mg.pos(m, gto);
          
          matrix_forms.push_back(new DiagonalStreamingAndReactions::Jacobian(m, gto, *mp, geom_type));
          vector_forms.push_back(new DiagonalStreamingAndReactions::Residual(m, gto, *mp, geom_type));
          
          if (include_fission && chi_nnz[gto])
            vector_forms.push_back(new FissionYield::Residual(m, N_odd, gto, *mp, geom_type));
          
          vector_forms.push_back(new OffDiagonalReactions::Residual(m, N_odd, gto, *mp, geom_type));
          
          if (G > 1)
            vector_forms.push_back(new OffDiagonalStreaming::Residual(m, gto, *mp, geom_type));
          
          for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          {
            if (gfrom != gto)
              matrix_forms.push_back(new OffDiagonalStreaming::Jacobian(m, gto, gfrom, *mp, geom_type));
            
            for (unsigned int n = 0; n < N_odd; n++)
            {
              unsigned int j = mg.pos(n, gfrom);
              
              if (include_fission && chi_nnz[gto])
                matrix_forms.push_back( new FissionYield::Jacobian(m, n, gto, gfrom, *mp, geom_type) );
              
              //// cout << "(" << i << "," << j << ") : P" << present[i][j] << " S" << sym[i][j] << endl;
              
              if (i != j)
              {
                if (present[i][j])
                  matrix_forms.push_back( new OffDiagonalReactions::Jacobian(m, n, gto, gfrom, *mp, geom_type, 
                                                                      sym[i][j] ? HERMES_SYM : HERMES_NONSYM) );
              }
            }
          }
        }
      }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                           GeomType geom_type) 
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type)
    {
      unsigned int N_odd = (N+1)/2;
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new ExternalSources::LinearForm(m, gto, matprop, geom_type);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           HermesFunction<double> *minus_isotropic_source, std::string src_area,
                                           GeomType geom_type  )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type)
    {
      unsigned int N_odd = (N+1)/2;
      MomentGroupFlattener mg(G);
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, minus_isotropic_source, geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           HermesFunction<double> *minus_isotropic_source,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type  )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type)
    {
      unsigned int N_odd = (N+1)/2;
      MomentGroupFlattener mg(G);
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, minus_isotropic_source, geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<HermesFunction<double>*>& minus_isotropic_sources,
                                           std::string src_area, 
                                           GeomType geom_type )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type)
    {
      if (minus_isotropic_sources.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
      unsigned int N_odd = (N+1)/2;
      MomentGroupFlattener mg(G);
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_area, minus_isotropic_sources[gto], geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }                                                                                   
    
    FixedSourceProblem::FixedSourceProblem(const MaterialPropertyMaps& matprop, unsigned int N, 
                                           const std::vector<HermesFunction<double>*>& minus_isotropic_sources,
                                           Hermes::vector<std::string> src_areas,
                                           GeomType geom_type )
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type)
    {
      if (minus_isotropic_sources.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
      unsigned int N_odd = (N+1)/2;
      MomentGroupFlattener mg(G);
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, true);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
        for (unsigned int gto = 0; gto < G; gto++)
        {
          VectorFormVol<double> *src = new WeakFormsH1::DefaultVectorFormVol<double>(mg.pos(m,gto), src_areas, minus_isotropic_sources[gto], geom_type);
          src->scaling_factor = Coeffs::even_moment(0, m);
          source_terms.push_back(src);
          add_vector_form(src);
        }
    }
    
    FixedSourceProblem::~FixedSourceProblem()
    {
      std::vector<VectorFormVol<double>*>::const_iterator it = source_terms.begin();
      for ( ; it != source_terms.end(); ++it)
        delete *it;
      source_terms.clear();
    }
    
    KeffEigenvalueProblem::KeffEigenvalueProblem(const MaterialPropertyMaps& matprop, unsigned int N,
                                                 const Hermes::vector<Solution<double>*>& iterates, 
                                                 const Hermes::vector<std::string>& fission_regions,
                                                 double initial_keff_guess, 
                                                 GeomType geom_type )
      : Common::KeffEigenvalueProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type, initial_keff_guess, fission_regions)
    {      
      SupportClasses::SPN::MomentFilter::get_scalar_fluxes_with_derivatives(iterates, &scalar_flux_iterates, G);
      
      unsigned int N_odd = (N+1)/2;
      MomentGroupFlattener mg(G);
      
      homogeneous_part = new HomogeneousPart(N_odd, &matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      for (unsigned int m = 0; m < N_odd; m++)
      {
        for (unsigned int gto = 0; gto < G; gto++)
        {            
          FissionYield::OuterIterationForm* keff_iteration_form;
          
          if (fission_regions.size() > 0) 
            keff_iteration_form = new FissionYield::OuterIterationForm( m, gto, fission_regions, matprop, scalar_flux_iterates, 
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
      
      std::vector<FissionYield::OuterIterationForm*>::const_iterator it = keff_iteration_forms.begin();
      for ( ; it != keff_iteration_forms.end(); ++it)
        (*it)->update_keff(new_keff); 
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter()
    {
      source_filters.push_back(new SupportClasses::SPN::SourceFilter(*matprop, fission_regions, geom_type));
      return source_filters.back();
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter(const Hermes::vector<Solution<double>*>& solutions)
    {
      source_filters.push_back(new SupportClasses::SPN::SourceFilter(solutions, *matprop, fission_regions, geom_type));
      return source_filters.back();
    }
    
    SupportClasses::Common::SourceFilter* KeffEigenvalueProblem::get_new_source_filter(const Hermes::vector<MeshFunction<double>*>& solutions)
    {
      source_filters.push_back(new SupportClasses::SPN::SourceFilter(solutions, *matprop, fission_regions, geom_type));
      return source_filters.back();
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