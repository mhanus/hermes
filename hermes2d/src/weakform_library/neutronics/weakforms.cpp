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
  
  namespace Common { namespace WeakForms
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
    
  /* WeakForms */
  }
  /* Common */
  }
  
  namespace Diffusion { namespace WeakForms
  {   
    HomogeneousPart::HomogeneousPart(const Common::MaterialProperties::MaterialPropertyMaps* matprop,
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
                                                const Hermes::vector<Solution<double>*>& iterates,
                                                double initial_keff_guess, 
                                                GeomType geom_type ) 
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess)
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
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G(), &matprop, geom_type, initial_keff_guess, fission_regions)
    {
      homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
      add_forms_from_homogeneous_part();
      
      init_rhs(iterates);
    }
    
    void KeffEigenvalueProblem::init_rhs(const Hermes::vector<Solution<double>*>& iterates)
    {
      const MaterialPropertyMaps *mp = static_cast<const MaterialPropertyMaps*>(matprop);
      
      stored_flux_solutions.reserve(G);
      scalar_flux_iterates.reserve(G);
      for (unsigned int gto = 0; gto < G; gto++)
      { 
        stored_flux_solutions.push_back(iterates[gto]);
        scalar_flux_iterates.push_back(static_cast<MeshFunction<double>*>(stored_flux_solutions.back()));
      }
      
      for (unsigned int gto = 0; gto < G; gto++)
      { 
        FissionYield::OuterIterationForm* keff_iteration_form;
        
        if (fission_regions.size() > 0) 
          keff_iteration_form = new FissionYield::OuterIterationForm( gto, fission_regions, *mp, scalar_flux_iterates, keff, geom_type );
        else
          keff_iteration_form = new FissionYield::OuterIterationForm( gto, *mp, scalar_flux_iterates, keff, geom_type );
        
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
    HomogeneousPart::HomogeneousPart(unsigned int N_odd,
                                     const Common::MaterialProperties::MaterialPropertyMaps* matprop,
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
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {
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
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {      
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
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {      
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
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {
      if (minus_isotropic_sources.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
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
      : NeutronicsProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type),
        SPNWeakForm(N, matprop.get_G())
    {
      if (minus_isotropic_sources.size() != G)
        error_function(Messages::E_INVALID_SIZE);
      
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
      : Common::WeakForms::KeffEigenvalueProblem(matprop.get_G()*(N+1)/2, &matprop, geom_type, initial_keff_guess, fission_regions),
        SPNWeakForm(N, matprop.get_G())
    { 
      stored_flux_solutions.reserve(iterates.size());
      for (Hermes::vector<Solution<double>*>::const_iterator it = iterates.begin(); it != iterates.end(); ++it)
        stored_flux_solutions.push_back(*it);
      
      SupportClasses::MomentFilter::get_scalar_fluxes_with_derivatives(stored_flux_solutions, &scalar_flux_iterates, G);
      
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
      
      SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_flux_iterates);
    }
    
    void KeffEigenvalueProblem::update_keff(double new_keff) 
    { 
      keff = new_keff;
      
      std::vector<FissionYield::OuterIterationForm*>::const_iterator it = keff_iteration_forms.begin();
      for ( ; it != keff_iteration_forms.end(); ++it)
        (*it)->update_keff(new_keff); 
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