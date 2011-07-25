#include "neutronics/support_classes.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics { namespace SupportClasses
{  
  namespace Common
  {
    void SourceFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
    {
      int marker = this->get_active_element()->marker;
      std::string material = matprop.get_material(marker, mesh);
      
      memset(result, 0, n*sizeof(double));
      
      if (markers.empty() || markers.find(marker) != markers.end())
      {           
        DataStructures::rank1 Sigma_f = matprop.get_Sigma_f(material);
        DataStructures::rank1 nu = matprop.get_nu(material);
      
        for (int i = 0; i < n; i++) 
          for (unsigned int j = 0; j < values.size(); j++)
            result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
      }
    }
    
    void SourceFilter::pre_init()
    {
      num = matprop.get_G();
      if(num > 10)
        error_function("Unable to create an instance of SourceFilter: Hermes is currently able to handle"
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

    void SourceFilter::assign_solutions(const Hermes::vector< Solution<double>* >& solutions)
    {
      if (solutions.size() != (unsigned) num)
        error_function("SourceFilter: Number of solutions does not match the size of data.");
      
      free();
      for (int i = 0; i < num; i++)
        sln[i] = solutions[i];
      init();
      post_init();
    }

    void SourceFilter::assign_solutions(const Hermes::vector< MeshFunction<double>* >& solutions)
    {
      if (solutions.size() != (unsigned) num)
        error_function("SourceFilter: Number of solutions does not match the size of data.");
      
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
      SimpleFilter<double>::set_active_element(e);
      
      order = sln[0]->get_fn_order();
      for (int i = 1; i < num; i++)
        if (sln[i]->get_fn_order() > order)
          order = sln[i]->get_fn_order();
    }

    double SourceFilter::integrate()
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
          double *uval = this->get_fn_values();
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
        error_function("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
              "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
      
      if (m > n) std::swap(m,n);
      
      if (k_of_Sigma_t2k > n)
        error_function("In the m-th SPN equation, m = %d, the coefficients at the unknown generalized fluxes involve"
              "Sigma_tk with k up to %d. Entered k = %d.", n, 2*n, 2*k_of_Sigma_t2k);
      
      return SYSTEM_MATRIX[m][n-m][k_of_Sigma_t2k];
    }
    
    double Coeffs::D_grad_F(unsigned int m, unsigned int n)
    {
      if (m >= N_MAX || n >= N_MAX)
        error_function("For the maximum implemented SPN order (N = %d), both m, n must lie in the range [0,%d]."
              "Entered (m,n) = (%d,%d).", 2*N_MAX-1, N_MAX-1, m, n);
      
      if (m > n) std::swap(m,n); 
      return D_GRAD_F[m][n-m];
    }
    
    double Coeffs::even_moment(unsigned int m, unsigned int n)
    {
      if (m >= N_MAX)
        error_function("For the maximum implemented SPN order (N = %d), there are %d even moment(s)."
              "Tried to access moment #%d.", 2*N_MAX-1, N_MAX, m+1);
      if (n > N_MAX || n < m)
        error_function("The even moment #%d must be expressed in terms of %d odd moment(s), starting with %d."
              "Tried to use odd moment #%d.", m+1, N_MAX-m, m+1, n+1);
              
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
    
    void MomentFilter::Val::filter_fn(int n, Hermes::vector< double* > values, double* result)
    {
      if (odd_req_mom)
        memcpy(result, values.at(req_mom_idx), n*sizeof(double));
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
      SimpleFilter<double>::set_active_element(e);

      order = -1;
      
      unsigned int exp_mom_idx = req_mom_idx;
      for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    void MomentFilter::ValDxDy::filter_fn(int n, Hermes::vector<double *> values, 
                                          Hermes::vector<double *> dx, Hermes::vector<double *> dy, 
                                          double* rslt, double* rslt_dx, double* rslt_dy)
    {
      if (odd_req_mom)
      {
        memcpy(rslt, values.at(req_mom_idx), n*sizeof(double));
        memcpy(rslt_dx, dx.at(req_mom_idx), n*sizeof(double));
        memcpy(rslt_dy, dy.at(req_mom_idx), n*sizeof(double));
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
      DXDYFilter<double>::set_active_element(e);
      
      order = -1;
      
      unsigned int exp_mom_idx = req_mom_idx;
      for (int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < num; sol_idx = mg.pos(++exp_mom_idx,g))
        if (sln[sol_idx]->get_fn_order() > order)
          order = sln[sol_idx]->get_fn_order();
    }
    
    // TODO: Templatize.
    void MomentFilter::get_scalar_fluxes(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::Val(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::get_scalar_fluxes_with_derivatives(const Hermes::vector< Solution<double>* >& angular_fluxes, 
                                                          Hermes::vector< MeshFunction<double>* >* scalar_fluxes,
                                                          unsigned int G)
    {          
      scalar_fluxes->reserve(G);
      for (unsigned int g = 0; g < G; g++)
        scalar_fluxes->push_back(new MomentFilter::ValDxDy(0, g, G, angular_fluxes));
    }
    
    void MomentFilter::clear_scalar_fluxes(Hermes::vector< MeshFunction<double>* >* scalar_fluxes)
    {
      Hermes::vector< MeshFunction<double>* >::const_iterator it = scalar_fluxes->begin();
      for( ; it != scalar_fluxes->end(); ++it)
        delete *it;
      scalar_fluxes->clear();
    }
    
    void SourceFilter::filter_fn(int n, Hermes::vector< double* > values, double* result)
    {
      int marker = this->get_active_element()->marker;
      std::string material = matprop.get_material(marker, mesh);
      
      memset(result, 0, n*sizeof(double));
      
      if (markers.empty() || markers.find(marker) != markers.end())
      {           
        DataStructures::rank1 Sigma_f = matprop.get_Sigma_f(material);
        DataStructures::rank1 nu = matprop.get_nu(material);
        
        for (int i = 0; i < n; i++) 
        {
          for (unsigned int g = 0; g < G; g++)
          {
            double group_scalar_flux = 0;
            
            unsigned int exp_mom_idx = 0;
            for (unsigned int sol_idx = mg.pos(exp_mom_idx,g); sol_idx < values.size(); sol_idx = mg.pos(++exp_mom_idx,g))
              group_scalar_flux += Coeffs::even_moment(0, exp_mom_idx) * values.at(sol_idx)[i];
            
            result[i] += nu[g] * Sigma_f[g] * group_scalar_flux;
          }
        }
      }
    }
  }
  
  double PostProcessor::integrate(MeshFunction<double>* solution, const Hermes::vector<std::string>& areas) const
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
        double *uval = solution->get_fn_values();
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

  void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution<double>* >* solutions, double integrated_fission_source) const
  {
    if (integrated_fission_source < 1e-12)
      error_function("PostProcessor::normalize_to_unit_fission_source : Invalid fission source.");
 // TODO (missing Solution::multiply)
 /*
    Hermes::vector< Solution<double>* >::iterator sln = solutions->begin();
    for ( ; sln != solutions->end(); ++sln)
      (*sln)->multiply(1./integrated_fission_source);
*/
  }

  void PostProcessor::normalize_to_unit_fission_source(Hermes::vector< Solution<double>* >* solutions, 
                                                        const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                        const Hermes::vector< std::string >& src_areas) const
  {
    Common::SourceFilter *sf;
    
    if (method == NEUTRONICS_DIFFUSION)
      sf = new Common::SourceFilter(*solutions, matprop, src_areas, geom_type);
    else if (method == NEUTRONICS_SPN)
      sf = new SPN::SourceFilter(*solutions, matprop, src_areas, geom_type);
    
    normalize_to_unit_fission_source(solutions, sf->integrate());
    
    delete sf;
  }

  void PostProcessor::normalize_to_unit_power(Hermes::vector< Solution<double>* >* solutions, 
                                              const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                              double power_per_fission, const Hermes::vector< std::string >& src_areas) const
  {
    // TODO
  }
  
  double PostProcessor::get_integrated_group_reaction_rates_internal( ReactionType reaction, MeshFunction<double>* solution, 
                                                                      const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                                      const Hermes::vector< std::string >& regions,
                                                                      unsigned int this_group, int other_group) const
  {    
    if (this_group > matprop.get_G())
      error_function(Messages::E_INVALID_GROUP_INDEX);
    
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
        double *uval = solution->get_fn_values();
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
              error_function(Messages::E_INVALID_GROUP_INDEX);
            
            xsec = matprop.compute_Sigma_s(mat)[this_group][other_group];
            break;
          case SELF_SCATTERING:
            xsec = matprop.compute_Sigma_s(mat)[this_group][this_group];
            break;
          case OUT_SCATTERING:
            if (other_group > (int) matprop.get_G() || other_group < 0)
              error_function(Messages::E_INVALID_GROUP_INDEX);
            
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
                                                          const Hermes::vector< Solution<double>* >& solutions, Hermes::vector< double >* results, 
                                                          const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                          unsigned int group, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
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
  
  double PostProcessor::get_integrated_group_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                            const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                            unsigned int group, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
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
  
  void PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                    Hermes::vector< double >* results, 
                                                    const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                    const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, matprop.get_G());
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
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
  
  double PostProcessor::get_integrated_reaction_rates(ReactionType reaction, const Hermes::vector< Solution<double>* >& solutions, 
                                                      const MaterialProperties::Common::MaterialPropertyMaps& matprop, 
                                                      const Hermes::vector< std::string >& regions) const
  {
    double result = 0.0;
    for (unsigned int group = 0; group < matprop.get_G(); group++)
      result += get_integrated_group_reaction_rates(reaction, solutions, matprop, group, regions);
    return result;
  }

  void PostProcessor::get_integrated_group_scalar_fluxes( const Hermes::vector< Solution<double>* >& solutions, 
                                                          Hermes::vector< double >* results,
                                                          unsigned int group, unsigned int G,
                                                          const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    for ( ; region != regions.end(); ++region)
      results->push_back(integrate(scalar_fluxes[group], *region));
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
  }
  
  double PostProcessor::get_integrated_group_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, 
                                                            unsigned int group, unsigned int G,
                                                            const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
        scalar_fluxes.push_back(*it);
    
    double result = integrate(scalar_fluxes[group], regions);
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);
    
    return result;
  }
  
  void PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, Hermes::vector< double >* results, 
                                                    unsigned int G, const Hermes::vector< std::string >& regions) const
  {
    Hermes::vector<MeshFunction<double>*> scalar_fluxes;
    
    if (method == NEUTRONICS_SPN)
      SPN::MomentFilter::get_scalar_fluxes(solutions, &scalar_fluxes, G);
    else if (method == NEUTRONICS_DIFFUSION)
      for (Hermes::vector<Solution<double>*>::const_iterator it = solutions.begin(); it != solutions.end(); ++it)
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

  double PostProcessor::get_integrated_scalar_fluxes(const Hermes::vector< Solution<double>* >& solutions, unsigned int G,
                                                      const Hermes::vector< std::string >& regions) const
  {
    double result = 0.0;
    for (unsigned int group = 0; group < G; group++)
      result += get_integrated_group_scalar_fluxes(solutions, group, G, regions);
    return result;
  }
  
  void PostProcessor::get_areas(Mesh *mesh, const Hermes::vector<std::string>& regions, Hermes::vector<double>* results) const
  {
    Solution<double> unity(mesh, 1);
    
    Hermes::vector<std::string>::const_iterator region = regions.begin(); 
    for ( ; region != regions.end(); ++region)
      results->push_back(integrate(&unity, *region));
  }

  double PostProcessor::get_area(Mesh* mesh, const Hermes::vector< std::string >& regions) const
  {
    Solution<double> unity(mesh, 1);
    return integrate(&unity, regions);
  }    
  
  namespace Common
  {
    const std::string Visualization::base_title_flux = "Neutron flux: group ";
    const std::string Visualization::base_title_order = "Polynomial orders: group ";
    const std::string Visualization::base_title_mesh = "Core mesh for group ";
    
    void Visualization::init(unsigned int nu, unsigned int ne, unsigned int ng) 
    {
      n_unknowns = nu; n_equations = ne; n_groups = ng;
      
      if (nu > 0 && ne > 0 && ng > 0)
      {
        sviews = new Views::ScalarView<double>* [n_unknowns];
        oviews = new Views::OrderView<double>* [n_equations];
        if (display_meshes)
          mviews = new Views::MeshView* [n_equations];
      }
      else
      {
        sviews = NULL;
        oviews = NULL;
        mviews = NULL;
      }
    }
    
    Visualization::~Visualization()
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
    
    void Visualization::inspect_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
      {
        show_meshes(meshes);
        Views::View::wait();
        
        for (unsigned int i = 0; i < n_equations; i++)
          delete mviews[i];
        delete [] mviews;
        
        mviews = NULL;
      }
    }
    
    void Visualization::inspect_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      show_solutions(solutions);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_unknowns; i++)
        delete sviews[i];
      delete [] sviews;
      
      sviews = NULL;
    }
    
    void Visualization::inspect_orders(Hermes::vector< Space<double>* > spaces)
    {
      show_orders(spaces);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_equations; i++)
        delete oviews[i];
      delete [] oviews;
      
      oviews = NULL;
    }
  }
  
  namespace Diffusion
  {
    Visualization::Visualization(unsigned int G, bool display_meshes) : Common::Visualization(G, G, G, display_meshes)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string title_flux = base_title_flux + itos(g);
        std::string title_order = base_title_order + itos(g);
        
        sviews[g] = new Views::ScalarView<double>(title_flux.c_str(), new Views::WinGeom(0, g*452, 450, 450));
        sviews[g]->show_mesh(false);
        sviews[g]->set_3d_mode(true);
        oviews[g] = new Views::OrderView<double>(title_order.c_str(), new Views::WinGeom(0, n_groups*452 + g*452, 450, 450));
      }
      
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
        {
          std::string title = base_title_mesh + itos(g);
          mviews[g] = new Views::MeshView(title.c_str(), new Views::WinGeom(0, g*352, 350, 350));
        }
    }
    
    void Visualization::show_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
          mviews[g]->show(meshes[g]);
    }
    
    void Visualization::show_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        sviews[g]->show(solutions[g]);
    }
    
    void Visualization::show_orders(Hermes::vector< Space<double>* > spaces)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        oviews[g]->show(spaces[g]);
    }
    
    void Visualization::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                    Hermes::vector< Solution<double>* > solutions, bool mode_3D)
    {
      Views::Linearizer<double> lin;
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string appendix = std::string("_group_") + itos(g);
        std::string file = base_filename + appendix + std::string(".vtk");
        std::string var = base_varname + appendix;
        lin.save_solution_vtk(solutions[g], file.c_str(), var.c_str(), mode_3D);
        info("Scalar flux in group %d saved in VTK format to file %s.", g, file.c_str());
      }
    }
    
    void Visualization::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space<double>* > spaces)
    {
      Views::Orderizer ord;
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
    Visualization::Visualization(unsigned int spn_order, unsigned int G, bool display_meshes) : Common::Visualization(display_meshes), mg(G)
    {
      n_moments = spn_order+1;
      n_odd_moments = (n_moments+1)/2;
      Common::Visualization::init(G * n_moments, G * n_odd_moments, G);
      
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string title_flux = base_title_flux + itos(g) + std::string(", moment ");
        std::string title_order = base_title_order + itos(g) + std::string(", moment ");
        for (unsigned int m = 0; m < n_moments; m++)
        {
          unsigned int i = mg.pos(m,g);
          
          sviews[i] = new Views::ScalarView<double>((title_flux + itos(m)).c_str(), new Views::WinGeom(m*452, g*452, 450, 450));
          sviews[i]->show_mesh(false);
          sviews[i]->set_3d_mode(true);
          
          if (m%2) 
            oviews[mg.pos((m-1)/2,g)] = new Views::OrderView<double>((title_order + itos(m)).c_str(), new Views::WinGeom(m*452, n_groups*452 + g*452, 450, 450));
        }
      }
      
      if (display_meshes)
      {
        for (unsigned int g = 0; g < n_groups; g++)
        {
          std::string title = base_title_mesh + itos(g) + std::string(", moment ");
          for (unsigned int m = 0; m < n_odd_moments; m++)
            mviews[mg.pos(m,g)] = new Views::MeshView((title + itos(m)).c_str(), new Views::WinGeom(m*352, g*352, 350, 350));
        }
      }
    }
            
    void Visualization::show_meshes(Hermes::vector< Mesh* > meshes)
    {
      if (display_meshes)
        for (unsigned int g = 0; g < n_groups; g++)
          for (unsigned int m = 0; m < n_odd_moments; m++)
            mviews[mg.pos(m,g)]->show(meshes[mg.pos(m,g)]);
    }
    
    void Visualization::show_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      for (unsigned int g = 0; g < n_groups; g++)
      {
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          unsigned int i = mg.pos(m,g);
          unsigned int j = mg.pos(2*m,g);
          unsigned int k = mg.pos(2*m+1,g);
          
          MomentFilter::Val mf(2*m, g, n_groups, solutions);
          
          sviews[j]->show(&mf);
          sviews[k]->show(solutions[i]);
        }
      }
    }
    
    void Visualization::save_solutions_vtk(const std::string& base_filename, const std::string& base_varname, 
                                    Hermes::vector< Solution<double>* > solutions, bool mode_3D)
    { 
      Views::Linearizer<double> lin;
      for (unsigned int g = 0; g < n_groups; g++)
      {
        std::string appendix = std::string("_group_") + itos(g);
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          MomentFilter::Val mf(2*m, g, n_groups, solutions);
          
          std::string file = base_filename + std::string("_moment_") + itos(2*m) + appendix + std::string(".vtk");
          std::string var = base_varname + std::string("_moment_") + itos(2*m) + appendix;
          lin.save_solution_vtk(&mf, file.c_str(), var.c_str(), mode_3D);
          info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m, g, file.c_str());
          
          file = base_filename + std::string("_moment_") + itos(2*m+1) + appendix + std::string(".vtk");
          var = base_varname + std::string("_moment_") + itos(2*m+1) + appendix;
          lin.save_solution_vtk(solutions[mg.pos(m,g)], file.c_str(), var.c_str(), mode_3D);
          info("SP%d moment #%d of solution in group %d saved in VTK format to file %s.", n_moments-1, 2*m+1, g, file.c_str());
        }
      }
    }
    
    void Visualization::save_orders_vtk(const std::string& base_filename, Hermes::vector< Space<double>* > spaces)
    {
      Views::Orderizer ord;
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_odd_moments; m++)
        {
          std::string file = base_filename + std::string("_moment_") + itos(2*m+1) + std::string("_group_") + itos(g) + std::string(".vtk");
          ord.save_orders_vtk(spaces[mg.pos(m,g)], file.c_str());
          info("Information about approximation space for moment %d, group %d saved in VTK format to file %s.", 2*m+1, g, file.c_str());
        }
    }
    
    void Visualization::inspect_solutions(Hermes::vector< Solution<double>* > solutions)
    {
      show_solutions(solutions);
      Views::View::wait();
      
      for (unsigned int i = 0; i < n_unknowns; i++)
        delete sviews[i];
      delete [] sviews;
      
      sviews = NULL;
    }

    void Visualization::show_orders(Hermes::vector< Space<double>* > spaces)
    {
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int m = 0; m < n_odd_moments; m++)
          oviews[mg.pos(m,g)]->show(spaces[mg.pos(m,g)]);
    }
  }


  
/* SupportClasses */
}
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
} 