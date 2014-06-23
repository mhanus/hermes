#include "definitions.h"
#include <cstdlib>

IsotropicScatteringAndFissionMatrixForms::IsotropicScatteringAndFissionMatrixForms(const MaterialProperties::MaterialPropertyMaps& matprop, const char* out_tensor) : WeakForm<double>(matprop.get_G())
{
  bool assemble_S = false;
  bool assemble_F = false;
  
  if (!strcmp(out_tensor, "S"))
    assemble_S = true;
  else if (!strcmp(out_tensor, "F"))
    assemble_F = true;
  else
    error("Wrong specification of matrix to save - use either S or F");
    
  bool1 chi_nnz = matprop.get_fission_nonzero_structure();
  
  std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
  for ( ; material != matprop.get_materials_list().end(); ++material)
  {
    Hermes::vector<std::string> regions = matprop.get_regions(*material);
          
    rank3 Sigma_sn = matprop.get_Sigma_sn(*material);
    rank1 Sigma_f = matprop.get_Sigma_f(*material);
    rank1 chi = matprop.get_chi(*material);
    rank1 nu = matprop.get_nu(*material);
    
    for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
    {
      for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
      {
        if (!Sigma_sn.empty() && assemble_S)                      
          add_matrix_form(new Diffusion::WeakFormParts::Scattering::Jacobian(regions, gto, gfrom, -Sigma_sn[0][gto][gfrom]));
        if (chi_nnz[gto] && assemble_F)
          add_matrix_form( new Diffusion::WeakFormParts::FissionYield::Jacobian(regions, gto, gfrom, 
                                                                                chi[gto], -nu[gfrom], Sigma_f[gfrom]) );
      }
      
      add_vector_form(new Diffusion::WeakFormParts::ExternalSources::LinearForm(regions, gto, -1.0));
    }   
  }
}

SNWeakForm::SNWeakForm(unsigned int N, const MaterialProperties::MaterialPropertyMaps& matprop,
                       const Hermes::vector<std::string>& reflective_boundaries,
                       const char* out_tensor) 
  : WeakForm<double>(matprop.get_G()*N*(N+2)/2), N(N), M(N*(N+2)/2), G(matprop.get_G()), odata(SupportClasses::OrdinatesData(N, "../lgvalues.txt"))
{  
    bool1 chi_nnz = matprop.get_fission_nonzero_structure();
   
    //std::cout << odata;
    
    bool assemble_all = strlen(out_tensor) == 0;    
    bool assemble_L = !strcmp(out_tensor, "L");
    
    for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
    {
      for (unsigned int n = 0; n < M; n++)
      {
        if (assemble_all || assemble_L)
        {
          if (!reflective_boundaries.empty())
          {
            add_matrix_form_surf(new SpecularReflectionMF_X(odata, n, gto, G, reflective_boundaries));
            add_matrix_form_surf(new SpecularReflectionMF_Y(odata, n, gto, G, reflective_boundaries));
          }
        }
      }
    }
    
    std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
    for ( ; material != matprop.get_materials_list().end(); ++material)
    {
      Hermes::vector<std::string> regions = matprop.get_regions(*material);
            
      rank1 Sigma_t = matprop.get_Sigma_t(*material);
      rank3 Sigma_sn = matprop.get_Sigma_sn(*material);
      rank1 Sigma_f = matprop.get_Sigma_f(*material);
      rank1 chi = matprop.get_chi(*material);
      rank1 nu = matprop.get_nu(*material);
     
      for (unsigned int gto = 0; gto < G; gto++)
      {
        for (unsigned int n = 0; n < M; n++)
        { 
          if (assemble_all || assemble_L)
            add_matrix_form(new VolumetricStreamingAndReactionsMF(regions, n, gto, G, Sigma_t[gto]));

          for (unsigned int gfrom = 0; gfrom < G; gfrom++)
          {
            if (chi_nnz[gto] && assemble_all)
              add_vector_form(new VolumetricFissionSourceVF(odata, regions, n, gto, gfrom, chi[gto], nu, Sigma_f, Sigma_sn, Sigma_t[gto]));
            if (!Sigma_sn.empty() && assemble_all)
              for (unsigned int m = 0; m < M; m++)
                add_matrix_form(new VolumetricScatteringSourceMF(odata, regions, m, n, gto, gfrom, Sigma_sn, Sigma_t[gto]));
          }
        }
      }   
    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
    result += wt[quad_pt] * ( static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt]) +
                        Sigma_t * u->val[quad_pt] )
                    * ( static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
                        Sigma_t * v->val[quad_pt] );
  
  //std::cout << "VolumetricStreamingAndReactionsMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::VolumetricStreamingAndReactionsMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                  Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return (u->dx[0]+u->dy[0]+u->val[0]) * (v->dx[0]+v->dy[0]+v->val[0]);
}


template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  Real result(0.0);
  Real *sln_moment_values_at_quad_pts = new Real [n];
  Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  
  for (unsigned int l = 0; l <= L; l++)
  {
    Real deg_l_result(0.0);
    double C = (2*l + 1) / (4*M_PI);
    
    for (int m = -l; m <= l; m++)
    {
      if ( ((l + m) % 2) == 0 )
      {
        odata.ordinate_to_moment(dfrom, l, m, u, n, sln_moment_values_at_quad_pts);
        odata.ordinate_to_moment(dto, l, m, v, n, testfn_moment_values_at_quad_pts);
        
        SupportClasses::SphericalHarmonic Rlm(l, m);
        double rlm_in_dto = Rlm(odata.xi[dto], odata.eta[dto], odata.mu[dto]);
        double rlm_in_dfrom = Rlm(odata.xi[dfrom], odata.eta[dfrom], odata.mu[dfrom]);
          
        for (int quad_pt = 0; quad_pt < n; quad_pt++)
        {
          Real sln_group_source_moment = sln_moment_values_at_quad_pts[quad_pt] * rlm_in_dfrom;
          Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dto;
          
          Real tmp = - sln_group_source_moment * C * Sigma_sn[l][gto][gfrom]
                                                * ( _wf->calculate_a_dot_v(dto, v->dx[quad_pt], v->dy[quad_pt]) +
                                                    Sigma_t * v->val[quad_pt] );
          tmp -= testfn_group_source_moment * C * Sigma_sn[l][gto][gfrom]
                                            * ( _wf->calculate_a_dot_v(dto, u->dx[quad_pt], u->dy[quad_pt]) +
                                                Sigma_t * u->val[quad_pt] );
          tmp += testfn_group_source_moment * C * Sigma_sn[l][gto][gto] * sln_group_source_moment * C * Sigma_sn[l][gto][gfrom];
          
          // Contribution to spatial integral.
          if (geom_type == HERMES_AXISYM_X)
            deg_l_result += tmp * wt[quad_pt] * e->y[quad_pt];
          else if (geom_type == HERMES_AXISYM_Y)
            deg_l_result += tmp * wt[quad_pt] * e->x[quad_pt];
          else
            deg_l_result += tmp * wt[quad_pt];
        }
      }
    }
      
      result += deg_l_result;
    }
  }
    
  delete [] sln_moment_values_at_quad_pts;
  delete [] testfn_moment_values_at_quad_pts;
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return result;
}


template<typename Real>
Real SNWeakForm::VolumetricFissionSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  Real *testfn_moment_values_at_quad_pts = new Real [n];
  Real *group_scalar_flux = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (int quad_pt = 0; quad_pt < n; quad_pt++)
    {
      // Angular integration to get scalar flux at quadrature point.
      Real group_scalar_flux(0.0);
      for (int dir = 0; dir < odata.M; dir++)
        group_scalar_flux[quad_pt] += odata.pw[dir] * u_ext[ag.pos(dir,gfrom)]->val[quad_pt];
      //group_scalar_flux *= 2 / (4*M_PI) * 2*M_PI;
      group_scalar_flux[quad_pt] *= 2;
      
      // Compute fission source with element-wise constant nuSigma_f.
      group_scalar_flux[quad_pt] *= nu[gfrom] * Sigma_f[gfrom];
    }
    
    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      double C = (2*l + 1) / (4*M_PI);
      
      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
          odata.ordinate_to_moment(direction, l, m, gto, G, v, n, testfn_moment_values_at_quad_pts);
          
          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
            
          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
            Func<Real>* u = u_ext[ag.pos(direction,gfrom)];
            
            Real tmp = group_scalar_flux[quad_pt] * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) + Sigma_t * v->val[quad_pt]
                             - C * Sigma_sn[l][gto][gto] * testfn_group_source_moment );
            
            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              deg_l_result += tmp * wt[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              deg_l_result += tmp * wt[quad_pt] * e->x[quad_pt];
            else
              deg_l_result += tmp * wt[quad_pt];
          }
        }
      }
    }
  }
    
  delete [] testfn_moment_values_at_quad_pts;
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result * chi_to;
}

double SNWeakForm::SpecularReflectionMF_X::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
  
    if (a_dot_n < 0)
      if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        result += -wt[quad_pt] * u->val[quad_pt] * (-a_dot_n) * v->val[quad_pt];
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_X::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::SpecularReflectionMF_Y::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
      if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
        result += -wt[quad_pt] * u->val[quad_pt] * (-a_dot_n) * v->val[quad_pt];
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_Y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::calculate_a_dot_v(int n, double vx, double vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return odata.xi[n]*vx + odata.eta[n]*vy;
}

Ord SNWeakForm::calculate_a_dot_v(int n, Ord vx, Ord vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return vx + vy;
}
