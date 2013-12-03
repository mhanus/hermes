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
                       const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries,
                       const char* out_tensor) 
  : WeakForm<double>(matprop.get_G()*N*(N+2)/2), N(N), M(N*(N+2)/2), G(matprop.get_G()), odata(SupportClasses::OrdinatesData(N, "../lgvalues.txt"))
{  
    bool1 chi_nnz = matprop.get_fission_nonzero_structure();
   
    //std::cout << odata;
    
    bool assemble_all = strlen(out_tensor) == 0;    
    bool assemble_A = !strcmp(out_tensor, "A");
    bool assemble_Q = !strcmp(out_tensor, "Q");
    bool assemble_S = !strcmp(out_tensor, "S");

    if (!assemble_all && !strcmp(out_tensor, "AQ"))
    {
      assemble_A = true;
      assemble_Q = true;
    }
    
    for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
    {
      for (unsigned int n = 0; n < M; n++)
      {
        if (assemble_all || assemble_A)
        {
          if (!reflective_boundaries.empty())
          {
            add_matrix_form_surf(new SpecularReflectionMF_X(odata, n, gto, G, reflective_boundaries));
            add_matrix_form_surf(new SpecularReflectionMF_Y(odata, n, gto, G, reflective_boundaries));
          }
          add_matrix_form_surf(new BoundaryStreamingMF(n, gto, G));
        }
        if (!inflow_boundaries.empty() && (assemble_all || assemble_Q))
          add_vector_form_surf(new BoundaryStreamingVF(n, gto, G, inflow_boundaries));
      }
    }
    
    std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
    for ( ; material != matprop.get_materials_list().end(); ++material)
    {
      Hermes::vector<std::string> regions = matprop.get_regions(*material);
            
      rank1 Sigma_t = matprop.get_Sigma_t(*material);
      rank3 Sigma_sn = matprop.get_Sigma_sn(*material);
      rank1 chi = matprop.get_chi(*material);
      rank1 nu = matprop.get_nu(*material);
      rank1 src_data = matprop.get_iso_src(*material);
     
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
      {
        for (unsigned int n = 0; n < M; n++)
        { 
          if (assemble_all || assemble_A)
            add_matrix_form(new VolumetricStreamingAndReactionsMF(regions, n, gto, G, Sigma_t[gto]));
          
          if (assemble_S)
        	for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
        		for (unsigned int m = 0; m < M; m++)
        			add_matrix_form(new VolumetricScatteringMF(odata, regions, m, n, gto, gfrom, G, Sigma_sn, Sigma_t[gto]));

          if (!Sigma_sn.empty() && assemble_all)
            add_vector_form(new VolumetricScatteringSourceVF(odata, regions, n, gto, G, Sigma_sn, Sigma_t[gto]));
          if (src_data[gto] > 0 && (assemble_all || assemble_Q)) 
            add_vector_form(new VolumetricExternalSourceVF(odata, regions, n, gto, G, Sigma_sn, Sigma_t[gto], src_data[gto], true));
        }
      }   
    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
/*	  result += wt[quad_pt] * ( 1/Sigma_t
			  	  	  	  	   * static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt])
			  	  	  	  	   * static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt])
			  	  	  	  	   + Sigma_t * u->val[quad_pt] * v->val[quad_pt] ); */
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
/*
template<typename Real>
Real SNWeakForm::VolumetricStreamingAndReactionsMF::b(Real x, Real y) const
{
  return Real(1);
}
*/

template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *sln_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      double C = (2*l + 1) / (4*M_PI);
      
      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
          odata.ordinates_to_moment(l, m, gfrom, G, u_ext, n, sln_moment_values_at_quad_pts);
          
          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
            
          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real tmp = sln_moment_values_at_quad_pts[quad_pt] * rlm_in_dir * C * Sigma_sn[l][gto][gfrom]
                                               * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
                                                   Sigma_t * v->val[quad_pt] );            
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
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return /*1/Sigma_t**/result;
}


template<typename Real>
Real SNWeakForm::VolumetricScatteringMF::matrix_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *trialfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);

    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      double C = (2*l + 1) / (4*M_PI);

      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
        	odata.ordinate_to_moment(direction, l, m, gfrom, G, u, n, trialfn_moment_values_at_quad_pts);


          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);

          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real tmp = trialfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir * C * Sigma_sn[l][gto][gfrom]
                                               * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
                                                   Sigma_t * v->val[quad_pt] );
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

  delete [] trialfn_moment_values_at_quad_pts;

  return /*1/Sigma_t**/result;
}

/*
template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, Func< Real >** ext) const
{
  Real result(0.0);
  Real *sln_moment_values_at_quad_pts = new Real [n];
  Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      double C = (2*l + 1) / (4*M_PI);
      
      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
          odata.ordinates_to_moment(l, m, gfrom, G, u_ext, n, sln_moment_values_at_quad_pts);
          odata.ordinate_to_moment(direction, l, m, gto, G, v, n, testfn_moment_values_at_quad_pts);
          
          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
            
          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real sln_group_source_moment = sln_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
            Real testfn_group_source_moment = testfn_moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
            Func<Real>* u = u_ext[ag.pos(direction,gfrom)];
            
            Real tmp = sln_group_source_moment * C * Sigma_sn[l][gto][gfrom]
                                               * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) +
                                                   Sigma_t * v->val[quad_pt] );
            tmp += testfn_group_source_moment * C * Sigma_sn[l][gto][gto]
                                              * ( _wf->calculate_a_dot_v(direction, u->dx[quad_pt], u->dy[quad_pt]) +
                                                  Sigma_t * u->val[quad_pt] );
            tmp -= testfn_group_source_moment * C * Sigma_sn[l][gto][gto] * sln_group_source_moment * C * Sigma_sn[l][gto][gfrom];
            
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
*/
template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
            
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {            
    Real tmp = Q * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) + Sigma_t * v->val[quad_pt] );
    
    // Contribution to spatial integral.
    if (geom_type == HERMES_AXISYM_X)
      result += tmp * wt[quad_pt] * e->y[quad_pt];
    else if (geom_type == HERMES_AXISYM_Y)
      result += tmp * wt[quad_pt] * e->x[quad_pt];
    else
      result += tmp * wt[quad_pt];
  }

  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return /*1/Sigma_t**/result;
}

/*
template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Real> **ext) const
{
  Real result = Real(0.0);
  Real *testfn_moment_values_at_quad_pts = new Real [n];
  SNWeakForm* _wf = static_cast<SNWeakForm*>(wf);
  
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
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
            
            Real tmp = Q * ( _wf->calculate_a_dot_v(direction, v->dx[quad_pt], v->dy[quad_pt]) + Sigma_t * v->val[quad_pt]
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
      
      result += deg_l_result;
    }
  }
    
  delete [] testfn_moment_values_at_quad_pts;
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result;
}
*/
/*
double SNWeakForm::VolumetricExternalSourceVF::Q(double x, double y) const
{
  //return 1 + (3*M_PI*(1 + x)*(1 + y)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/20. + (M_PI*std::pow(1 + y,2)*std::cos((M_PI*(1 + x)*std::pow(1 + y,2))/8.))/10. + std::sin((M_PI*(1 + x)*std::pow(1 + y,2))/8.);
 
  return 6;
}

Ord SNWeakForm::VolumetricExternalSourceVF::Q(Ord x, Ord y) const
{
  return Ord(0);
}
*/

double SNWeakForm::SpecularReflectionMF_X::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
  
    if (a_dot_n < 0)
      if (fabs(e->ny[quad_pt] - 1.0) < eps || (fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps))
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

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, Func<double> **ext) const
{
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
      result += wt[quad_pt] * u->val[quad_pt] * (-a_dot_n) * v->val[quad_pt];
  }
  
  //std::cout << "BoundaryStreamingMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::BoundaryStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, Func<Ord> **ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
    {
      double boundary_data = influx<double>(e->x[quad_pt],e->y[quad_pt]);
      result += wt[quad_pt] * boundary_data * (-a_dot_n) * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::BoundaryStreamingVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  //std::cout << "BoundaryStreamingVF :: " << v->val[0] * influx<Ord>(e->x[0], e->y[0]) << std::endl;
  return v->val[0] * influx<Ord>(e->x[0], e->y[0]);
}

double SNWeakForm::SpecularReflectionVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double eps = 1e-14;
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
   
    if (a_dot_n < 0)
    {
      double boundary_data = 0.0;
      
      if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
        boundary_data = u_ext[ag.pos(odata.reflections_about_x[direction],g)]->val[quad_pt];
      else if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        boundary_data = u_ext[ag.pos(odata.reflections_about_y[direction],g)]->val[quad_pt];
      else
        Hermes::Mixins::Loggable::Static::warn("Only horizontal or vertical boundaries are currently supported for specular reflection b. c.");

      result += wt[quad_pt] * boundary_data * (-a_dot_n) * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::SpecularReflectionVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  int refl_ord = std::max(u_ext[ag.pos(odata.reflections_about_x[direction],g)]->val[0].get_order(),u_ext[ag.pos(odata.reflections_about_y[direction],g)]->val[0].get_order());
  
  //std::cout << "SpecularReflectionVF :: " << static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), Ord(refl_ord), Ord(1)) * v->val[0] << std::endl;
  
  return Ord(refl_ord) * v->val[0];
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
