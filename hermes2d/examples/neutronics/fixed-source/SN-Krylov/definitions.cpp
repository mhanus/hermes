#include "definitions.h"
#include <cstdlib>

DiffusionWeakForm::DiffusionWeakForm(const Mesh* mesh, bool DG, const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop, 
                                     const Hermes::vector< std::string >& void_boundaries) : Common::WeakForms::NeutronicsProblem(matprop.get_G(), &matprop, HERMES_PLANAR)
{
  homogeneous_part = new HomogeneousPart(&matprop, geom_type, false);
  add_forms_from_homogeneous_part();
        
  std::set<std::string>::const_iterator material = matprop.get_materials_list().begin();
  for ( ; material != matprop.get_materials_list().end(); ++material)
  {
    Hermes::vector<std::string> regions = matprop.get_regions(*material);
    
    rank1 src_data = matprop.get_iso_src(*material);
    
    for (unsigned int gto = 0; gto < G; gto++)
      add_vector_form(new Diffusion::WeakFormParts::ExternalSources::LinearForm(regions, gto, -src_data[gto], geom_type));
  }
  
  if (DG) 
  {  
    int theta = -1;
    int C_W = 100;
    
    for (unsigned int g = 0; g < G; g++)
    {
      if (!void_boundaries.empty())
        add_matrix_form_surf(new Diffusion::WeakFormParts::VacuumBoundaryCondition::Jacobian(void_boundaries, g));    
      
      //add_matrix_form_surf(new BoundaryJacobian (mesh, &matprop, theta, C_W));
      /*add_vector_form_surf(new BoundaryResidual (theta, C_W));*/
      add_matrix_form_surf(new InterfaceJacobian(mesh, &matprop, theta, C_W));
    }
  }
}

DiffusionWeakForm::HomogeneousPart::HomogeneousPart(const Common::MaterialProperties::MaterialPropertyMaps* matprop, GeomType geom_type, bool include_fission)
{
  const Diffusion::MaterialProperties::MaterialPropertyMaps *mp = static_cast<const Diffusion::MaterialProperties::MaterialPropertyMaps*>(matprop);
  
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
      matrix_forms.push_back(new Diffusion::WeakFormParts::DiffusionReaction::Jacobian(regions, gto, D[gto], Sigma_r[gto], geom_type));
      
      for (unsigned int gfrom = 0; gfrom < mp->get_G(); gfrom++)
      {
        if (Ss_nnz[gto][gfrom] && gto != gfrom)
        {
          matrix_forms.push_back(new Diffusion::WeakFormParts::Scattering::Jacobian(regions, gto, gfrom, Sigma_s[gto][gfrom], geom_type));
        }
        
        if (include_fission && chi_nnz[gto])
        {
          matrix_forms.push_back(new Diffusion::WeakFormParts::FissionYield::Jacobian(regions, gto, gfrom, 
                                                            chi[gto], nu[gfrom], Sigma_f[gfrom], 
                                                            geom_type));
        }
      }
    }
  }
}

#define JUMP(w)       ( w->get_val_central(i) - w->get_val_neighbor(i) )
#define AVG_GRAD(w)   ( 0.5 * ( (w->get_dx_central(i) + w->get_dx_neighbor(i))*e->nx[i] + \
                                (w->get_dy_central(i) + w->get_dy_neighbor(i))*e->ny[i] ) )

double DiffusionWeakForm::BoundaryJacobian::value(int n, double* wt, Func< double >* u_ext[], Func< double >* u, Func< double >* v, Geom< double >* e, ExtData< double >* ext) const
{
  double result(0.0);
  
  std::string marker_e = matprop->get_material(e->elem_marker, mesh);
  double D = matprop->get_D(marker_e)[0];
  
  double edge_len (0.);
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  double sigma = C_W * D / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * D * ( /*-dot2<double>(u->dx[i], u->dy[i], e->nx[i], e->ny[i]) * v->val[i] + */
                             dot2<double>(v->dx[i], v->dy[i], e->nx[i], e->ny[i]) * theta * u->val[i] );   
    //result += wt[i] * sigma * u->val[i] * v->val[i];                    
  }
  
  return result;
}

double DiffusionWeakForm::BoundaryResidual::value(int n, double* wt, Func< double >* u_ext[], Func< double >* v, Geom< double >* e, ExtData< double >* ext) const
{
  double result(0.0);
  double sigma = C_W * D / e->diam;
  //double edge_len = 0.;
  //for (int i = 0; i < n; i++)
  //  edge_len += wt[i];
  
  //double sigma = C_W * D / (0.5*edge_len);
  
  //std::string marker = wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker);
  
  for (int i = 0; i < n; i++)
  {
    //WeaklyImposableBC* bc = static_cast<WeaklyImposableBC*>(boundary_values.get_boundary_condition(marker));
    double bnd_diff;// = u_ext[0]->val[i] - bc->value(e->x[i], e->y[i]);
    result += wt[i] * D * ( -dot2<double>(u_ext[0]->dx[i], u_ext[0]->dy[i], e->nx[i], e->ny[i]) * v->val[i]
                                  +dot2<double>(v->dx[i], v->dy[i], e->nx[i], e->ny[i]) * theta * bnd_diff );   
    result += wt[i] * sigma * bnd_diff * v->val[i];                        
  }
  
  return result;
}

double DiffusionWeakForm::InterfaceJacobian::value(int n, double* wt, Hermes::Hermes2D::Func< double >* u_ext[], Func< double >* u, Func< double >* v, Geom< double >* e, Hermes::Hermes2D::ExtData< double >* ext) const
{
  double result(0.0);
  
  std::string marker_e = matprop->get_material(e->elem_marker, mesh);
//  std::cout << "E2  " << marker_e << std::endl;
  std::string marker_n = matprop->get_material(e->get_neighbor_marker(), mesh);
//  std::cout << "N2" << marker_n << std::endl;
  
  double D = 2 * matprop->get_D(marker_e)[0] * matprop->get_D(marker_n)[0] / (matprop->get_D(marker_n)[0] + matprop->get_D(marker_e)[0]);
  
  double edge_len(0.);
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  double sigma = C_W * D / (0.5*edge_len);
  
  //double sigma = D * C_W / (e->diam + e->get_neighbor_diam());
  
//  std::cout << "S-" << sigma << std::endl;
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * D * (-AVG_GRAD(u) * JUMP(v) + theta * AVG_GRAD(v) * JUMP(u)); // diffusion
    result += wt[i] * sigma * JUMP(u) * JUMP(v);                                          // interior discontinuity penalization
  }
  
//  std::cout << "R-" << result << std::endl;
  
  return result;
}







SNWeakForm::SNWeakForm(unsigned int N, const MaterialProperties::MaterialPropertyMaps& matprop, const Hermes::vector<Solution<double>*>& iterates,
                       const Hermes::vector<std::string>& reflective_boundaries, const Hermes::vector<std::string>& inflow_boundaries) 
  : WeakForm<double>(N*(N+2)/2), mesh(iterates[0]->get_mesh()), N(N), M(N*(N+2)/2), G(matprop.get_G()), odata(SupportClasses::OrdinatesData(N, "lgvalues.txt"))
{  
    bool1 chi_nnz = matprop.get_fission_nonzero_structure();
   
    std::cout << odata;
    
    Hermes::vector<MeshFunction<double>*> mfn_iterates;
    mfn_iterates.resize(M*G);
    std::copy(iterates.begin(), iterates.end(), mfn_iterates.begin());
    
    for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
    {
      for (unsigned int n = 0; n < M; n++)
      {
        add_matrix_form_surf(new InterfaceStreamingMF(n, gto, G));
        if (!reflective_boundaries.empty())
        {
          add_matrix_form_surf(new SpecularReflectionMF_X(odata, n, gto, G));
          add_matrix_form_surf(new SpecularReflectionMF_Y(odata, n, gto, G));
        }
        add_matrix_form_surf(new BoundaryStreamingMF(n, gto, G));
        if (!inflow_boundaries.empty())
          add_vector_form_surf(new BoundaryStreamingVF(n, gto, G, inflow_boundaries));
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
      rank1 src_data = matprop.get_iso_src(*material);
     
      for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
      {
        for (unsigned int n = 0; n < M; n++)
        {        
          add_matrix_form(new VolumetricStreamingAndReactionsMF(regions, n, gto, G, Sigma_t[gto]));
          
          if (src_data[gto] > 0) 
            add_vector_form(new VolumetricExternalSourceVF(regions, n, gto, G, src_data[gto], true));
        }
      }
    }
}

double SNWeakForm::VolumetricStreamingAndReactionsMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0.0;
  //std::cout << "VSR " <<  n << std::endl;
  for (int i = 0; i < n; i++) {
    result += wt[i] * u->val[i] * ( -static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, v->dx[i], v->dy[i]) + /*b<Real,Real>(e->x[i],e->y[i])*/
                                    Sigma_t * v->val[i] );
  
 // std::cout << "VolumetricStreamingAndReactionsMF :: " << i << "(" << v->dx[i] << "," << v->dy[i] << ")" << std::endl;
  }
  return result;
}

Ord SNWeakForm::VolumetricStreamingAndReactionsMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                  Geom<Ord> *e, ExtData<Ord> *ext) const
{ 
  return u->val[0] * (v->dx[0]+v->dy[0]+v->val[0]);
}
/*
template<typename Real>
Real SNWeakForm::VolumetricStreamingAndReactionsMF::b(Real x, Real y) const
{
  return Real(1);
}
*/

template<typename Real>
Real SNWeakForm::VolumetricScatteringSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Real >* ext) const
{
  Real result(0.0);
  Real *moment_values_at_quad_pts = new Real [n];

  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (unsigned int l = 0; l <= L; l++)
    {
      Real deg_l_result(0.0);
      
      for (int m = -l; m <= l; m++)
      {
        if ( ((l + m) % 2) == 0 )
        {
          odata.ordinates_to_moment(l, m, gfrom, G, ext, n, moment_values_at_quad_pts);
          SupportClasses::SphericalHarmonic Rlm(l, m);
          double rlm_in_dir = Rlm(odata.xi[direction], odata.eta[direction], odata.mu[direction]);
            
          for (int quad_pt = 0; quad_pt < n; quad_pt++)
          {
            Real group_source_moment = moment_values_at_quad_pts[quad_pt] * rlm_in_dir;
            
            // Contribution to spatial integral.
            if (geom_type == HERMES_AXISYM_X)
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt] * e->y[quad_pt];
            else if (geom_type == HERMES_AXISYM_Y)
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt] * e->x[quad_pt];
            else
              deg_l_result += group_source_moment * wt[quad_pt] * v->val[quad_pt];
          }
        }
      }
      
      deg_l_result *= (2*l + 1) / (4*M_PI) * Sigma_sn[l][gto][gfrom];
      result += deg_l_result;
    }
  }
    
  delete [] moment_values_at_quad_pts;
  
  //std::cout << "VolumetricScatteringSourceVF :: " << result << std::endl;
  
  return result;
}


template<typename Real>
Real SNWeakForm::VolumetricFissionSourceVF::vector_form(int n, double* wt, Func< Real >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Real >* ext) const
{
  Real result(0.0);
    
  for (int gfrom = 0; gfrom < G; gfrom++)
  {
    for (int quad_pt = 0; quad_pt < n; quad_pt++)
    {
      // Angular integration to get scalar flux at quadrature point.
      Real group_scalar_flux(0.0);
      for (int dir = 0; dir < odata.M; dir++)
        group_scalar_flux += odata.pw[dir] * ext->fn[ag.pos(dir,gfrom)]->val[quad_pt];
      //group_scalar_flux *= 2 / (4*M_PI) * 2*M_PI;
      group_scalar_flux *= 2;
      
      // Compute fission source with element-wise constant nuSigma_f.
      group_scalar_flux *= nu[gfrom] * Sigma_f[gfrom];
      
      // Contribution to spatial integral.
      if (geom_type == HERMES_AXISYM_X)
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt] * e->y[quad_pt];
      else if (geom_type == HERMES_AXISYM_Y)
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt] * e->x[quad_pt];
      else
        result += group_scalar_flux * wt[quad_pt] * v->val[quad_pt];
    }
  }
  
  return result * chi_to;
}


template<typename Real>
Real SNWeakForm::VolumetricExternalSourceVF::vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, ExtData<Real> *ext) const
{
  Real result = Real(0.0);
  for (int i = 0; i < n; i++)
    result += wt[i] * Q * v->val[i]; //F(e->x[i], e->y[i]) * v->val[i];
  
  //std::cout << "VolumetricExternalSourceVF :: " << result << std::endl;
  
  return result;
}

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
                                                      Geom<double> *e, ExtData<double> *ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
  
    if (a_dot_n < 0)
      if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
      {
//        std::cout << a_dot_n << " --- " << -static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(odata.reflections_about_x[direction], e->nx[quad_pt], e->ny[quad_pt]) << std::endl;
        result += wt[quad_pt] * u->val[quad_pt] * a_dot_n * v->val[quad_pt];      
      }
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_X::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, ExtData<Ord> *ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::SpecularReflectionMF_Y::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, ExtData<double> *ext) const
{
  double eps = 1e-14;
  double result = 0.0;
  //std::cout << "(" << odata.xi[direction] << "," << odata.eta[direction] << ") / ";
  //std::cout << "(" << odata.xi[odata.reflections_about_y[direction]] << "," << odata.eta[odata.reflections_about_y[direction]] << ")" << std::endl;
  for (int quad_pt = 0; quad_pt < n; quad_pt++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
      if (fabs(e->nx[quad_pt] - 1.0) < eps || fabs(e->nx[quad_pt] + 1.0) < eps && fabs(e->ny[quad_pt]) < eps)
      {
//        std::cout << a_dot_n << " | " << -static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(odata.reflections_about_y[direction], e->nx[quad_pt], e->ny[quad_pt]) << std::endl;
        result += wt[quad_pt] * u->val[quad_pt] * a_dot_n * v->val[quad_pt];
      }
  }
  
  return result;
}

Ord SNWeakForm::SpecularReflectionMF_Y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, ExtData<Ord> *ext) const
{ 
  return u->val[0]*v->val[0];
}

double SNWeakForm::InterfaceStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                        Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0.0;

  //std::cout << "IFS " <<  n << std::endl;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);
    //if (std::abs(a_dot_n) < 1e-10)
    //  std::cout << "a_dot_n = " << a_dot_n << std::endl;
    double jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
    result += wt[i] * static_cast<SNWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
  }
  
  //std::cout << "InterfaceStreamingMF :: " << result << std::endl;
  
  return result;
}

double SNWeakForm::BoundaryStreamingMF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                      Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0.0;
  //std::cout << "BSV " <<  n << std::endl;
  
  for (int i = 0; i < n; i++)
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[i], e->ny[i]);
    
    if (a_dot_n >= 0)
      result += wt[i] * u->val[i] * a_dot_n * v->val[i];
  }
  
  //std::cout << "BoundaryStreamingMF :: " << result << std::endl;
  
  return result;
}

Ord SNWeakForm::BoundaryStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, ExtData<Ord> *ext) const
{ 
  return u->val[0]*v->val[0];
}

Ord SNWeakForm::InterfaceStreamingMF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                      Geom<Ord> *e, ExtData<Ord> *ext) const
{ 
  return static_cast<SNWeakForm*>(wf)->upwind_flux(u->get_val_central(0), u->get_val_neighbor(0)) * (v->get_val_central(0) - v->get_val_neighbor(0));
}

double SNWeakForm::BoundaryStreamingVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int quad_pt = 0; quad_pt < n; quad_pt++) 
  {
    double a_dot_n = static_cast<SNWeakForm*>(wf)->calculate_a_dot_v(direction, e->nx[quad_pt], e->ny[quad_pt]);
    
    if (a_dot_n < 0)
    {
      double x = e->x[quad_pt], y = e->y[quad_pt];
      double boundary_data = influx<double>(x,y);
      result += -wt[quad_pt] * boundary_data * a_dot_n * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::BoundaryStreamingVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  //std::cout << "BoundaryStreamingVF :: " << v->val[0] * influx<Ord>(e->x[0], e->y[0]) << std::endl;
  return v->val[0] * influx<Ord>(e->x[0], e->y[0]);
}

double SNWeakForm::SpecularReflectionVF::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, ExtData<double> *ext) const
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
        boundary_data = ext->fn[ag.pos(odata.reflections_about_x[direction],g)]->val[quad_pt];
      else if (fabs(e->ny[quad_pt] - 1.0) < eps || fabs(e->ny[quad_pt] + 1.0) < eps && fabs(e->nx[quad_pt]) < eps)
        boundary_data = ext->fn[ag.pos(odata.reflections_about_y[direction],g)]->val[quad_pt];
      else
        Hermes::Mixins::Loggable::Static::warn("Only horizontal or vertical boundaries are currently supported for specular reflection b. c.");

      result += -wt[quad_pt] * boundary_data * a_dot_n * v->val[quad_pt];
    }
  }
  return result;
}

Ord SNWeakForm::SpecularReflectionVF::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  int refl_ord = std::max(ext->fn[ag.pos(odata.reflections_about_x[direction],g)]->val[0].get_order(),ext->fn[ag.pos(odata.reflections_about_y[direction],g)]->val[0].get_order());
  
  //std::cout << "SpecularReflectionVF :: " << static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), Ord(refl_ord), Ord(1)) * v->val[0] << std::endl;
  
  return static_cast<SNWeakForm*>(wf)->upwind_flux(Ord(0), Ord(refl_ord)) * v->val[0];
}

double SNWeakForm::calculate_a_dot_v(int n, double vx, double vy) const
{
  //return std::cos((n+1)*M_PI/(2.1*N))*vx + std::sin((n+1)*M_PI/(2.1*N))*vy;
  return odata.xi[n]*vx + odata.eta[n]*vy;
}

double SNWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  if (fabs(a_dot_n) < 1e-14)
    return 0;
  else
    return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord SNWeakForm::upwind_flux(Ord u_cent, Ord u_neib) const
{
  return (u_cent + u_neib);
}
