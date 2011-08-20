#ifndef ___NEUTRONICS_WEAK_FORM_PARTS_IMPLEMENTATION_H_
#define ___NEUTRONICS_WEAK_FORM_PARTS_IMPLEMENTATION_H_

#include "neutronics/weakform_parts.h"

namespace Hermes { namespace Hermes2D { namespace Neutronics 
{  
  namespace Diffusion { namespace WeakFormParts
  { 
    template<typename Real>
    Real VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result;
      
      if (geom_type == HERMES_PLANAR) 
        result = 0.5 * int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = 0.5 * int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = 0.5 * int_x_u_v<Real, Real>(n, wt, u, v, e);
      
      return result;
    }
    
    template<typename Real>
    Real VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result;
      
      if (geom_type == HERMES_PLANAR) 
        result = 0.5 * int_u_ext_v<Real, Real>(n, wt, u_ext[g], v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = 0.5 * int_y_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
      else 
        result = 0.5 * int_x_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
      
      return result;
    }
     
    template<typename Real>
    Real DiffusionReaction::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    {
      Real result;
      
      std::string mat = matprop.get_material(e->elem_marker, wf);     
      rank1 D_elem = matprop.get_D(mat);
      rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
      
      if (geom_type == HERMES_PLANAR) 
      {
        result = D_elem[g] * int_grad_u_grad_v<Real, Real>(n, wt, u, v) +
                  Sigma_r_elem[g] * int_u_v<Real, Real>(n, wt, u, v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D_elem[g] * int_y_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                    Sigma_r_elem[g] * int_y_u_v<Real, Real>(n, wt, u, v, e);
        }
        else 
        {
          result = D_elem[g] * int_x_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                    Sigma_r_elem[g] * int_x_u_v<Real, Real>(n, wt, u, v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real DiffusionReaction::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result;
      
      std::string mat = matprop.get_material(e->elem_marker, wf);        
      rank1 D_elem = matprop.get_D(mat);
      rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
      
      if (geom_type == HERMES_PLANAR) 
      {
        result = D_elem[g] * int_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v) +
                  Sigma_r_elem[g] * int_u_ext_v<Real, Real>(n, wt, u_ext[g], v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D_elem[g] * int_y_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v, e) + 
                    Sigma_r_elem[g] * int_y_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
        }
        else 
        {
          result = D_elem[g] * int_x_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v, e) + 
                    Sigma_r_elem[g] * int_x_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const 
    {
      if (!matprop.get_fission_nonzero_structure()[gto])
        return Real(0);
      
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Real>(n, wt, u, v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Real>(n, wt, u, v, e);
        else result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      }
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
      
      return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
    }
    
    template<typename Real>
    Real FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      if (!matprop.get_fission_nonzero_structure()[g])
        return Real(0);
        
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
      
      if ((unsigned)ext->nf != nu_elem.size() || (unsigned)ext->nf != Sigma_f_elem.size())
        error_function(Messages::E_INVALID_GROUP_INDEX);
      
      Real result(0);
      for (int i = 0; i < n; i++) 
      {
        Real local_res(0);
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
    
    template<typename Real>
    Real FissionYield::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      if (!matprop.get_fission_nonzero_structure()[gto])
        return Real(0);
      
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
        else result = int_x_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
      }
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
      
      return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
    }
                                                    
    template<typename Real>
    Real Scattering::Jacobian::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const  
    {
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Real>(n, wt, u, v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Real>(n, wt, u, v, e);
        else result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      }
      
      return result * matprop.get_Sigma_s(matprop.get_material(e->elem_marker, wf))[gto][gfrom];
    }
    
    template<typename Real>
    Real Scattering::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                              Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
        else result = int_x_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
      }
      
      return result * matprop.get_Sigma_s(matprop.get_material(e->elem_marker, wf))[gto][gfrom];
    }    
   
    template<typename Real>
    Real ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
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
    
  /* WeakFormParts */
  }
  /* Diffusion */
  }
  
  namespace SPN { namespace WeakFormParts
  {        
    template<typename Real>
    Real VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result;
      
      if (geom_type == HERMES_PLANAR) 
        result = int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      
      return Coeffs::D_grad_F(mrow, mcol) * result;
    }
    
    template<typename Real>
    Real VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result(0);
      
      for (unsigned int mcol = 0; mcol < N_odd; mcol++)
      {
        double coeff = Coeffs::D_grad_F(mrow, mcol);

        unsigned int i = mg.pos(mcol,g);
        
        if (geom_type == HERMES_PLANAR) 
          result += coeff * int_u_ext_v<Real, Real>(n, wt, u_ext[i], v);
        else if (geom_type == HERMES_AXISYM_X) 
          result += coeff * int_y_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
        else 
          result += coeff * int_x_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
      }
      
      return result;
    }
                                                        
    template<typename Real>
    Real DiagonalStreamingAndReactions::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const
    {
      Real result;
      
      std::string mat = matprop.get_material(e->elem_marker, wf);     
      
      double Sigma_r_elem = 0.;
      for (unsigned int k = 0; k <= mrow; k++)
        Sigma_r_elem += Coeffs::system_matrix(mrow, mrow, k) * matprop.get_Sigma_rn(mat)[2*k][g][g];
      
      double D_elem = -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][g][g];
      
      // cout << "DiagonalStreamingAndReactions::Jacobian (mom. #" << mrow << ") | " << mat << " | Sigma_r = " << Sigma_r_elem << " | D = " << D_elem << endl;

      if (geom_type == HERMES_PLANAR) 
      {
        result = D_elem * int_grad_u_grad_v<Real, Real>(n, wt, u, v) +
                  Sigma_r_elem * int_u_v<Real, Real>(n, wt, u, v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D_elem * int_y_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                    Sigma_r_elem * int_y_u_v<Real, Real>(n, wt, u, v, e);
        }
        else 
        {
          result = D_elem * int_x_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                    Sigma_r_elem * int_x_u_v<Real, Real>(n, wt, u, v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real DiagonalStreamingAndReactions::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const
    { 
      Real result;
      
      std::string mat = matprop.get_material(e->elem_marker, wf);     
      
      double Sigma_r_elem = 0.;
      for (unsigned int k = 0; k <= mrow; k++)
        Sigma_r_elem += Coeffs::system_matrix(mrow, mrow, k) * matprop.get_Sigma_rn(mat)[2*k][g][g];
      
      double D_elem = -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][g][g];
      
      // cout << "DiagonalStreamingAndReactions::Residual (mom. #" << mrow << ") | " << mat << " | Sigma_r = " << Sigma_r_elem << " | D = " << D_elem << endl;          
      
      unsigned int i = mg.pos(mrow,g);
      
      if (geom_type == HERMES_PLANAR) 
        result = D_elem * int_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v) +
                  Sigma_r_elem * int_u_ext_v<Real, Real>(n, wt, u_ext[i], v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = D_elem * int_y_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e) + 
                  Sigma_r_elem * int_y_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
      else 
        result = D_elem * int_x_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e) + 
                  Sigma_r_elem * int_x_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
      
      return result;
    }
    
    template<typename Real>
    Real FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const 
    {
      if (!matprop.get_fission_nonzero_structure()[gto])
        return Real(0);
      
      Real result;
      
      if (geom_type == HERMES_PLANAR) 
        result = int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
      
      return result * (-Coeffs::system_matrix(mrow, mcol, 0)) * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
    }
    
    template<typename Real>
    Real FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {  
      if (!matprop.get_fission_nonzero_structure()[g])
        return Real(0);
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
                
      Real result(0);
      for (int i = 0; i < n; i++) 
      {
        Real local_res(0);
        for (int gfrom = 0; gfrom < ext->nf; gfrom++)
          local_res += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext->fn[gfrom]->val[i]; // double flux in group 'gfrom'
        
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
    
    template<typename Real>
    Real FissionYield::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      if (!matprop.get_fission_nonzero_structure()[gto])
        return Real(0);
                
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 nu_elem = matprop.get_nu(mat);
      rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
      rank1 chi_elem = matprop.get_chi(mat);
      
      Real result(0);
      for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
      {
        double nSf = nu_elem[gfrom] * Sigma_f_elem[gfrom];
        
        for (unsigned int mcol = 0; mcol < N_odd; mcol++)
        {
          if (geom_type == HERMES_PLANAR) 
            result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_u_ext_v<Real, Real>(n, wt, u_ext[mg.pos(mcol,gfrom)], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_y_u_ext_v<Real, Real>(n, wt, u_ext[mg.pos(mcol,gfrom)], v, e);
          else 
            result += nSf * (-Coeffs::system_matrix(mrow, mcol, 0)) * int_x_u_ext_v<Real, Real>(n, wt, u_ext[mg.pos(mcol,gfrom)], v, e);
        }
      }
      
      return result * chi_elem[gto];
    }
    
    template<typename Real>
    Real OffDiagonalStreaming::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const
    {
      if (gfrom == gto)
        return Real(0);
      
      Real result(0);
      
      if (geom_type == HERMES_PLANAR) 
        result = int_grad_u_grad_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_grad_u_grad_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_grad_u_grad_v<Real, Real>(n, wt, u, v, e);
      
      std::string mat = matprop.get_material(e->elem_marker, wf);     
      
      // cout << "OffDiagonalStreaming::Jacobian (mom. #" << mrow << ") | " << mat << " | D = " << -Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto][gfrom] << endl;          
      
      return -result * Coeffs::D(mrow) * matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto][gfrom];
    }
    
    template<typename Real>
    Real OffDiagonalStreaming::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      Real result(0);
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank1 D_elem = matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto];
      
      for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
      { 
        if (gfrom != gto)
        {
          unsigned int i = mg.pos(mrow, gfrom);
          
          if (geom_type == HERMES_PLANAR) 
            result += D_elem[gfrom] * int_u_ext_v<Real, Real>(n, wt, u_ext[i], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result += D_elem[gfrom] * int_y_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
          else 
            result += D_elem[gfrom] * int_x_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
        }
      }
      
      // cout << "OffDiagonalStreaming::Residual (mom. #" << mrow << ") | " << mat << " | D = " << -Coeffs::D(mrow) * D_elem[0] << endl;          
      
      return -result * Coeffs::D(mrow);
    }
    
    template<typename Real>
    Real OffDiagonalReactions::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const
    {
      if (mrow == mcol)
        return Real(0);
      
      Real result(0);
                                    
      if (geom_type == HERMES_PLANAR)
        result = int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      
      std::string mat = matprop.get_material(e->elem_marker, wf);
      
      double Sigma_rn_elem = 0.;
      for (unsigned int k = 0; k <= mrow; k++)
        Sigma_rn_elem += Coeffs::system_matrix(mrow, mcol, k) * matprop.get_Sigma_rn(mat)[2*k][gto][gfrom];
      
      // cout << "OffDiagonalReactions::Jacobian (mom. #(" << mrow << "," << mcol << ") | " << mat << " | Sigma_r = " << Sigma_rn_elem << endl;
      
      return result * Sigma_rn_elem;
    }
    
    template<typename Real>
    Real OffDiagonalReactions::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      std::string mat = matprop.get_material(e->elem_marker, wf);
      rank3 Sigma_rn_elem = matprop.get_Sigma_rn(mat);
      
      Real result(0);
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
              result += coeff * int_u_ext_v<Real, Real>(n, wt, u_ext[j], v);
            else if (geom_type == HERMES_AXISYM_X) 
              result += coeff * int_y_u_ext_v<Real, Real>(n, wt, u_ext[j], v, e);
            else 
              result += coeff * int_x_u_ext_v<Real, Real>(n, wt, u_ext[j], v, e);
          }
        }
      }
      
      return result;
    }
    
    template<typename Real>
    Real ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      std::string mat = matprop.get_material(e->elem_marker, wf);
      
      if (geom_type == HERMES_PLANAR) 
        return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_v<Real>(n, wt, v);
      else if (geom_type == HERMES_AXISYM_X) 
        return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_y_v<Real>(n, wt, v, e);
      else 
        return Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[g] * int_x_v<Real>(n, wt, v, e);
    }
  
  /* WeakFormParts */
  }
  /* SPN */
  }
    
/* Neutronics */
}
/* Hermes2D */
}
/* Hermes */
}

#endif
