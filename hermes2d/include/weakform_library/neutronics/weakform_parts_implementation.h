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
            
      if (geom_type == HERMES_PLANAR) 
      {
        result = D * int_grad_u_grad_v<Real, Real>(n, wt, u, v) +
                 Sigma_r * int_u_v<Real, Real>(n, wt, u, v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D * int_y_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                   Sigma_r * int_y_u_v<Real, Real>(n, wt, u, v, e);
        }
        else 
        {
          result = D * int_x_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                   Sigma_r * int_x_u_v<Real, Real>(n, wt, u, v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real DiffusionReaction::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    { 
      Real result;
            
      if (geom_type == HERMES_PLANAR) 
      {
        result = D * int_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v) +
                 Sigma_r * int_u_ext_v<Real, Real>(n, wt, u_ext[g], v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D * int_y_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v, e) + 
                   Sigma_r * int_y_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
        }
        else 
        {
          result = D * int_x_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[g], v, e) + 
                   Sigma_r * int_x_u_ext_v<Real, Real>(n, wt, u_ext[g], v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const 
    { 
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Real>(n, wt, u, v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Real>(n, wt, u, v, e);
        else result = int_x_u_v<Real, Real>(n, wt, u, v, e);
      }
            
      return result * chi_to * nu_from * Sigma_f_from;
    }
    
    template<typename Real>
    Real FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {         
      
      if ((unsigned)ext->get_nf() != nu.size() || (unsigned)ext->get_nf() != Sigma_f.size())
        ErrorHandling::error_function(Messages::E_INVALID_GROUP_INDEX);
      
      Real result(0);
      for (int i = 0; i < n; i++) 
      {
        Real local_res(0);
        for (int gfrom = 0; gfrom < ext->get_nf(); gfrom++)
          local_res += nu[gfrom] * Sigma_f[gfrom] * ext->fn[gfrom]->val[i];
                  
        local_res = local_res * wt[i] * v->val[i];
        
        if (geom_type == HERMES_AXISYM_X)
          local_res = local_res * e->y[i];
        else if (geom_type == HERMES_AXISYM_Y)
          local_res = local_res * e->x[i];
        
        result += local_res;
      }
      
      return result * chi_to / keff;
    }
    
    template<typename Real>
    Real FissionYield::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {       
      Real result(0);
      if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
        else result = int_x_u_ext_v<Real, Real>(n, wt, u_ext[gfrom], v, e);
      }
            
      return result * chi_to * nu_from * Sigma_f_from;
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
      
      return result * Sigma_s_to_from;
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
      
      return result * Sigma_s_to_from;
    }    
   
    template<typename Real>
    Real ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const 
    {       
      if (geom_type == HERMES_PLANAR) 
        return src * int_v<Real>(n, wt, v);
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
          return src * int_y_v<Real>(n, wt, v, e);
        else 
          return src * int_x_v<Real>(n, wt, v, e);
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

      if (geom_type == HERMES_PLANAR) 
      {
        result = D * int_grad_u_grad_v<Real, Real>(n, wt, u, v) +
                 Sigma_r * int_u_v<Real, Real>(n, wt, u, v);
      }
      else 
      {
        if (geom_type == HERMES_AXISYM_X) 
        {
          result = D * int_y_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                   Sigma_r * int_y_u_v<Real, Real>(n, wt, u, v, e);
        }
        else 
        {
          result = D * int_x_grad_u_grad_v<Real, Real>(n, wt, u, v, e) + 
                   Sigma_r * int_x_u_v<Real, Real>(n, wt, u, v, e);
        }
      }
      return result;
    }
    
    template<typename Real>
    Real DiagonalStreamingAndReactions::Residual::vector_form(int n, double *wt, Func<Real> *u_ext[],
                                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext) const
    { 
      Real result;
      
      unsigned int i = mg.pos(mrow,g);
      
      if (geom_type == HERMES_PLANAR) 
        result = D * int_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v) +
                 Sigma_r * int_u_ext_v<Real, Real>(n, wt, u_ext[i], v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = D * int_y_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e) + 
                 Sigma_r * int_y_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
      else 
        result = D * int_x_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e) + 
                 Sigma_r * int_x_u_ext_v<Real, Real>(n, wt, u_ext[i], v, e);
      
      return result;
    }
    
    template<typename Real>
    Real FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const 
    {      
      Real result;
      
      if (geom_type == HERMES_PLANAR) 
        result = int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_u_v<Real, Real>(n, wt, u, v, e);
            
      return result * (-Coeffs::system_matrix(mrow, mcol, 0)) * chi_to * nu_from * Sigma_f_from;
    }
    
    template<typename Real>
    Real FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                          Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {                 
      Real result(0);
      for (int i = 0; i < n; i++) 
      {
        Real local_res(0);
        for (int gfrom = 0; gfrom < ext->get_nf(); gfrom++)
          local_res += nu[gfrom] * Sigma_f[gfrom] * ext->fn[gfrom]->val[i]; // scalar flux in group 'gfrom'
                
        local_res = local_res * wt[i] * v->val[i];
        
        if (geom_type == HERMES_AXISYM_X)
          local_res = local_res * e->y[i];
        else if (geom_type == HERMES_AXISYM_Y)
          local_res = local_res * e->x[i];
        
        result += local_res;
      }
      
      return result * Coeffs::even_moment(0, mrow) * chi_to / keff;
    }
    
    template<typename Real>
    Real FissionYield::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {       
      Real result(0);
      for (unsigned int gfrom = 0; gfrom < G; gfrom++)
      {
        double nSf = nu[gfrom] * Sigma_f[gfrom];
        
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
      
      return result * chi_to;
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
  
      return result * D;
    }
    
    template<typename Real>
    Real OffDiagonalStreaming::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    { 
      Real result(0);

      for (unsigned int gfrom = 0; gfrom < G; gfrom++)
      { 
        if (gfrom != gto)
        {
          unsigned int i = mg.pos(mrow, gfrom);
          
          if (geom_type == HERMES_PLANAR) 
            result += D[gfrom] * int_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result += D[gfrom] * int_y_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e);
          else 
            result += D[gfrom] * int_x_grad_u_ext_grad_v<Real, Real>(n, wt, u_ext[i], v, e);
        }
      }
            
      return -result * Coeffs::D(mrow);
    }
    
    template<typename Real>
    Real OffDiagonalReactions::Jacobian::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext  ) const
    {
      if (mg.pos(mrow, gto) == mg.pos(mcol, gfrom))
        return Real(0);
      
      Real result(0);
                                    
      if (geom_type == HERMES_PLANAR)
        result = int_u_v<Real, Real>(n, wt, u, v);
      else if (geom_type == HERMES_AXISYM_X) 
        result = int_y_u_v<Real, Real>(n, wt, u, v, e);
      else 
        result = int_x_u_v<Real, Real>(n, wt, u, v, e);
            
      return result * Sigma_rn;
    }
    
    template<typename Real>
    Real OffDiagonalReactions::Residual::vector_form( int n, double *wt, Func<Real> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Real> *ext ) const 
    {       
      Real result(0);
      unsigned int i = mg.pos(mrow, gto);
      for (unsigned int gfrom = 0; gfrom < G; gfrom++)
      {            
        for (unsigned int mcol = 0; mcol < N_odd; mcol++)
        {
          unsigned int j = mg.pos(mcol, gfrom);
          
          if (i != j)
          {
            double coeff = 0.;
            for (unsigned int k = 0; k <= std::min(mrow, mcol); k++)
              coeff += Sigma_rn[2*k][gto][gfrom] * Coeffs::system_matrix(mrow, mcol, k);
            
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
      if (geom_type == HERMES_PLANAR) 
        return Coeffs::even_moment(0, mrow) * src * int_v<Real>(n, wt, v);
      else if (geom_type == HERMES_AXISYM_X) 
        return Coeffs::even_moment(0, mrow) * src * int_y_v<Real>(n, wt, v, e);
      else 
        return Coeffs::even_moment(0, mrow) * src * int_x_v<Real>(n, wt, v, e);
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
