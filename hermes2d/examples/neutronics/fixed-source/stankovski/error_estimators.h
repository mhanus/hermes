#include "hermes2d.h"
using namespace Hermes::Hermes2D; 

#include "neutronics/support_classes.h"
using namespace Neutronics;
using namespace SPN::SupportClasses;

///////////////////////////////////// SPN /////////////////////////////////////

///////////////////////// PROJECTION BASED ADAPTIVITY 

template <typename Scalar>
class ErrorFormSPN : public Adapt<Scalar>::MatrixFormVolError
{
  private:
    double factor;
    
  public:
    
    ErrorFormSPN(unsigned int moment1, unsigned int moment2, ProjNormType norm) : Adapt<Scalar>::MatrixFormVolError(0,0,norm)
    { 
      factor = Coeffs::even_moment(0, moment1) * Coeffs::even_moment(0, moment2);
    }

    /// Evaluate value of the error norm.
    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
                          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
                          ExtData<Scalar> *ext) const
    {
      return factor * Adapt<Scalar>::MatrixFormVolError::value(n, wt, u_ext, u, v, e, ext);
    }
    
    /// Use the default form from Adapt<Scalar>::MatrixFormVolError to evaluate the quadrature order.
};

//////////////////// ADAPTIVITY BASED ON JUMPS OF SOLUTION 

class InterfaceEstimatorSPN : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
  
public:
  
  InterfaceEstimatorSPN(unsigned int m, unsigned int gto, const MomentGroupFlattener& mg,
                        const SPN::MaterialProperties::MaterialPropertyMaps& matprop, bool multiply_by_scalar_flux_factor) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(mg.pos(m,gto), Hermes::H2D_DG_INNER_EDGE),
      mrow(m), gto(gto), mg(mg), matprop(matprop), 
      scalar_flux_factor(multiply_by_scalar_flux_factor ? Coeffs::even_moment(0, m) : 1.0)
  { };

  /// Evaluate value of the error norm.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
  {
    double result = 0.;
    
    std::string mat = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker));
    std::string mat_neighbor = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->get_neighbor_marker()));
    rank1 odd_Sigma_rn_inv = matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto];
    rank1 odd_Sigma_rn_inv_neighbor = matprop.get_odd_Sigma_rn_inv(mat_neighbor)[mrow][gto];
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {
      double coeff = Coeffs::D(mrow) * (odd_Sigma_rn_inv[gfrom] + odd_Sigma_rn_inv_neighbor[gfrom]) * 0.5;
      unsigned int sln_idx = mg.pos(mrow, gfrom);
      
      for (int i = 0; i < n; i++)
        result += wt[i] * Hermes::sqr( coeff * ( e->nx[i] * (u_ext[sln_idx]->get_dx_central(i) - u_ext[sln_idx]->get_dx_neighbor(i)) +
                                                 e->ny[i] * (u_ext[sln_idx]->get_dy_central(i) - u_ext[sln_idx]->get_dy_neighbor(i)) ) );
    }
    
    return result * scalar_flux_factor;
  }
  
  /// Evaluate quadrature order required for accurate integration of the error norm.
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                          Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                          ExtData<Hermes::Ord> *ext) const
  {
    Hermes::Ord result(0);
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {
      unsigned int sln_idx = mg.pos(mrow, gfrom);
      result += Hermes::sqr( (u_ext[sln_idx]->get_dx_central(0) - u_ext[sln_idx]->get_dx_neighbor(0)) +
                             (u_ext[sln_idx]->get_dy_central(0) - u_ext[sln_idx]->get_dy_neighbor(0)) );
    }
    
    return result;
  }

private:
  
    unsigned int mrow;
    unsigned int gto;
    MomentGroupFlattener mg;
    const SPN::MaterialProperties::MaterialPropertyMaps& matprop;
    double scalar_flux_factor;
};

class BoundaryEstimatorSPN : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
  
public:
  
  BoundaryEstimatorSPN(unsigned int m, unsigned int gto, const MomentGroupFlattener& mg,
                        const SPN::MaterialProperties::MaterialPropertyMaps& matprop, bool multiply_by_scalar_flux_factor) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(mg.pos(m,gto), Hermes::HERMES_ANY),
      mrow(m), gto(gto), mg(mg), matprop(matprop),
      scalar_flux_factor(multiply_by_scalar_flux_factor ? Coeffs::even_moment(0, m) : 1.0)
  { };

  /// Evaluate value of the error norm.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
  {
    double result = 0.;
    
    std::string mat = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker));
    rank1 odd_Sigma_rn_inv = matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto];
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {
      double coeff = Coeffs::D(mrow) * odd_Sigma_rn_inv[gfrom];
      unsigned int sln_idx = mg.pos(mrow, gfrom);
      
      for (int i = 0; i < n; i++)
        result += wt[i] * Hermes::sqr( coeff * ( e->nx[i] * u_ext[sln_idx]->dx[i] +
                                                 e->ny[i] * u_ext[sln_idx]->dy[i] ) );
    }
    
    return result * e->diam * scalar_flux_factor;
  }
  
  /// Evaluate quadrature order required for accurate integration of the error norm.
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                          Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                          ExtData<Hermes::Ord> *ext) const
  {
    Hermes::Ord result(0);
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {
      unsigned int sln_idx = mg.pos(mrow, gfrom);
      result += Hermes::sqr( e->nx[0] * u_ext[sln_idx]->dx[0] +
                             e->ny[0] * u_ext[sln_idx]->dy[0] );
    }
    
    return result;
  }

private:
  
    unsigned int mrow;
    unsigned int gto;
    MomentGroupFlattener mg;
    const SPN::MaterialProperties::MaterialPropertyMaps& matprop;
    double scalar_flux_factor;
};

class VolumetricEstimatorSPN : public KellyTypeAdapt<double>::ErrorEstimatorForm
{  
public:
  
  VolumetricEstimatorSPN(unsigned int m, unsigned int gto, const MomentGroupFlattener& mg,
                        const SPN::MaterialProperties::MaterialPropertyMaps& matprop, bool multiply_by_scalar_flux_factor) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(mg.pos(m,gto), Hermes::HERMES_ANY),
      mrow(m), gto(gto), mg(mg), matprop(matprop),
      scalar_flux_factor(multiply_by_scalar_flux_factor ? Coeffs::even_moment(0, m) : 1.0)
  { 
    chi_nnz = matprop.get_fission_nonzero_structure();
  };

  /// Evaluate value of the error norm.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
  {
    std::string mat = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker));
    
    rank1 odd_Sigma_rn_inv = matprop.get_odd_Sigma_rn_inv(mat)[mrow][gto];
    
    rank1 nu_elem = matprop.get_nu(mat);
    rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
    double chi = matprop.get_chi(mat)[gto];
    if (chi_nnz[gto])
    {
      nu_elem = matprop.get_nu(mat);
      Sigma_f_elem = matprop.get_Sigma_f(mat);
      chi = matprop.get_chi(mat)[gto];
    }
    
    double *result = new double [n];
    std::fill(result, result+n, Coeffs::even_moment(0, mrow) * matprop.get_iso_src(mat)[gto]);
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {
      double nSf = nu_elem[gfrom] * Sigma_f_elem[gfrom];
      double D = Coeffs::D(mrow) * odd_Sigma_rn_inv[gfrom];
      
      for (unsigned int mcol = 0; mcol < matprop.get_N_odd(); mcol++)
      {
        double Sigma_rn_elem = 0.;
        for (unsigned int k = 0; k <= mrow; k++)
          Sigma_rn_elem += Coeffs::system_matrix(mrow, mcol, k) * matprop.get_Sigma_rn(mat)[2*k][gto][gfrom];
              
        unsigned int sln_idx = mg.pos(mcol, gfrom);
      
        for (int l = 0; l < n; l++)
        {          
          if (mrow == mcol)
            result[l] += D * u_ext[sln_idx]->laplace[l];
          if (chi_nnz[gto])
            result[l] += -chi * nSf * Coeffs::system_matrix(mrow, mcol, 0) * u_ext[sln_idx]->val[l];
          result[l] += Sigma_rn_elem * u_ext[sln_idx]->val[l];
        }
      }
    }
    
    double final_result = 0.0;
    for (int l = 0; l < n; l++)
      final_result += wt[l] * Hermes::sqr(result[l]);
    
    delete [] result;
    
    return final_result * Hermes::sqr(e->diam) * scalar_flux_factor;
  }
  
  /// Evaluate quadrature order required for accurate integration of the error norm.
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                          Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                          ExtData<Hermes::Ord> *ext) const
  {
    Hermes::Ord result(0);
    
    for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
    {      
      for (unsigned int mcol = 0; mcol < matprop.get_N_odd(); mcol++)
      {
        unsigned int sln_idx = mg.pos(mcol, gfrom);
        result += Hermes::sqr( u_ext[sln_idx]->val[0] + u_ext[sln_idx]->laplace[0] );
      }
    }
    
    return result;
  }

private:
  
    unsigned int mrow;
    unsigned int gto;
    MomentGroupFlattener mg;
    const SPN::MaterialProperties::MaterialPropertyMaps& matprop;
    double scalar_flux_factor;
    bool1 chi_nnz;
};

///////////////////////////////////// DIFFUSION /////////////////////////////////////

////////////////////// ADAPTIVITY BASED ON JUMPS OF SOLUTION 

class InterfaceEstimatorDiffusion : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  
  InterfaceEstimatorDiffusion(unsigned int gto, const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(gto, Hermes::H2D_DG_INNER_EDGE),
      gto(gto), matprop(matprop)
  { };

  /// Evaluate value of the error norm.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
  {
    double result = 0.;
    
    std::string mat = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker));
    std::string mat_neighbor = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->get_neighbor_marker()));
    double D = matprop.get_D(mat)[gto];
    double D_neighbor = matprop.get_D(mat_neighbor)[gto];
    double coeff = (D + D_neighbor) * 0.5;
    
    for (int i = 0; i < n; i++)
      result += wt[i] * Hermes::sqr( coeff * ( e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                                               e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i)) ) );
    
    return result;
  }
  
  /// Evaluate quadrature order required for accurate integration of the error norm.
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                          Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                          ExtData<Hermes::Ord> *ext) const
  {
    return Hermes::sqr( (u->get_dx_central(0) - u->get_dx_neighbor(0)) +
                        (u->get_dy_central(0) - u->get_dy_neighbor(0)) );
  }

private:
  
    unsigned int gto;
    const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop;
};

class BoundaryEstimatorDiffusion : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
  
public:
  
  BoundaryEstimatorDiffusion(unsigned int gto, const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(gto, Hermes::HERMES_ANY),
      gto(gto), matprop(matprop)
  { };

  /// Evaluate value of the error norm.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
  {
    double result = 0.;
    
    std::string mat = matprop.get_material(adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker));
    double D = matprop.get_D(mat)[gto];
    
    for (int i = 0; i < n; i++)
      result += wt[i] * Hermes::sqr( D * ( e->nx[i] * u->dx[i] + e->ny[i] * u->dy[i] ) );
    
    return result * e->diam;
  }
  
  /// Evaluate quadrature order required for accurate integration of the error norm.
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                          Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                          ExtData<Hermes::Ord> *ext) const
  {
    return Hermes::sqr( e->nx[0] * u->dx[0] + e->ny[0] * u->dy[0] );
  }

private:
  
    unsigned int gto;
    const Diffusion::MaterialProperties::MaterialPropertyMaps& matprop;
};
