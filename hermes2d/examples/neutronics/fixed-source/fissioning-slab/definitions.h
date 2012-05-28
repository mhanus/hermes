#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 

#include "../../utils.h"

// Choose one of the three following models, leave the other ones commented
#define USE_SPN
//#define USE_DIFFUSION_WITH_TRANSPORT_CORRECTION
//#define USE_SIMPLE_DIFFUSION

#ifdef USE_SPN 
  using namespace SPN;
#else // DIFFUSION
  using namespace Diffusion;
#endif

using namespace SupportClasses;

class CustomWeakForm : public WeakForms::FixedSourceProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N);
};

#ifdef USE_SPN
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
#endif
