#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 
  using namespace SPN;
    using namespace SupportClasses;

class CustomWeakForm : public WeakForms::FixedSourceProblem
{
  public:
    CustomWeakForm(const MaterialProperties::MaterialPropertyMaps& matprop, unsigned int N);
};

template <typename Scalar>
class ErrorFormSPN : public Adapt<Scalar>::MatrixFormVolError
{
  private:
    double factor;
    
  public:
    
    ErrorFormSPN(unsigned int moment1, unsigned int moment2, ProjNormType norm) : Adapt<Scalar>::MatrixFormVolError(norm)
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

void report_num_dof(const std::string& msg, const Hermes::vector<Space<double> *> spaces);
void report_errors(const std::string& msg, const Hermes::vector< double > errors);
std::string itos(int t);
