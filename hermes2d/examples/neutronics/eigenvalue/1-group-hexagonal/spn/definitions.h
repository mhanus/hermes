#include "hermes2d.h"
#include <sstream>

using namespace WeakFormsNeutronics::Multigroup;
using namespace CompleteWeakForms::SPN;
using namespace ElementaryForms::SPN;
using namespace SupportClasses::SPN;

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop, unsigned int N,
                   const Hermes::vector<Solution*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, std::string bdy_vacuum);
};

// Integral over the active core.
double integrate(MeshFunction* sln, Mesh* mesh, const std::vector<std::string>& fission_regions);

void report_num_dof(const std::string& msg, const Hermes::vector<Space *> spaces);


/// \brief Power iteration. 
///
/// Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
/// has converged, also updating the global eigenvalue 'k_eff'.
///
/// \param[in]     hermes2d     Class encapsulating global Hermes2D functions.
/// \param[in]     spaces       Pointers to spaces on which the solutions are defined (one space for each energy group).
/// \param[in]     wf           Pointer to the weak form of the problem.
/// \param[in,out] solution     A set of Solution* pointers to solution components (neutron fluxes in each group). 
///                             Initial guess for the iteration on input, converged result on output.
/// \param[in] fission_regions  Strings specifiying the parts of the solution domain where fission occurs.
/// \param[in]     tol          Relative difference between two successive eigenvalue approximations that stops the iteration.
/// \param[in,out] mat          Pointer to a matrix to which the system associated with the power iteration will be assembled.
/// \param[in,out] rhs          Pointer to a vector to which the right hand sides of the power iteration will be successively assembled.
/// \param[in]     solver       Solver for the resulting matrix problem (specified by \c mat and \c rhs).
///
/// \return  number of iterations needed for convergence within the specified tolerance.
///
int power_iteration(const Hermes2D& hermes2d, const MaterialPropertyMaps& matprop, 
                    const Hermes::vector<Space *>& spaces, DefaultWeakFormSourceIteration* wf, 
                    const Hermes::vector<Solution *>& solution, 
                    const std::vector<std::string>& fission_regions, 
                    double tol, SparseMatrix *mat, Vector* rhs, Solver *solver);
                    
class Views
{
  unsigned int n_moments, n_odd_moments, n_unknowns, n_equations, n_groups;
  bool display_meshes;
  
  ScalarView** sviews;
  OrderView** oviews;
  MeshView** mviews;
  
  static const std::string  base_title_flux;
  static const std::string  base_title_order;
  static const std::string  base_title_mesh;
  static const unsigned int MAX_SOLUTIONS_SETS = 10;
  
  MomentGroupFlattener mg;
  MomentFilter::Val*** moment_filters[MAX_SOLUTIONS_SETS];
  
  std::string itos(int t)
  {
    std::stringstream ss; ss << t;
    return ss.str();
  }
  
  public:
    Views(unsigned int spn_order, unsigned int G, bool display_meshes = false);
    ~Views();
    
    void show_meshes(Hermes::vector<Mesh*> meshes);
    void show_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set = 0);
    void show_orders(Hermes::vector<Space*> spaces);
    
    void inspect_meshes(Hermes::vector<Mesh*> meshes);
    void inspect_solutions(Hermes::vector< Solution* > solutions, unsigned int solutions_set = 0);
    void inspect_orders(Hermes::vector<Space*> spaces);
};