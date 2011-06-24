#define HERMES_REPORT_ALL
#include "definitions.h"
#include "problem_data.h"

using namespace RefinementSelectors;
using namespace WeakFormsNeutronics::Multigroup::MaterialProperties::Diffusion;

const unsigned int N_GROUPS = 1;  // Monoenergetic (single group) problem.

const unsigned int N_EQUATIONS = N_GROUPS;

const int INIT_REF_NUM[N_EQUATIONS] = {  // Initial uniform mesh refinement for the individual solution components.
  2
};
const int P_INIT[N_EQUATIONS] = {        // Initial polynomial orders for the individual solution components. 
  2
}; 
        
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = -1;                 // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  
const bool display_meshes = false;
const bool use_transport_corrected_cross_sections = true;

// Power iteration control.
double k_eff = 1.0;         // Initial eigenvalue approximation.
double TOL_PIT_CM = 1e-7;   // Tolerance for eigenvalue convergence on the coarse mesh.
double TOL_PIT_RM = 1e-8;   // Tolerance for eigenvalue convergence on the fine mesh.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  MaterialPropertyMaps *matprop;
  
  if (use_transport_corrected_cross_sections)
  {
    MaterialPropertyMap2 Ss1;
    MaterialPropertyMap3::const_iterator it = Ssn.begin();
    for ( ; it != Ssn.end(); it++)
      Ss1[it->first] = it->second[1];
    
    matprop = new TransportCorrectedMaterialPropertyMaps(N_GROUPS, Ss1, rm_map);
  }
  else
  {
    matprop = new MaterialPropertyMaps(N_GROUPS, rm_map);
  }
  
  MaterialPropertyMap2 Ss0;
  MaterialPropertyMap3::const_iterator it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss0[it->first] = it->second[0];
  
  matprop->set_nuSigma_f(nSf);
  matprop->set_nu(nu);
  matprop->set_Sigma_t(St);
  matprop->set_Sigma_s(Ss0);
  
  matprop->validate();
  
  cout << *matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  H2DReader mloader;
  mloader.load(mesh_file.c_str(), meshes[0]);
  
  // Convert the mesh so that it has one type of elements (optional). 
  //meshes[0]->convert_quads_to_triangles();
  meshes[0]->convert_triangles_to_quads();
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  
  Views views(N_GROUPS, display_meshes);
  views.inspect_meshes(meshes);

  // Create pointers to solutions on coarse and fine meshes and from the latest power iteration, respectively.
  Hermes::vector<Solution*> coarse_solutions, fine_solutions, power_iterates;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
  {
    coarse_solutions.push_back(new Solution());
    fine_solutions.push_back(new Solution());
    power_iterates.push_back(new Solution(meshes[i], 1.0));   
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space *> spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces.push_back(new H1Space(meshes[i], P_INIT[i]));
  
  // Initialize the weak formulation.
  CustomWeakForm wf(*matprop, power_iterates, fission_regions, k_eff, bdy_vacuum);
  
  // Initialize the discrete algebraic representation of the problem.
  DiscreteProblem dp(&wf, spaces);
  
  // Initialize the eigenvalue iterator.
  SupportClasses::SourceIteration si(NEUTRONICS_DIFFUSION, *matprop, fission_regions, hermes2d, dp);
  
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  report_num_dof("Coarse mesh power iteration, ", spaces);
  si.eigenvalue_iteration(power_iterates, TOL_PIT_CM, matrix_solver);
  
  if (STRATEGY >= 0)
  {
    // DOF and CPU convergence graphs
    GnuplotGraph graph_dof("Error convergence", "NDOF", "log(error)");
    graph_dof.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_dof.add_row("L2 err. est. [%]", "g", "-", "s");
    graph_dof.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_dof.set_log_y();
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Error convergence", "CPU time [s]", "log(error)");
    graph_cpu.add_row("H1 err. est. [%]", "r", "-", "o");
    graph_cpu.add_row("L2 err. est. [%]", "g", "-", "s");
    graph_cpu.add_row("Keff err. est. [milli-%]", "b", "-", "d");
    graph_cpu.set_log_y();
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    Hermes::vector<ProjNormType> proj_norms_h1, proj_norms_l2;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
    {
      proj_norms_h1.push_back(HERMES_H1_NORM);
      proj_norms_l2.push_back(HERMES_L2_NORM);
    }
  }
  else
  {
    views.show_solutions(power_iterates);
    views.show_orders(spaces);
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double keff_err = 1e5*fabs(wf.get_keff() - REF_K_EFF)/REF_K_EFF;
    info("K_eff error = %g pcm", keff_err);
  }
  
  delete matprop;
  
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for the view to be closed.  
  View::wait();
  
  return 0;
}
