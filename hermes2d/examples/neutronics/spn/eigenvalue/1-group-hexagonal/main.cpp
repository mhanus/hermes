#define HERMES_REPORT_ALL
#include "definitions.h"
#include "problem_data.h"

using namespace RefinementSelectors;

const unsigned int N_ODD_MOMENTS = (N_MOMENTS+1)/2;
const unsigned int N_TOTAL = N_GROUPS * N_ODD_MOMENTS;

const int INIT_REF_NUM[N_TOTAL] = {      // Initial uniform mesh refinement for the individual solution components.
  2, 1                              
};
const int P_INIT[N_TOTAL] = {            // Initial polynomial orders for the individual solution components. 
  1, 1                              
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

// Power iteration control.
double k_eff = 1.0;         // Initial eigenvalue approximation.
double TOL_PIT_CM = 5e-5;   // Tolerance for eigenvalue convergence on the coarse mesh.
double TOL_PIT_RM = 5e-6;   // Tolerance for eigenvalue convergence on the fine mesh.

// Macros for simpler reporting (four group case).
#define report_num_dofs(spaces) "%d + %d = %d",\
                                Space::get_num_dofs(spaces[0]), Space::get_num_dofs(spaces[1]),\
                                Space::get_num_dofs(spaces)
#define report_errors(errors) errors[0],errors[1]

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
   
  MaterialPropertyMaps matprop(N_GROUPS, N_MOMENTS, rm_map);
  matprop.set_nuSigma_f(nSf);
  matprop.set_nu(nu);
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
  
  matprop.validate();
  
  cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
  for (unsigned int i = 0; i < N_TOTAL; i++) 
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  H2DReader mloader;
  mloader.load(mesh_file.c_str(), meshes[0]);
  
  // Convert the mesh so that it has one type of elements (optional). 
  //meshes[0].convert_quads_to_triangles();
  //meshes[0].convert_triangles_to_quads();
  
  for (unsigned int i = 1; i < N_TOTAL; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  
  MomentGroupFlattener mg(N_GROUPS);
   
  // Display the meshes.
  if (display_meshes)
  {
    MeshView* mviews[N_TOTAL];  
    std::string base_title = "Core mesh for group ";
    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::string title = base_title + itos(g) + std::string(", moment ");
      for (unsigned int m = 0; m < N_ODD_MOMENTS; m++)
      {
        unsigned int i = mg.pos(m,g);
        mviews[i] = new MeshView((title + itos(m)).c_str(), new WinGeom(m*352, g*352, 350, 350));
        mviews[i]->show(meshes[i]);
      }
    }
    // Wait for the view to be closed.
    View::wait();
    
    for (unsigned int i = 0; i < N_TOTAL; i++)
      delete mviews[i];
  }
  
  ScalarView* sviews[N_TOTAL];
  OrderView* oviews[N_TOTAL];
  std::string base_title_flux = "Neutron flux: group ";
  std::string base_title_order = "Polynomial orders: group ";
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    std::string title_flux = base_title_flux + itos(g) + std::string(", moment ");
    std::string title_order = base_title_order + itos(g) + std::string(", moment ");
    for (unsigned int m = 0; m < N_ODD_MOMENTS; m++)
    {
      unsigned int i = mg.pos(m,g);
      sviews[i] = new ScalarView((title_flux + itos(m)).c_str(), new WinGeom(m*452, g*452, 450, 450));
      oviews[i] = new OrderView((title_order + itos(m)).c_str(), new WinGeom(m*452, N_GROUPS*452 + g*452, 450, 450));
      sviews[i]->show_mesh(false);
      sviews[i]->set_3d_mode(true);
    }
  }
  
  // Create pointers to solutions on coarse and fine meshes and from the latest power iteration, respectively.
  Hermes::vector<Solution*> coarse_solutions, fine_solutions, power_iterates;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_TOTAL; i++) 
  {
    coarse_solutions.push_back(new Solution());
    fine_solutions.push_back(new Solution());
    power_iterates.push_back(new Solution(meshes[i], 1.0));   
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space *> spaces;
  for (unsigned int i = 0; i < N_TOTAL; i++) 
    spaces.push_back(new H1Space(meshes[i], P_INIT[i]));
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, N_MOMENTS, power_iterates, fission_regions, k_eff, bdy_vacuum);
  
  // Initialize the discrete algebraic representation of the problem and its solver.
  //
  // Create the matrix and right-hand side vector for the solver.
  SparseMatrix* mat = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  // Instantiate the solver itself.
  Solver* solver = create_linear_solver(matrix_solver, mat, rhs);
  
  // Initial power iteration to obtain a coarse estimate of the eigenvalue and the fission source.
  info("Coarse mesh power iteration, " report_num_dofs(spaces));
  power_iteration(hermes2d, matprop, spaces, &wf, power_iterates, fission_regions, TOL_PIT_CM, mat, rhs, solver);
  
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
    for (unsigned int i = 0; i< N_TOTAL; i++)
    {
      proj_norms_h1.push_back(HERMES_H1_NORM);
      proj_norms_l2.push_back(HERMES_L2_NORM);
    }
  }
  else
  {
    for (unsigned int i = 0; i < N_TOTAL; i++)
    {
      sviews[i]->show(power_iterates[i]);
      oviews[i]->show(spaces[i]);
    }
  }
  
  // Wait for the view to be closed.
  View::wait();
  return 0;
}
