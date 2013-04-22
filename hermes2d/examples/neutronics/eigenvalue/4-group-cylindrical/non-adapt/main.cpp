#define HERMES_REPORT_ALL
#include "../weak_formulation.h"
#include "../problem_data.h"

// This example solves a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. The corresponding diffusion operator is given by (r = x, z = y):
//
//	\nabla \cdot D \nabla \phi = \frac{1}{x} (x D \phi_x)_x  + (D \phi_y)_y 
//
// BC:
//
// Homogeneous neumann on symmetry axis,
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.
const int P_INIT_1 = 1,                           // Initial polynomial degree for approximation of group 1 fluxes.
          P_INIT_2 = 1,                           // Initial polynomial degree for approximation of group 2 fluxes.
          P_INIT_3 = 1,                           // Initial polynomial degree for approximation of group 3 fluxes.
          P_INIT_4 = 1;                           // Initial polynomial degree for approximation of group 4 fluxes.
const double ERROR_STOP = 1e-5;                   // Tolerance for the eigenvalue.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)

// Initial eigenvalue approximation.
double k_eff = 1.0;         

int main(int argc, char* argv[])
{  
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 2);

  // Load the mesh.
  MeshSharedPtr mesh = new Mesh();
  MeshReaderH2D mesh_reader;
  mesh_reader.load((std::string("../") + mesh_file).c_str(), mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Solution variables.
  MeshFunctionSharedPtr<double> sln1 = new Solution<double>();
  MeshFunctionSharedPtr<double> sln2 = new Solution<double>();
  MeshFunctionSharedPtr<double> sln3 = new Solution<double>();
  MeshFunctionSharedPtr<double> sln4 = new Solution<double>();
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions(sln1, sln2, sln3, sln4);
  
  // Create H1 spaces with default shapesets.
  SpaceSharedPtr<double> space1 = new H1Space<double>(mesh, P_INIT_1);
  SpaceSharedPtr<double> space2 = new H1Space<double>(mesh, P_INIT_1);
  SpaceSharedPtr<double> space3 = new H1Space<double>(mesh, P_INIT_1);
  SpaceSharedPtr<double> space4 = new H1Space<double>(mesh, P_INIT_1);
  Hermes::vector<SpaceSharedPtr<double> > spaces(space1, space2, space3, space4);
  
  int ndof = Space<double>::get_num_dofs(spaces);
  Loggable::Static::info("ndof = %d.", ndof);
  
  // Initialize views.
  Views::ScalarView view1("Neutron flux 1", new Views::WinGeom(0, 0, 320, 600));
  Views::ScalarView view2("Neutron flux 2", new Views::WinGeom(350, 0, 320, 600));
  Views::ScalarView view3("Neutron flux 3", new Views::WinGeom(700, 0, 320, 600));
  Views::ScalarView view4("Neutron flux 4", new Views::WinGeom(1050, 0, 320, 600));
  
  // Do not show meshes.
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);
  
  // Load physical data of the problem for the 4 energy groups.
  MaterialProperties::MaterialPropertyMaps matprop(4, std::set<std::string>(all_regions.begin(), all_regions.end()));
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_Sigma_a(Sa);
  matprop.set_Sigma_f(Sf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.set_fission_materials(fission_regions);
  matprop.validate();
  
  std::cout << matprop;
  
  // Time measurement.
  TimeMeasurable cpu_time, solver_time;
  cpu_time.tick(); 
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, bdy_vacuum);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);
  
  Neutronics::KeffEigenvalueIteration keff_eigenvalue_iteration(&wf, spaces);
  keff_eigenvalue_iteration.set_keff_tol(ERROR_STOP);
  keff_eigenvalue_iteration.set_max_allowed_iterations(100);

  keff_eigenvalue_iteration.output_matrix(1);
  keff_eigenvalue_iteration.set_matrix_E_matrix_dump_format(DF_HERMES_BIN);
  keff_eigenvalue_iteration.set_matrix_number_format("%1.15f");
  
  keff_eigenvalue_iteration.solve();
  
  Solution<double>::vector_to_solutions(keff_eigenvalue_iteration.get_sln_vector(), spaces, solutions);
  
  // Time measurement.
  cpu_time.tick();
  solver_time.tick(TimeMeasurable::HERMES_SKIP);
  
   
  // Show solutions.
  view1.show(sln1);
  view2.show(sln2);
  view3.show(sln3);    
  view4.show(sln4);
  
  // Skip visualization time.
  cpu_time.tick(TimeMeasurable::HERMES_SKIP);

  // Print timing information.
  Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
    
  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
