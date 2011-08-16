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

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT_1 = 1,                           // Initial polynomial degree for approximation of group 1 fluxes.
          P_INIT_2 = 1,                           // Initial polynomial degree for approximation of group 2 fluxes.
          P_INIT_3 = 2,                           // Initial polynomial degree for approximation of group 3 fluxes.
          P_INIT_4 = 2;                           // Initial polynomial degree for approximation of group 4 fluxes.
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
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mesh_reader;
  mesh_reader.load((std::string("../") + mesh_file).c_str(), &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Solution variables.
  Solution<double> sln1, sln2, sln3, sln4;
  Hermes::vector<Solution<double>*> solutions(&sln1, &sln2, &sln3, &sln4);
  
  // Define initial conditions.
  info("Setting initial conditions.");
  ConstantSolution<double> iter1(&mesh, 1.0), 
                           iter2(&mesh, 1.0), 
                           iter3(&mesh, 1.0), 
                           iter4(&mesh, 1.0);
  Hermes::vector<Solution<double>*> iterates(&iter1, &iter2, &iter3, &iter4);

  // Create H1 spaces with default shapesets.
  H1Space<double> space1(&mesh, P_INIT_1);
  H1Space<double> space2(&mesh, P_INIT_2);
  H1Space<double> space3(&mesh, P_INIT_3);
  H1Space<double> space4(&mesh, P_INIT_4);
  Hermes::vector<Space<double>*> spaces(&space1, &space2, &space3, &space4);
  
  int ndof = Space<double>::get_num_dofs(spaces);
  info("ndof = %d.", ndof);
  
  // Initialize views.
  Views::ScalarView<double> view1("Neutron flux 1", new Views::WinGeom(0, 0, 320, 600));
  Views::ScalarView<double> view2("Neutron flux 2", new Views::WinGeom(350, 0, 320, 600));
  Views::ScalarView<double> view3("Neutron flux 3", new Views::WinGeom(700, 0, 320, 600));
  Views::ScalarView<double> view4("Neutron flux 4", new Views::WinGeom(1050, 0, 320, 600));
  
  // Do not show meshes.
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);
  
  // Load physical data of the problem for the 4 energy groups.
  MaterialProperties::MaterialPropertyMaps matprop(4);
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_Sigma_a(Sa);
  matprop.set_Sigma_f(Sf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.validate();
  
  std::cout << matprop;
  
  // Time measurement.
  Hermes::TimePeriod cpu_time, solver_time;
  cpu_time.tick(); 
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, iterates, k_eff, bdy_vacuum);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);
  
  // Initial coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];
  
  NewtonSolver<double> solver(&dp, matrix_solver);
  solver.set_verbose_output(false);
  
  if (matrix_solver == Hermes::SOLVER_AZTECOO) 
  {
    solver.set_iterative_method(iterative_method);
    solver.set_preconditioner(preconditioner);
  }
  
  // Main power iteration loop:
  int it = 1; bool done = false;
  do
  {
    solver_time.tick();
    
    info("------------ Power iteration %d:", it);
    
    info("Newton's method (matrix problem solved by %s).", Hermes::MatrixSolverNames[matrix_solver].c_str());
    
    memset(coeff_vec, 0, ndof*sizeof(double));
    
    if (!solver.solve_keep_jacobian(coeff_vec)) 
      error_function("Newton's iteration failed.");
    else
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, solutions);
    
    solver_time.tick();
    cpu_time.tick();
    
    // Show intermediate solutions.
    view1.show(&sln1);    
    view2.show(&sln2);
    view3.show(&sln3);    
    view4.show(&sln4);
    
    solver_time.tick(Hermes::HERMES_SKIP);
    cpu_time.tick(Hermes::HERMES_SKIP);
    
    // Compute eigenvalue.    
    SupportClasses::SourceFilter source(solutions, matprop, fission_regions, HERMES_AXISYM_Y);
    SupportClasses::SourceFilter source_prev(iterates, matprop, fission_regions, HERMES_AXISYM_Y);
    
    double k_new = k_eff * (source.integrate() / source_prev.integrate());
    info("Largest eigenvalue: %.8g, rel. difference from previous it.: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;

    // Update eigenvalue.
    k_eff = k_new;
    wf.update_keff(k_eff);
    
    if (!done)
    {
      // Save solutions for the next iteration.
      iter1.copy(&sln1);    
      iter2.copy(&sln2);
      iter3.copy(&sln3);    
      iter4.copy(&sln4);
      
      it++;
    }
  }
  while (!done);
  
  delete [] coeff_vec;
  
  // Time measurement.
  cpu_time.tick();
  solver_time.tick(Hermes::HERMES_SKIP);
  
  // Print timing information.
  verbose("Average solver time for one power iteration: %g s", solver_time.accumulated() / it);
  
  // Show solutions.
  view1.show(&sln1);
  view2.show(&sln2);
  view3.show(&sln3);    
  view4.show(&sln4);
  
  // Skip visualization time.
  cpu_time.tick(Hermes::HERMES_SKIP);

  // Print timing information.
  verbose("Total running time: %g s", cpu_time.accumulated());
    
  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
