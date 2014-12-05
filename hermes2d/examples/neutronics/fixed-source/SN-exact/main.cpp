#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"
#include "weakforms_neutronics.h"

#include <iterator>
#include <iomanip>
const bool SAVE_SLN_VECTOR = true;

const bool HERMES_ANG_VISUALIZATION = false;
const bool VTK_ANG_VISUALIZATION = true;
const bool HERMES_SCAL_VISUALIZATION = true;
const bool VTK_SCAL_VISUALIZATION = true;
const bool HERMES_MESH_VISUALIZATION = false;

const bool MULTIMESH = false;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 4;
const int REF_TO_BND = 0;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 1;

const int N = 14;          
const int M = N*(N+2)/2;

//
// Parameters for the source iteration via the Picard's method.
//
// Number of last iterations used.
// 1... standard fixed point.
// >1... Anderson acceleration.
const int PICARD_NUM_LAST_ITER_USED = 3;
// 0 <= beta <= 1... parameter for the Anderson acceleration.
const double PICARD_ANDERSON_BETA = 0;
// Stopping criterion for the Picard's method.
const double PICARD_TOL = 1e-5;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 1000; 
// Value for constant initial condition.
const double INIT_COND_CONST = 1.0; 

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;

int main(int argc, char* args[])
{
  // Set the number of threads used in Hermes.
  
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::matrixSolverType, matrix_solver_type);
  //Hermes::HermesCommonApi.set_integral_param_value(Hermes::numThreads, 2);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  cpu_time.tick();
  
  Hermes::vector<MeshSharedPtr > meshes;
  for (int i = 0; i < M; i++)
    meshes.push_back(MeshSharedPtr(new Mesh()));
  
  MeshReaderH2D mloader;
  mloader.load(mesh_file.c_str(), meshes[0]);
  
  if (MULTIMESH)
  {  
    for (int i = 1; i < M; i++) 
    {
      // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
      meshes[i]->copy(meshes[0]);
          
      // Initial uniform refinements.
      for (int j = 0; j < INIT_REF_NUM; j++) 
        meshes[i]->refine_all_elements();
      
      meshes[i]->refine_towards_boundary("vacuum", REF_TO_BND);
    }
  }
  for (int j = 0; j < INIT_REF_NUM; j++) 
    meshes[0]->refine_all_elements();
  meshes[0]->refine_towards_boundary("vacuum", REF_TO_BND);
  
  if (HERMES_MESH_VISUALIZATION)
  {
    MeshView mview("Coarse mesh", new WinGeom(0, 0, 440, 350));
    mview.show(meshes[0]);
    Views::View::wait();
  }
  
  Loggable::Static::info("%d elements, %d vertices", meshes[0]->get_num_active_elements(), meshes[0]->get_num_vertex_nodes() );
  
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for (int i = 0; i < M; i++)
    slns.push_back(new Solution<double>(MULTIMESH ? meshes[i] : meshes[0]));
 
  // Initialize the weak formulation.
  Hermes::vector<std::string> reflective_boundaries;
  Hermes::vector<std::string> vacuum_boundaries(1);
  vacuum_boundaries.push_back("vacuum");
  
  switch (argc)
  {
    case 2:
    {
      if (!strcmp(args[1], "const"))
      {
        L2Space<double> sp1(meshes[0], 2);
        int ndof = sp1.get_num_dofs();
        
        Loggable::Static::info("NDOF = %d", ndof);
        
        double *ones_v = new double [ndof];
        for (int i = 0; i < 1; i++)
          ones_v[i] = 1.0;
        Solution<double> ones(meshes[0]);
        Solution<double>::vector_to_solution(ones_v, &sp1, &ones);
        
        if (HERMES_SCAL_VISUALIZATION)
        {
          BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
          bview.show(&sp1);

          ScalarView view("Const", new WinGeom(0, 0, 450, 350));
          view.fix_scale_width(60);
          view.show(&ones);
          Views::View::wait();
        }
        
        delete [] ones_v;
        return 0;
      }
      
      bool assemble_matrix = strcmp(args[1], "Q");
      bool assemble_Q = !strcmp(args[1], "Q") || !strcmp(args[1], "LQ");
      
      WeakForm<double> *wf;
      Hermes::vector<SpaceSharedPtr<double> > spaces;
      
      if (!strcmp(args[1], "S") || !strcmp(args[1], "F"))
      {
        wf = new IsotropicScatteringMatrixForm(scattering_ratio, sigma_t);
        spaces.push_back(new L2Space<double>(meshes[0], P_INIT));
      }
      else
      {
        wf = new SNWeakForm(N, scattering_ratio, sigma_t, reflective_boundaries, vacuum_boundaries, args[1]);  
        for (int n = 0; n < M; n++)
          spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
      }
      
      Loggable::Static::info("Saving %s. NDOF = %d", args[1], Space<double>::get_num_dofs(spaces));
      
      DiscreteProblem<double> dp(wf, spaces);
      SourceIteration solver(&dp);
      
      // Perform the source iteration (by Picard's method with Anderson acceleration).
      solver.set_max_allowed_iterations(1);
      
      if (assemble_Q)  
      {
        solver.output_rhs(1);
        solver.set_rhs_export_format(EXPORT_FORMAT_MATLAB_MATIO);
        solver.set_rhs_filename("Q");
        solver.set_rhs_number_format("%1.15f");
        solver.set_rhs_varname("Q");
      }
      if (assemble_matrix)
      {
        solver.output_matrix(1);
        solver.set_matrix_export_format(EXPORT_FORMAT_MATLAB_MATIO);
        solver.set_matrix_filename(!strcmp(args[1], "LQ") ? "L" : args[1]);
        solver.set_matrix_number_format("%1.15f");
        solver.set_matrix_varname(!strcmp(args[1], "LQ") ? "L" : args[1]);
      }
      
      try 
      { 
        solver.solve(); 
      } 
      catch(std::exception& e) { }
      
      delete wf;
      
      return 0;
    }
    case 3:
    {
      if (!strcmp(args[1], "-sln"))
      {
        Hermes::vector<SpaceSharedPtr<double> > spaces;
  
        for (int n = 0; n < M; n++)
          spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
        
        int ndof =  Space<double>::get_num_dofs(spaces);
    
        std::vector<double> x_ext;
        x_ext.reserve(ndof);
        std::ifstream ifs( args[2] , std::ifstream::in );
        read_solution_from_file(ifs, std::back_inserter(x_ext));
        ifs.close();
        
        Hermes::vector<MeshFunctionSharedPtr<double> > sol_ext;
        for (unsigned int i = 0; i < M; i++) 
          sol_ext.push_back(new Solution<double>()); 
        Solution<double>::vector_to_solutions(x_ext.data(), spaces, sol_ext);
        
        Hermes::vector<MeshFunctionSharedPtr<double> > scalar_fluxes;
        SupportClasses::OrdinatesData odata(N, "lgvalues.txt");
        SupportClasses::MomentFilter::get_scalar_fluxes(sol_ext, &scalar_fluxes, 1, odata);
        
        if (HERMES_SCAL_VISUALIZATION)
        {
          ScalarView view("Scalar flux", new WinGeom(0, 0, 450, 350));
          view.fix_scale_width(60);
          view.show(scalar_fluxes[0]);
          Views::View::wait();
        }
        
        // VTK output.
        if(VTK_SCAL_VISUALIZATION)
        {
          // Output solution in VTK format.
          Linearizer lin(FileExport);
          bool mode_3D = false;
          lin.save_solution_vtk(scalar_fluxes[0], "scalar_flux.vtk", "Solution", mode_3D);
        }
      }
      
      return 0;
    }
  }
  
  SNWeakForm wf(N, scattering_ratio, sigma_t, reflective_boundaries, vacuum_boundaries);
 
  // Initialize the FE problem.
  
  // Approximation spaces.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  
  for (int n = 0; n < M; n++)
    spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
  
  int ndof =  Space<double>::get_num_dofs(spaces);
  
  // Display the mesh.
//  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
//  oview.show(spaces[0]);
//  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
//  bview.show(spaces[0]);

  // Discrete formulation.
  DiscreteProblem<double> dp(&wf, spaces);
  
  // Algebraic solver.
  SourceIteration solver(&dp);
  
  solver.use_Anderson_acceleration(false);
  solver.set_tolerance(PICARD_TOL, Solvers::NonlinearConvergenceMeasurementType::SolutionChangeRelative);
  solver.set_max_allowed_iterations(PICARD_MAX_ITER);
  solver.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
  solver.set_anderson_beta(PICARD_ANDERSON_BETA);
  solver.set_verbose_output(true);
  
  Loggable::Static::info("Solving. NDOF = %d", ndof);
  cpu_time.tick();
  
  try
  {
    solver.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
  
  double *sln_vector = solver.get_sln_vector();
  Solution<double>::vector_to_solutions(sln_vector, spaces, slns);
  
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());

  // View the coarse mesh solution and/or save it to .vtk files.
  for (int n = 0; n < M; n++)
  {
    if(HERMES_ANG_VISUALIZATION)
    {
      ScalarView view("Solution", new WinGeom(0+450*n, 400, 450, 350));
      view.fix_scale_width(60);
      view.show(slns[n]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_ANG_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin(FileExport);
      bool mode_3D = false;
      lin.save_solution_vtk(slns[n], (std::string("sln_") + tostr(n) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
  // View the coarse mesh scalar flux and/or save it to .vtk files.
  
  Hermes::vector<MeshFunctionSharedPtr<double> > scalar_fluxes;
  SupportClasses::MomentFilter::get_scalar_fluxes(slns, &scalar_fluxes, 1, wf.get_ordinates_data());
  
  if (HERMES_SCAL_VISUALIZATION)
  {
    ScalarView view("Scalar flux", new WinGeom(0, 0, 450, 350));
    view.fix_scale_width(60);
    view.show(scalar_fluxes[0]);
    Views::View::wait();
  }
  
  // VTK output.
  if(VTK_SCAL_VISUALIZATION)
  {
    // Output solution in VTK format.
    Linearizer lin(FileExport);
    bool mode_3D = false;
    lin.save_solution_vtk(scalar_fluxes[0], "scalar_flux.vtk", "Solution", mode_3D);
  }

  if (SAVE_SLN_VECTOR)
  {
    std::string file = "x-R"+tostr(INIT_REF_NUM)+"P"+tostr(P_INIT)+"-S"+tostr(N)+".dat";

    std::ofstream fs(file.c_str());
    Loggable::Static::info("Saving the solution vector to %s", file.c_str());

      fs << std::setprecision(16);
    std::copy(sln_vector, sln_vector+ndof, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
  }
  
  return 0;
}
