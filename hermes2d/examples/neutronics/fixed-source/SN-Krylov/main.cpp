#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"
#include "weakforms_neutronics.h"


const bool HERMES_ANG_VISUALIZATION = false;
const bool VTK_ANG_VISUALIZATION = false;
const bool HERMES_SCAL_VISUALIZATION = true;
const bool VTK_SCAL_VISUALIZATION = true;


const bool MULTIMESH = false;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
const int N = 24;                    
const int M = N_GROUPS * N*(N+2)/2;

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
const double PICARD_TOL = 5e-3;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 150; 
// Value for constant initial condition.
const double INIT_COND_CONST = 1.0; 

MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;
// Only for matrix_solver_type == SOLVER_AZTECOO :
const char* iterative_method = "bicgstab";
const char* preconditioner = "jacobi";


#define SN
//#define DIFFUSION

int main(int argc, char* args[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.setParamValue(Hermes::exceptionsPrintCallstack, 1);
  Hermes::HermesCommonApi.setParamValue(Hermes::matrixSolverType, matrix_solver_type);
  Hermes::Hermes2D::Hermes2DApi.setParamValue(Hermes::Hermes2D::numThreads, 1);
  
  // Time measurement.
  TimeMeasurable cpu_time;
  MeshReaderH2D mloader;
  
#ifdef SN  
  cpu_time.tick();
  
  Hermes::vector<Mesh *> meshes;
  for (int i = 0; i < M; i++)
    meshes.push_back(new Mesh());
  
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
    }
  }
  for (int j = 0; j < INIT_REF_NUM; j++) 
    meshes[0]->refine_all_elements();


  MeshView mview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  mview.show(meshes[0]);
  Views::View::wait();
  
  Loggable::Static::info("%d elements", meshes[0]->get_num_active_elements());
  
  
  Hermes::vector<const Space<double> *> spaces;
  
  for (int n = 0; n < M; n++)
    spaces.push_back(new L2Space<double>(MULTIMESH ? meshes[n] : meshes[0], P_INIT));
  
  int ndof =  Space<double>::get_num_dofs(spaces);

  // Display the mesh.
//  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
//  oview.show(spaces[0]);
//  BaseView<double> bview("Shape functions", new WinGeom(450, 0, 440, 350));
//  bview.show(spaces[0]);

  Hermes::vector<Solution<double>* > slns;
  for (int i = 0; i < M; i++)
    slns.push_back(new Solution<double>(MULTIMESH ? meshes[i] : meshes[0]));
  
    // Load material data.
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, rm_map);
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_sn(Ssn);
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;

  // Initialize the weak formulation.
  Hermes::vector<std::string> reflective_boundaries;
  reflective_boundaries.push_back("reflective");
  SNWeakForm wf(N, matprop, slns, reflective_boundaries);
 
  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, spaces);
  dp.set_fvm();
  
  LinearSolver<double> solver(&dp);
  
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
  
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());
  
  Loggable::Static::info("num_assemble_DG = %d", dp.num_assemble_DG);
  
  // Translate the resulting coefficient vector into instances of Solution.
  Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, slns);
  
    // View the coarse mesh solution.
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
      Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(slns[n], (std::string("sln_") + itos(n) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
    
  
  Hermes::vector<MeshFunction<double>*> scalar_fluxes;
  SupportClasses::MomentFilter::get_scalar_fluxes(slns, &scalar_fluxes, N_GROUPS, wf.get_ordinates_data());
  
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    if (HERMES_SCAL_VISUALIZATION)
    {
      ScalarView view("Scalar flux", new WinGeom(0+450*g, 0, 450, 350));
      view.fix_scale_width(60);
      view.show(scalar_fluxes[g]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_SCAL_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(scalar_fluxes[g], (std::string("scalar_flux_g_") + itos(g) + std::string(".vtk")).c_str(), "Solution", mode_3D);
    }
  }
  
  SupportClasses::MomentFilter::clear_scalar_fluxes(&scalar_fluxes);

  
  const SupportClasses::OrdinatesData& odata = wf.get_ordinates_data();
  Node *ed;
  for_all_edge_nodes(ed, meshes[0])
  {
    if (ed->bnd)
    {
      // n->elem[0] != NULL
      Element *el = ed->elem[0];
      
      int ied = 0;
      for ( ; ied < el->get_nvert(); ied++)
        if (el->en[ied]->id == ed->id)
          break;
      
      double p1x = el->vn[ied]->x;
      double p1y = el->vn[ied]->y;
      double p2x = el->vn[(ied + 1) % el->get_nvert()]->x;
      double p2y = el->vn[(ied + 1) % el->get_nvert()]->y;
      double p3x = el->vn[(ied + 2) % el->get_nvert()]->x;
      double p3y = el->vn[(ied + 2) % el->get_nvert()]->y;

      double tx = p2x - p1x;
      double ty = p2y - p1y;
      double nx = -ty;
      double ny = tx;
      nx /= (Hermes::sqrt(nx*nx + ny*ny));
      ny /= (Hermes::sqrt(nx*nx + ny*ny));
      
      double cx = (p1x + p2x + p3x)/3.;
      double cy = (p1y + p2y + p3y)/3.;
      std::cout << "Elem: " << el->id << " T=[" << cx << "," << cy << "];  t=(" << tx << "," << ty << ")\n";
      
      double sum_pos = 0.;
      double sum_neg = 0.;
      
      int cp = 0;
      int cn = 0;
      for (int dir = 0; dir < odata.M; dir++)
      {  
        std::cout << "D" << dir << " : (" << odata.xi[dir] << "," << odata.eta[dir] << ")" << std::endl;
        double a_dot_n = odata.xi[dir]*nx + odata.eta[dir]*ny;
        std::cout << a_dot_n << std::endl;
        
        if (a_dot_n > 0)
        {
          cp++;
          sum_pos += /*a_dot_n * odata.pw[dir] * */slns[dir]->get_pt_value(cx,cy);
        }
        else if (a_dot_n < 0)
        {
          cn++;
          sum_neg += /*a_dot_n * odata.pw[dir] * */slns[dir]->get_pt_value(cx,cy);
        }
        else
          std::cout << "warning" << std::endl;
      }
      
      std::cout << std::endl << sum_pos << " ---- " << sum_neg << " (" << cp << "/" << cn << ")" << std::endl;
    }
  }
  
#endif
  
  

  
  
  // diffusion part
  
  
  
  
 
#ifdef DIFFUSION  
  cpu_time.tick();
  
  Hermes::vector<Mesh *> meshes_d;
  for (int i = 0; i < N_GROUPS; i++)
    meshes_d.push_back(new Mesh());
  
  mloader.load(mesh_file.c_str(), meshes_d[0]);
  
  if (MULTIMESH)
  {  
    for (int i = 1; i < N_GROUPS; i++) 
    {
      // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
      meshes_d[i]->copy(meshes_d[0]);
          
      // Initial uniform refinements.
      for (int j = 0; j < INIT_REF_NUM; j++) 
        meshes_d[i]->refine_all_elements();
    }
  }
  for (int j = 0; j < INIT_REF_NUM; j++) 
    meshes_d[0]->refine_all_elements();
  
  MaterialPropertyMap2 Ss;
  MaterialPropertyMap3::const_iterator it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss[it->first] = it->second[0];
  
  Diffusion::MaterialProperties::MaterialPropertyMaps matprop_d(N_GROUPS, rm_map);
  matprop_d.set_Sigma_t(St);
  matprop_d.set_Sigma_s(Ss);
  matprop_d.set_iso_src(src);
  matprop_d.validate();
  
  // Print material data.
  std::cout << matprop_d;
  
  Hermes::vector<const Space<double> *> spaces_d;
  
  for (int n = 0; n < N_GROUPS; n++)
    //spaces_d.push_back(new H1Space<double>(MULTIMESH ? meshes_d[n] : meshes_d[0], 1));
    spaces_d.push_back(new L2Space<double>(MULTIMESH ? meshes_d[n] : meshes_d[0], P_INIT));
  
  int ndof_d =  Space<double>::get_num_dofs(spaces_d);
  
  Hermes::vector<Solution<double>* > slns_d;
  for (int i = 0; i < N_GROUPS; i++)
    slns_d.push_back(new Solution<double>(MULTIMESH ? meshes_d[i] : meshes_d[0]));
  
  DiffusionWeakForm wf_d(meshes_d[0], true, matprop_d);
  DiscreteProblemLinear<double> dp_d(&wf_d, spaces_d);
  dp_d.set_fvm();
  LinearSolver<double> solver_d(&dp_d);
  solver_d.set_verbose_output(true);
  
  Loggable::Static::info("Solving. NDOF = %d", ndof_d);
  cpu_time.tick();

  try
  {
    solver_d.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
  
  Solution<double>::vector_to_solutions(solver_d.get_sln_vector(), spaces_d, slns_d);
  
  cpu_time.tick();
  Loggable::Static::info("Time taken: %lf s", cpu_time.last());
  
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    if (HERMES_SCAL_VISUALIZATION)
    {
      ScalarView view("Scalar flux", new WinGeom(0+450*g, 0, 450, 350));
      view.fix_scale_width(60);
      view.show(slns_d[g]);
      Views::View::wait();
    }
    
    // VTK output.
    if(VTK_SCAL_VISUALIZATION)
    {
      // Output solution in VTK format.
      Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(slns_d[g], (std::string("scalar_flux_g_") + itos(g) + std::string("_dif.vtk")).c_str(), "Solution", mode_3D);
    }
  }
#endif
  return 0;
}
