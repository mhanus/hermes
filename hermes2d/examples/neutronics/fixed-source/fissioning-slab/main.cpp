//
#define HERMES_REPORT_ALL
#include "definitions.h"
#include <iterator>

using namespace RefinementSelectors;

const unsigned int N_GROUPS = 1;    // Monoenergetic (single group) problem.
#ifdef USE_SPN
  const unsigned int SPN_ORDER = 9; // Currently implemented maximum is 9.
#else // DIFFUSION
  const unsigned int SPN_ORDER = 1;
#endif

const unsigned int N_MOMENTS = SPN_ORDER+1;
const unsigned int N_ODD_MOMENTS = N_MOMENTS/2;
const unsigned int N_EQUATIONS = N_GROUPS * N_ODD_MOMENTS;

// Initial uniform mesh refinement and initial polynomial orders for the individual solution components. 

#ifdef USE_SPN
  const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    5,//  5,     // SP1
    5,//  5,     // SP3
    5,//  5,     // SP5
    5,//  5,     // SP7
    5/*,  5/*      // SP9
  */};
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    2,//  2,     // SP1
    2,//  2,     // SP3
    2,//  2,     // SP5
    2,//  2,     // SP7
    2/*,  2/*      // SP9
  */};
#else // DIFFUSION
  const int INIT_REF_NUM[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    5//,  5
  };
  const int P_INIT[N_EQUATIONS] = {
  /* g1 g2 */
  /*-------*/
    2//,  2
  };
#endif

const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric. Absolute element errors should be
                                         //   used in Adapt::calc_err_est for optimal multi-mesh adaptivity.
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
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent). Setting a very small value
                                         // effectively disables this convergence test in favor of the following ones.
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit.
const int MAX_ADAPT_NUM = 60;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
                                         
// Solvers used for solving the discrete problems. Possibilities (depends on how you configured the hermes library): 
//  SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                         
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

// Visualisation options.
const bool HERMES_VISUALIZATION = true;         // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;           // Set to "true" to enable VTK output.
const bool DISPLAY_MESHES = false;              // Set to "true" to display initial mesh data. Requires HERMES_VISUALIZATION == true.
const bool SHOW_INTERMEDIATE_ORDERS = true;     // Set to "true" to display coarse mesh solutions during adaptivity.
const bool SHOW_INTERMEDIATE_SOLUTIONS = true;  // Set to "true" to display solutions on intermediate meshes during adaptivity.

const bool SAVE_FLUX_PROFILES = true;

// Problem parameters:

// Isotropic distributed sources:
/// MG
/*
const MaterialPropertyMap1 src = material_property_map<rank1>
(
  "mixture", row(1.0)(0.01)
);
*/
/// SG
const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "mixture", 1.0
);


// Total cross-sections:
/// MG
/*
const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "mixture", row(0.65)(1.30)
);
*/
/// SG
const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "mixture", row(0.65)
);


// Scattering cross-sections:
/// MG ANISO
/*
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "mixture", 
  page(
    matrix(
      row(0.15)(0.0)
    )(
      row(0.4)(0.18)
    )
  )(
    matrix(
      row(0.08)(0.0)
    )(
      row(0.09)(0.03)
    )
  )
);
*/
/// MG ISO
/*
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "mixture", 
  page(
    matrix(
      row(0.15)(0.0)
    )(
      row(0.4)(0.18)
    )
  )
);
*/
/// SG ANISO
/*
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "mixture", 
  page(
    matrix(
      row(0.15)
    )
  )(
    matrix(
      row(0.08)
    )
  )
);
*/
/// SG ISO
const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "mixture", 
  page(
    matrix(
      row(0.6)
    )
  )
);
// Reference solutions - average fluxes by composition:
const Hermes::vector<std::string> material_region = HermesMultiArray<std::string>
  ("mixture");
 
/// MG ANISO
/*
const double average_flux_dragon[N_GROUPS] = {
  1.3432932E+00, 4.2068257E-01
};
*/
/// MG ISO
/*
const double average_flux_dragon[N_GROUPS] = {
  1.3694813E+00, 4.3017557E-01
};
*/
/// SG ANISO
/*
const double average_flux_dragon[N_GROUPS] = {
  1.3432932E+00
};
*/
/// SG ISO
/*
const double average_flux_dragon[N_GROUPS] = {
  1.3694818E+00
};
*/
/// SG ISO Ss=0.60
const double average_flux_dragon[N_GROUPS] = {
  3.8764552E+00
};

const int n_diag_pts = 50;
const double slab_width = 2.5;
const double diag_step = slab_width / n_diag_pts;

/// MG ANISO
/*
const double diag_flux_dragon[N_GROUPS][n_diag_pts] = 
{
  { 
    1.6471964E+00, 1.6468259E+00, 1.6457715E+00, 1.6442480E+00, 1.6423539E+00, 1.6399459E+00,
    1.6368761E+00, 1.6334189E+00, 1.6295545E+00, 1.6249522E+00, 1.6198722E+00, 1.6143723E+00,
    1.6082629E+00, 1.6014713E+00, 1.5942188E+00, 1.5864602E+00, 1.5778728E+00, 1.5687313E+00,
    1.5590579E+00, 1.5486161E+00, 1.5374028E+00, 1.5255820E+00, 1.5130710E+00, 1.4995700E+00,
    1.4853394E+00, 1.4703689E+00, 1.4543649E+00, 1.4373800E+00, 1.4195155E+00, 1.4006340E+00,
    1.3804492E+00, 1.3591832E+00, 1.3367779E+00, 1.3128543E+00, 1.2874879E+00, 1.2606966E+00, 
    1.2322352E+00, 1.2017661E+00, 1.1694320E+00, 1.1350447E+00, 1.0980088E+00, 1.0582978E+00, 
    1.0157118E+00, 9.6959854E-01, 9.1918631E-01, 8.6413646E-01, 8.0340727E-01, 7.3442234E-01, 
    6.5438861E-01, 5.5598901E-01
  },
  { 
    5.7022607E-01, 5.6998647E-01, 5.6949741E-01, 5.6876373E-01, 5.6778827E-01, 5.6656355E-01, 
    5.6508311E-01, 5.6335257E-01, 5.6136750E-01, 5.5911463E-01, 5.5659618E-01, 5.5381003E-01, 
    5.5074404E-01, 5.4738997E-01, 5.4374789E-01, 5.3980842E-01, 5.3555489E-01, 5.3098533E-01, 
    5.2609159E-01, 5.2085638E-01, 5.1526812E-01, 5.0932010E-01, 5.0299647E-01, 4.9627616E-01, 
    4.8915056E-01, 4.8160423E-01, 4.7361274E-01, 4.6515875E-01, 4.5622661E-01, 4.4679190E-01, 
    4.3682639E-01, 4.2631175E-01, 4.1522229E-01, 4.0352329E-01, 3.9118786E-01, 3.7818760E-01, 
    3.6448537E-01, 3.5004016E-01, 3.3481911E-01, 3.1877979E-01, 3.0186950E-01, 2.8404229E-01, 
    2.6524572E-01, 2.4541293E-01, 2.2446457E-01, 2.0231624E-01, 1.7884221E-01, 1.5382997E-01, 
    1.2687769E-01, 9.6660183E-02 
  }
};
*/
/// MG ISO
/*
const double diag_flux_dragon[N_GROUPS][n_diag_pts] = 
{
  { 
    1.7068595E+00, 1.7063930E+00, 1.7052898E+00, 1.7036646E+00, 1.7015694E+00, 1.6989241E+00, 1.6956456E+00,
    1.6918782E+00, 1.6876071E+00, 1.6826498E+00, 1.6771390E+00, 1.6711002E+00, 1.6644249E+00, 1.6570652E+00,
    1.6491302E+00, 1.6405865E+00, 1.6312530E+00, 1.6212610E+00, 1.6106101E+00, 1.5991635E+00, 1.5869018E+00,
    1.5738941E+00, 1.5600801E+00, 1.5452853E+00, 1.5296195E+00, 1.5130549E+00, 1.4954181E+00, 1.4767037E+00,
    1.4569324E+00, 1.4360032E+00, 1.4137351E+00, 1.3901936E+00, 1.3653062E+00, 1.3388384E+00, 1.3107640E+00,
    1.2810284E+00, 1.2494447E+00, 1.2157629E+00, 1.1799510E+00, 1.1418105E+00, 1.1009309E+00, 1.0571260E+00,
    1.0101065E+00, 9.5933127E-01, 9.0409966E-01, 8.4380099E-01, 7.7733592E-01, 7.0250556E-01, 6.1597624E-01,
    5.0912603E-01
  },
  { 
    5.9032728E-01, 5.9006756E-01, 5.8954154E-01, 5.8875193E-01, 5.8769975E-01, 5.8637944E-01, 5.8478561E-01,
    5.8292045E-01, 5.8077957E-01, 5.7835257E-01, 5.7563877E-01, 5.7263435E-01, 5.6932906E-01, 5.6571471E-01,
    5.6178802E-01, 5.5753977E-01, 5.5295553E-01, 5.4803002E-01, 5.4275363E-01, 5.3711097E-01, 5.3108965E-01,
    5.2468001E-01, 5.1786660E-01, 5.1063020E-01, 5.0295894E-01, 4.9483645E-01, 4.8624063E-01, 4.7715306E-01,
    4.6755611E-01, 4.5742674E-01, 4.4673896E-01, 4.3547224E-01, 4.2360168E-01, 4.1109665E-01, 3.9793057E-01,
    3.8407584E-01, 3.6950001E-01, 3.5416748E-01, 3.3804744E-01, 3.2110388E-01, 3.0329507E-01, 2.8458275E-01,
    2.6492536E-01, 2.4427455E-01, 2.2257308E-01, 1.9976048E-01, 1.7575056E-01, 1.5039539E-01, 1.2338820E-01,
    9.3649836E-02
  }
};
*/
/// SG ANISO
/*
const double diag_flux_dragon[N_GROUPS][n_diag_pts] = 
{
  { 
    1.6471964E+00, 1.6468259E+00, 1.6457715E+00, 1.6442480E+00, 1.6423539E+00, 1.6399459E+00,
    1.6368761E+00, 1.6334189E+00, 1.6295545E+00, 1.6249522E+00, 1.6198722E+00, 1.6143723E+00,
    1.6082629E+00, 1.6014713E+00, 1.5942188E+00, 1.5864602E+00, 1.5778728E+00, 1.5687313E+00,
    1.5590579E+00, 1.5486161E+00, 1.5374028E+00, 1.5255820E+00, 1.5130710E+00, 1.4995700E+00,
    1.4853394E+00, 1.4703689E+00, 1.4543649E+00, 1.4373800E+00, 1.4195155E+00, 1.4006340E+00,
    1.3804492E+00, 1.3591832E+00, 1.3367779E+00, 1.3128543E+00, 1.2874879E+00, 1.2606966E+00, 
    1.2322352E+00, 1.2017661E+00, 1.1694320E+00, 1.1350447E+00, 1.0980088E+00, 1.0582978E+00, 
    1.0157118E+00, 9.6959854E-01, 9.1918631E-01, 8.6413646E-01, 8.0340727E-01, 7.3442234E-01, 
    6.5438861E-01, 5.5598901E-01
  }
};
*/
/// SG ISO
/*
const double diag_flux_dragon[N_GROUPS][n_diag_pts] = 
{
  { 
    1.7068534E+00, 1.7063869E+00, 1.7052839E+00, 1.7036589E+00, 1.7015640E+00, 1.6989190E+00, 1.6956408E+00,
    1.6918739E+00, 1.6876032E+00, 1.6826464E+00, 1.6771362E+00, 1.6710979E+00, 1.6644231E+00, 1.6570641E+00,
    1.6491297E+00, 1.6405866E+00, 1.6312537E+00, 1.6212623E+00, 1.6106120E+00, 1.5991659E+00, 1.5869048E+00,
    1.5738975E+00, 1.5600838E+00, 1.5452894E+00, 1.5296239E+00, 1.5130595E+00, 1.4954228E+00, 1.4767084E+00,
    1.4569370E+00, 1.4360077E+00, 1.4137394E+00, 1.3901975E+00, 1.3653098E+00, 1.3388414E+00, 1.3107666E+00,
    1.2810303E+00, 1.2494460E+00, 1.2157634E+00, 1.1799508E+00, 1.1418096E+00, 1.1009292E+00, 1.0571237E+00,
    1.0101037E+00, 9.5932800E-01, 9.0409610E-01, 8.4379733E-01, 7.7733236E-01, 7.0250234E-01, 6.1597358E-01,
    5.0912413E-01
  }
};
*/
/// SG ISO Ss=0.60
const double diag_flux_dragon[N_GROUPS][n_diag_pts] = 
{
  { 
    5.5233908E+00, 5.5198283E+00, 5.5124396E+00, 5.5014090E+00, 5.4868288E+00, 5.4685844E+00, 5.4465589E+00,
    5.4209975E+00, 5.3918980E+00, 5.3589988E+00, 5.3225314E+00, 5.2825692E+00, 5.2389748E+00, 5.1917078E+00,
    5.1409786E+00, 5.0867815E+00, 5.0288778E+00, 4.9675208E+00, 4.9027661E+00, 4.8344595E+00, 4.7626263E+00,
    4.6874455E+00, 4.6088931E+00, 4.5267700E+00, 4.4413238E+00, 4.3526011E+00, 4.2604219E+00, 4.1648649E+00,
    4.0660729E+00, 3.9639995E+00, 3.8584732E+00, 3.7497177E+00, 3.6377608E+00, 3.5223776E+00, 3.4036713E+00,
    3.2817318E+00, 3.1564557E+00, 3.0276293E+00, 2.8954212E+00, 2.7597658E+00, 2.6202818E+00, 2.4769527E+00,
    2.3296771E+00, 2.1780135E+00, 2.0212671E+00, 1.8590452E+00, 1.6903831E+00, 1.5127108E+00, 1.3220509E+00,
    1.1086001E+00
  }
};
int main(int argc, char* argv[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);

  // Time measurement.
  TimeMeasurable cpu_time;
  double total_cpu_time;
  cpu_time.tick();

  // Load material data (last argument specifies a list of material, which is required
  // for automatic fill-in of missing standard data (fission properties in this case)).
  
#ifndef USE_SPN
  // Extract the isotropic part from the scattering matrix.
  
  MaterialPropertyMap2 Ss0;
  MaterialPropertyMap3::const_iterator it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss0[it->first] = it->second[0];
#endif

#ifdef USE_SPN  

  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, SPN_ORDER, 
                                                   std::set<std::string>(material_region.begin(), material_region.end()));
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);

#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)  

  MaterialPropertyMap2 Ss1;
  it = Ssn.begin();
  for ( ; it != Ssn.end(); it++)
    Ss1[it->first] = it->second[1];
  
  MaterialProperties::TransportCorrectedMaterialPropertyMaps matprop(N_GROUPS, Ss1);
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_s(Ss0);
  
#else // SIMPLE DIFFUSION
  
  MaterialProperties::MaterialPropertyMaps matprop(N_GROUPS, 
                                                   std::set<std::string>(material_region.begin(), material_region.end()));
  matprop.set_Sigma_t(St);
  matprop.set_Sigma_s(Ss0);

#endif
    
  matprop.set_iso_src(src);
  matprop.validate();
  
  // Print material data.
  std::cout << matprop;
  
  // Use multimesh, i.e. create one mesh for each energy group and pseudo-flux.
  Hermes::vector<MeshSharedPtr > meshes;
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
    meshes.push_back(new Mesh());
  
  // Load the mesh on which the 1st solution component (1st group, 0th moment) will be approximated.
  MeshReaderH2D mesh_reader;
  mesh_reader.load("domain.mesh", meshes[0]);
  
  for (unsigned int i = 1; i < N_EQUATIONS; i++) 
  {
    // Obtain meshes for the subsequent components by cloning the mesh loaded for the 1st one.
    meshes[i]->copy(meshes[0]);
        
    // Initial uniform refinements.
    for (int j = 0; j < INIT_REF_NUM[i]; j++) 
      meshes[i]->refine_all_elements();
    //meshes[i]->refine_towards_boundary("vacuum", 4);
  }
  for (int j = 0; j < INIT_REF_NUM[0]; j++) 
    meshes[0]->refine_all_elements();
  //meshes[0]->refine_towards_boundary("vacuum", 4);

#ifdef USE_SPN
  SupportClasses::Visualization views(SPN_ORDER, N_GROUPS, DISPLAY_MESHES);
#else // DIFFUSION
  SupportClasses::Visualization views(N_GROUPS, DISPLAY_MESHES);
#endif

  if (DISPLAY_MESHES && HERMES_VISUALIZATION)
    views.inspect_meshes(meshes);
  
  // Create pointers to the coarse and fine mesh solutions.
  Hermes::vector<MeshFunctionSharedPtr<double> > coarse_solutions, solutions;
  
  // Initialize all the new solution variables.
  for (unsigned int i = 0; i < N_EQUATIONS; i++)
  {
    coarse_solutions.push_back(new Solution<double>());
    solutions.push_back(new Solution<double>());
  }
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  for (unsigned int i = 0; i < N_EQUATIONS; i++) 
    spaces.push_back(new H1Space<double>(meshes[i], P_INIT[i]));
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, SPN_ORDER);
  
  // Get material_region areas for flux averaging.
#ifdef USE_SPN
  PostProcessor pp(NEUTRONICS_SPN);      
#else // DIFFUSION
  PostProcessor pp(NEUTRONICS_DIFFUSION);
#endif

  double area = pp.get_area(meshes[0], material_region);
  Loggable::Static::info("material_region \"%s\": A=%gcm2", material_region[0].c_str(), area);
  
  if (STRATEGY >= 0)
  { 
    // DOF and CPU convergence graphs initialization.
    GnuplotGraph graph_dof("Convergence of rel. errors", "NDOF", "error [%]");
  #ifdef USE_SPN
    graph_dof.add_row("pseudo-fluxes (H1)", "r", "-", "+");
  #else // DIFFUSION
    graph_dof.add_row("fluxes (H1)", "r", "-", "+");
  #endif
    graph_dof.add_row("region-average total flux w.r.t. DRAGON", "b", "-", "+");
    graph_dof.add_row("region-average total flux est.", "b", "--", "o");
    graph_dof.add_row("diagonal total flux w.r.t. DRAGON", "g", "-", "+");
    graph_dof.add_row("diagonal total flux est.", "g", "--", "o");
    
    graph_dof.set_log_x();
    graph_dof.set_log_y();
    graph_dof.show_legend();
    graph_dof.show_grid();
    
    GnuplotGraph graph_cpu("Convergence of rel. errors", "CPU time [s]", "error [%]");
  #ifdef USE_SPN  
    graph_cpu.add_row("pseudo-fluxes (H1)", "r", "-", "+");
  #else // DIFFUSION
    graph_cpu.add_row("fluxes (H1)", "r", "-", "+");
  #endif
    graph_cpu.add_row("region-average total flux w.r.t. DRAGON", "b", "-", "+");
    graph_cpu.add_row("region-average total flux est.", "b", "--", "o");
    graph_cpu.add_row("diagonal total flux w.r.t. DRAGON", "g", "-", "+");
    graph_cpu.add_row("diagonal total flux est.", "g", "--", "o");
        
    graph_cpu.set_log_x();
    graph_cpu.set_log_y();
    graph_cpu.show_legend();
    graph_cpu.show_grid();
    
    // Initialize the refinement selectors.
    H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    Hermes::vector<RefinementSelectors::Selector<double>*> selectors;
    for (unsigned int i = 0; i < N_EQUATIONS; i++)
      selectors.push_back(&selector);
    
    // Adaptivity loop.
    int as = 1; bool done = false;
    do 
    {
      Loggable::Static::info("---- Adaptivity step %d:", as);
      
      // Initialize the fine mesh problem.
      Hermes::vector<SpaceSharedPtr<double> > fine_spaces = Space<double>::construct_refined_spaces(spaces);
      int ndof_fine = Space<double>::get_num_dofs(fine_spaces);
      
      report_num_dof("Solving on fine meshes, #DOF: ", fine_spaces);
    
      DiscreteProblem<double> dp(&wf, fine_spaces);

      // Perform Newton's iteration on reference mesh.
      NewtonSolver<double> *newton = new NewtonSolver<double>(&dp);
      newton->set_verbose_output(false);
    
      try
      {
        newton->solve();
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.print_msg();
        ErrorHandling::error_function("Newton's iteration failed.");
      }
      
      // Translate the resulting coefficient vector into instances of Solution.
      Solution<double>::vector_to_solutions(newton->get_sln_vector(), fine_spaces, solutions);

      // Solution of discrete problem completed, delete the solver object.
      delete newton;
      
      // Project the fine mesh solution onto the coarse mesh.
      report_num_dof("Projecting fine-mesh solutions onto coarse meshes, #DOF: ", spaces);
      OGProjection<double> ogProjection;
      ogProjection.project_global(spaces, solutions, coarse_solutions);

      // View the coarse-mesh solutions and polynomial orders.
      if (HERMES_VISUALIZATION)
      {
        cpu_time.tick();
        Loggable::Static::info("Visualizing.");
        if (SHOW_INTERMEDIATE_SOLUTIONS)
          views.show_solutions(coarse_solutions);
        if (SHOW_INTERMEDIATE_ORDERS)
          views.show_orders(spaces);
        cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      }   

      cpu_time.tick();
      
      // Calculate element errors.
      Loggable::Static::info("Calculating error estimate.");
      
      Loggable::Static::info("  --- Calculating relative errors of scalar flux approximation.");
      
  #ifdef USE_SPN    
      Hermes::vector< MeshFunctionSharedPtr<double> >* coarse_scalar_fluxes = new Hermes::vector< MeshFunctionSharedPtr<double> >();
      Hermes::vector< MeshFunctionSharedPtr<double> >* fine_scalar_fluxes = new Hermes::vector< MeshFunctionSharedPtr<double> >();
      
      MomentFilter::get_scalar_fluxes(coarse_solutions, coarse_scalar_fluxes, N_GROUPS);
      MomentFilter::get_scalar_fluxes(solutions, fine_scalar_fluxes, N_GROUPS);
  #else // DIFFUSION
      Hermes::vector< MeshFunctionSharedPtr<double> >* coarse_scalar_fluxes = &coarse_solutions;
      Hermes::vector< MeshFunctionSharedPtr<double> >* fine_scalar_fluxes = &solutions;
  #endif

      double avg_flux_err_est_rel = 0.0;
      double avg_flux_err_dragon_rel = 0.0;
      double diag_flux_err_est_rel = 0.0;
      double diag_flux_err_dragon_rel = 0.0;
        
      for (unsigned int g = 0; g < N_GROUPS; g++)
      {      
        // Calculate relative error (squared) of the material_region-averaged fluxes (in l2 norm).
        double avg_group_coarse_flux = pp.integrate(coarse_scalar_fluxes->at(g), material_region) / area;
        double avg_group_fine_flux = pp.integrate(fine_scalar_fluxes->at(g), material_region) / area;
        
        // Estimate.
        avg_flux_err_est_rel += Hermes::sqr(avg_group_coarse_flux - avg_group_fine_flux) / Hermes::sqr(avg_group_fine_flux);
        
        // "Exact".
        avg_flux_err_dragon_rel += Hermes::sqr(avg_group_fine_flux - average_flux_dragon[g]) / Hermes::sqr(average_flux_dragon[g]);
        
        double x = diag_step/2.;
        for (int i = 0; i < n_diag_pts; i++, x+=diag_step)
        {
          // Estimate.
          {
            double approx = coarse_scalar_fluxes->at(g)->get_pt_value(x,x);
            double ref = fine_scalar_fluxes->at(g)->get_pt_value(x,x);
            diag_flux_err_est_rel += Hermes::sqr((approx - ref) / ref);
          }
        
          // "Exact".
          {
            double approx = fine_scalar_fluxes->at(g)->get_pt_value(x,x);
            double ref = diag_flux_dragon[g][i];
            diag_flux_err_dragon_rel += Hermes::sqr((approx - ref) / ref);
          }
        }
      }  
      
      avg_flux_err_est_rel = sqrt(avg_flux_err_est_rel) * 100;
      avg_flux_err_dragon_rel = sqrt(avg_flux_err_dragon_rel) * 100;
      diag_flux_err_est_rel = sqrt(diag_flux_err_est_rel) * 100;
      diag_flux_err_dragon_rel = sqrt(diag_flux_err_dragon_rel) * 100;        

  #ifdef USE_SPN    
      delete coarse_scalar_fluxes;
      delete fine_scalar_fluxes;
  #endif

      cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      
      // Calculate error estimate for each solution component and the total error estimate.
      Loggable::Static::info("  --- Calculating total relative error of the solution approximation.");
          
      Adapt<double> adaptivity(spaces);  
      
  #ifdef USE_SPN    
      MomentGroupFlattener mg(N_GROUPS);  // Creates a single index from the moment-group pair.
      // Set the error estimation/normalization form.
      for (unsigned int g = 0; g < N_GROUPS; g++)
        for (unsigned int mrow = 0; mrow < N_ODD_MOMENTS; mrow++)
          for (unsigned int mcol = 0; mcol < N_ODD_MOMENTS; mcol++)
            adaptivity.set_error_form(mg.pos(mrow,g), mg.pos(mcol,g), new ErrorFormSPN<double>(mrow, mcol, HERMES_H1_NORM));
  #endif            
      
      // Calculate the element and total error estimates and make them available for mesh adaptation.
      double solution_err_est_rel = adaptivity.calc_err_est(coarse_solutions, solutions, NULL, true,
                                                            HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100; 
              
      // Report results.
      Loggable::Static::info("  --- solutions rel. error estimate: %g%%", solution_err_est_rel);
      Loggable::Static::info("  --- l2 norm of estimated rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_est_rel);    
      Loggable::Static::info("  --- l2 norm of \"exact\" rel. errors of region-averaged scalar fluxes: %g%%", avg_flux_err_dragon_rel);
      Loggable::Static::info("  --- l2 norm of estimated rel. errors of diagonal scalar fluxes: %g%%", diag_flux_err_est_rel);    
      Loggable::Static::info("  --- l2 norm of \"exact\" rel. errors of diagonal scalar fluxes: %g%%", diag_flux_err_dragon_rel);
      
      // Add the results to convergence graphs.
      cpu_time.tick();
      graph_dof.add_values(0, ndof_fine, solution_err_est_rel);
      graph_cpu.add_values(0, cpu_time.accumulated(), solution_err_est_rel);
      graph_dof.add_values(1, ndof_fine, avg_flux_err_dragon_rel);
      graph_cpu.add_values(1, cpu_time.accumulated(), avg_flux_err_dragon_rel);
      graph_dof.add_values(2, ndof_fine, avg_flux_err_est_rel);
      graph_cpu.add_values(2, cpu_time.accumulated(), avg_flux_err_est_rel);
      graph_dof.add_values(3, ndof_fine, diag_flux_err_est_rel);
      graph_cpu.add_values(3, cpu_time.accumulated(), diag_flux_err_est_rel);
      cpu_time.tick(TimeMeasurable::HERMES_SKIP);
      
      // If err_est is too large, adapt the mesh.
      if (solution_err_est_rel < ERR_STOP || as == MAX_ADAPT_NUM || ndof_fine >= NDOF_STOP) 
        done = true;
      else 
      {
        Loggable::Static::info("Adapting the coarse meshes.");
        done = adaptivity.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);
      }
      
      if (!done) // Increase the adaptivity step counter.
        as++;
      else
      {
        cpu_time.tick();
        total_cpu_time = cpu_time.accumulated();
        Loggable::Static::info("Total running time: %g s", total_cpu_time);
        cpu_time.reset();
        
        if (HERMES_VISUALIZATION)
        {
          Loggable::Static::info("Visualizing final solutions.");
          views.inspect_solutions(solutions);
          views.inspect_orders(fine_spaces);
        }
      }
    }
    while (done == false);
    
    // Save the convergence graphs.
  #ifdef USE_SPN  
    graph_dof.save(("conv_dof_sp"+itos(SPN_ORDER)+".gp").c_str());
    graph_cpu.save(("conv_cpu_sp"+itos(SPN_ORDER)+".gp").c_str());
  #elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
    graph_dof.save("conv_dof_diffusion_trc.gp");
    graph_cpu.save("conv_cpu_diffusion_trc.gp");
  #else // SIMPLE DIFFUSION
    graph_dof.save("conv_dof_diffusion.gp");
    graph_cpu.save("conv_cpu_diffusion.gp");
  #endif
  }
  else
  {
    report_num_dof("Solving - #DOF: ", spaces);
    
    cpu_time.tick();
       
    DiscreteProblem<double> dp(&wf, spaces);
  
    // Perform Newton's iteration on reference mesh.
    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(true);
  
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      ErrorHandling::error_function("Newton's iteration failed.");
    }
  
    // Translate the resulting coefficient vector into instances of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, solutions);
    
    cpu_time.tick();
    total_cpu_time = cpu_time.accumulated();
    Loggable::Static::info("Total running time: %g s", total_cpu_time);
    cpu_time.reset();
        
    if (HERMES_VISUALIZATION)
    {
      Loggable::Static::info("Visualizing solutions.");
      views.inspect_solutions(solutions);
      views.inspect_orders(spaces);
    }
  }
  
  //
  // Analysis of the solution.
  //
  
  // Visualization.

#ifdef USE_SPN
  if (HERMES_VISUALIZATION)
  {
    views.show_all_flux_moments(solutions, matprop);    
    // Wait for the view to be closed.  
    Views::View::wait();
  }
#endif

  if (VTK_VISUALIZATION)
  {
    views.save_solutions_vtk("flux", "flux", solutions);
    views.save_orders_vtk("space", spaces);
  }
  
  // Output file names.
#ifdef USE_SPN  
  std::string file = "integrated_flux-sp"+itos(SPN_ORDER)+".dat";
#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
  std::string file = "integrated_flux-diffusion_trc.dat";
#else // SIMPLE DIFFUSION
  std::string file = "integrated_flux-diffusion.dat";
#endif

  std::cout << std::endl;
  Loggable::Static::info("Calculating errors of integrated fluxes w.r.t. DRAGON and saving to %s", file.c_str());

  FILE *fp = fopen(file.c_str(), "wt");
  fprintf(fp, "Total running time: %g s\n\n", total_cpu_time);
  
  for (unsigned int g = 0; g < N_GROUPS; g++)
  {
    double average_flux = pp.get_integrated_group_scalar_fluxes(solutions, g, N_GROUPS, material_region) / area;
        
    Loggable::Static::info("--- Average scalar flux in group %d : %f (error w.r.t DRAGON = %g%%)", g, average_flux,
         fabs(average_flux - average_flux_dragon[g])/average_flux_dragon[g] * 100);
    fprintf(fp, "--- Average scalar flux in group %d : %f (error w.r.t DRAGON = %g%%)\n", g, average_flux,
         fabs(average_flux - average_flux_dragon[g])/average_flux_dragon[g] * 100);
  }
  
  fclose(fp);
  std::cout << std::endl;
  
  if (SAVE_FLUX_PROFILES)
  {
#ifdef USE_SPN
    Hermes::vector< Filter<double>* >* scalar_fluxes =  new Hermes::vector< Filter<double>* >();
    MomentFilter::get_scalar_fluxes(solutions, scalar_fluxes, N_GROUPS);
#else // DIFFUSION
    Hermes::vector< MeshFunctionSharedPtr<double> >* scalar_fluxes;
    scalar_fluxes = &solutions;
#endif

    // Output file names.
#ifdef USE_SPN  
    file = "diagonal_flux-sp"+itos(SPN_ORDER)+".dat";
    std::string file_err = "diagonal_flux-sp"+itos(SPN_ORDER)+".err";
#elif defined(USE_DIFFUSION_WITH_TRANSPORT_CORRECTION)
    file = "diagonal_flux-diffusion_trc.dat";
    std::string file_err = "diagonal_flux-sp"+itos(SPN_ORDER)+".err";
#else // SIMPLE DIFFUSION
    file = "diagonal_flux-diffusion.dat";
    std::string file_err = "diagonal_flux-sp"+itos(SPN_ORDER)+".err";
#endif
   
    double *res = new double [n_diag_pts*N_GROUPS];
    double *err = new double [n_diag_pts*N_GROUPS];
    
    std::ofstream fs(file.c_str());
    std::ofstream fs_err(file_err.c_str());
    Loggable::Static::info("Saving the diagonal scalar flux profile to %s, error w.r.t DRAGON to %s", file.c_str(), file_err.c_str());
    
    for (unsigned int g = 0; g < N_GROUPS; g++)
    {
      std::cout << std::endl << "GROUP " << g << std::endl;
      
      double x = diag_step/2.;
      
      for (int i = 0; i < n_diag_pts; i++, x+=diag_step)
      {
        res[i+g*n_diag_pts] = scalar_fluxes->at(g)->get_pt_value(x, x);
        err[i+g*n_diag_pts] = fabs(res[i+g*n_diag_pts] - diag_flux_dragon[g][i]) / diag_flux_dragon[g][i] * 100;
        Loggable::Static::info("Flux at (%.3f, %.3f) : %g, error w.r.t. DRAGON : %g%%", x, x, res[i+g*n_diag_pts], err[i+g*n_diag_pts]);
      }

      std::copy(res+g*n_diag_pts, res+(g+1)*n_diag_pts, std::ostream_iterator<double>(fs, "\n"));
      std::copy(err+g*n_diag_pts, err+(g+1)*n_diag_pts, std::ostream_iterator<double>(fs_err, "\n"));
      fs << std::endl;
      fs_err << std::endl;
    }
    
    fs.close();
    fs_err.close();
    delete [] res;
    delete [] err;

#ifdef USE_SPN
    delete scalar_fluxes;
#endif    
  }

  return 0;
}
