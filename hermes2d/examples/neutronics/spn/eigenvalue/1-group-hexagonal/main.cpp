#define HERMES_REPORT_ALL
#include "definitions.h"
#include "problem_data.h"

// TODO

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("core.mesh", &mesh);

  // Conversion between triangular and quadrilateral meshes (optional). 
  //mesh.convert_quads_to_triangles();
  //mesh.convert_triangles_to_quads();

  // Refine mesh uniformly (optional).
  mesh.refine_all_elements();          

  // Display the mesh.
  MeshView mview("Core mesh", new WinGeom(0, 0, 350, 350));
  mview.show(&mesh);
  
  MaterialPropertyMaps matprop(G, N, rm_map);
  matprop.set_nuSigma_f(nSf);
  matprop.set_nu(nu);
  matprop.set_Sigma_tn(St);
  matprop.set_Sigma_sn(Ssn);
  
  matprop.validate();
  
  cout << matprop;
  
  // Wait for the view to be closed.
  View::wait();
  return 0;
}
