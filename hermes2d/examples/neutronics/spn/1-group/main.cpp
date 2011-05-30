#define HERMES_REPORT_ALL
#include "hermes2d.h"

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

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
