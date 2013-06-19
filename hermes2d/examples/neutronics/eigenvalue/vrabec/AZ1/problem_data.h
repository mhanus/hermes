#include "neutronics/material_properties.h"
using namespace Hermes::Hermes2D::Neutronics; 

// Reference eigenvalue.
const double REF_K_EFF = 1.0;

//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain.mesh";

// Boundary markers.
const std::string bdy_vacuum = "Vokraj";

const Hermes::vector<std::string> reflective_boundaries = Hermes::vector<std::string>();
const Hermes::vector<std::string> vacuum_boundaries = HermesMultiArray<std::string>
  (bdy_vacuum);

const Hermes::vector<std::string> fission_materials = HermesMultiArray<std::string>
  ("palivo_501")("palivo_511");

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

struct ProblemData
{
  static const int G = 2;
  
  MaterialPropertyMap1 D, Sa, Sf, nu;
  MaterialPropertyMap2 Ss;
  
  std::set<std::string> mat_list;
  
  ProblemData(const std::string& datafile);
};
