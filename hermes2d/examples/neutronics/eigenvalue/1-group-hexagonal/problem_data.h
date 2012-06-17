#include "neutronics/material_properties.h"
using namespace Hermes::Hermes2D::Neutronics; 

// Reference k_effective reactor eigenvalue.
const double REF_K_EFF = 1.000332; // SP3
//const double REF_K_EFF = 1.001271; // SP5

//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "core.mesh";

// Boundary markers.
const std::string bdy_vacuum = "Vacuum";
const std::string bdy_symmetry = "Reflective";

const RegionMaterialMap rm_map = region_material_map
  ("Assembly 1", "Mixture 1")
  ("Assembly 2", "Mixture 1")
  ("Assembly 3", "Mixture 1")
  ("Assembly 4", "Mixture 1")
  ("Assembly 5", "Mixture 1")
  ("Assembly 6", "Mixture 1")
  ("Assembly 7", "Mixture 2")
  ("Assembly 8", "Mixture 1")
  ("Assembly 9", "Mixture 2")
  ("Assembly 10", "Mixture 2")
  ("Assembly 11", "Mixture 2")
  ("Assembly 12", "Mixture 2")
  ("Assembly 13", "Mixture 2")
  ("Assembly 14", "Mixture 2")
  ("Assembly 15", "Mixture 2")
  ("Assembly 16", "Mixture 3")
  ("Assembly 17", "Mixture 2")
  ("Assembly 18", "Mixture 2")
  ("Assembly 19", "Mixture 3")
  ("Assembly 20", "Mixture 3")
  ("Assembly 21", "Mixture 3")
  ("Assembly 22", "Mixture 3")
  ("Assembly 23", "Mixture 3")
  ("Assembly 24", "Mixture 3");

  
const Hermes::vector<std::string> fission_materials = HermesMultiArray<std::string>("Mixture 1");

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "Mixture 1", row(0.025)
)(
  "Mixture 2", row(0.025)
)(
  "Mixture 3", row(0.075)  
);

const MaterialPropertyMap1 nSf = material_property_map<rank1>
(
  "Mixture 1", row(0.0155)
)(
  "Mixture 2", row(0.0)
)(
  "Mixture 3", row(0.0)  
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "Mixture 1", page(matrix(row(0.013)))(matrix(row(0.0)))
)(
  "Mixture 2", page(matrix(row(0.024)))(matrix(row(0.006)))
)(
  "Mixture 3", page(matrix(row(0.0)))(matrix(row(0.0)))  
);

const rank0 nu = 2.43;