#include "neutronics/material_properties.h"
using namespace Hermes::Hermes2D::Neutronics; 

// Reference k_effective reactor eigenvalue.
const double REF_K_EFF = 1.00970; 

//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "core.mesh";

// Boundary markers.
const std::string bdy_vacuum = "Vacuum";
const std::string bdy_symmetry = "Reflective";

const RegionMaterialMap rm_map = region_material_map
  ("Assembly 1",  "Mixture 4")
  ("Assembly 2",  "Mixture 1")
  ("Assembly 3",  "Mixture 2")
  ("Assembly 4",  "Mixture 2")
  ("Assembly 5",  "Mixture 1")
  ("Assembly 6",  "Mixture 1")
  ("Assembly 7",  "Mixture 4")
  ("Assembly 8",  "Mixture 2")
  ("Assembly 9",  "Mixture 3")
  ("Assembly 10", "Mixture 2")
  ("Assembly 11", "Mixture 5")
  ("Assembly 12", "Mixture 2")
  ("Assembly 13", "Mixture 1")
  ("Assembly 14", "Mixture 2")
  ("Assembly 15", "Mixture 2")
  ("Assembly 16", "Mixture 2")
  ("Assembly 17", "Mixture 1")
  ("Assembly 18", "Mixture 1")
  ("Assembly 19", "Mixture 3")
  ("Assembly 20", "Mixture 3")
  ("Assembly 21", "Mixture 5")
  ("Assembly 22", "Mixture 2")
  ("Assembly 23", "Mixture 1")
  ("Assembly 24", "Mixture 1")
  ("Assembly 25", "Mixture 2")
  ("Assembly 26", "Mixture 2")
  ("Assembly 27", "Mixture 3")
  ("Assembly 28", "Mixture 3")
  ("Assembly 29", "Mixture 5")
  ("Assembly 30", "Mixture 2")
  ("Assembly 31", "Mixture 1")
  ("Assembly 32", "Mixture 1")
  ("Assembly 33", "Mixture 1")
  ("Assembly 34", "Mixture 3")
  ("Assembly 35", "Mixture 5")
  ("Assembly 36", "Mixture 2")
  ("Assembly 37", "Mixture 2")
  ("Assembly 38", "Mixture 3")
  ("Assembly 39", "Mixture 3")
  ("Assembly 40", "Mixture 5")
  ("Assembly 41", "Mixture 2")
  ("Assembly 42", "Mixture 3")
  ("Assembly 43", "Mixture 5")
  ("Assembly 44", "Mixture 5");

  
const Hermes::vector<std::string> fission_materials = HermesMultiArray<std::string>("Mixture 1")("Mixture 2")("Mixture 3");

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap1 D = material_property_map<rank1>
(
  "Mixture 1", row(1.3466)(0.37169)
)(
  "Mixture 2", row(1.3377)(0.36918)
)(
  "Mixture 3", row(1.3322)(0.36502)
)(
  "Mixture 4", row(1.1953)(0.19313)
)(
  "Mixture 5", row(1.4485)(0.25176)
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "Mixture 1", row(2.5255e-2)(6.4277e-2)
)(
  "Mixture 2", row(2.4709e-2)(7.9361e-2)
)(
  "Mixture 3", row(2.4350e-2)(1.0010e-1)
)(
  "Mixture 4", row(3.5636e-2)(1.3498e-1)
)(
  "Mixture 5", row(3.3184e-2)(3.2839e-2)
);

const MaterialPropertyMap1 Sf = material_property_map<rank1>
(
  "Mixture 1", row(2.21676e-3)(3.94368e-2)
)(
  "Mixture 2", row(2.79212e-3)(5.6572e-2)
)(
  "Mixture 3", row(3.59068e-3)(8e-2)
)(
  "Mixture 4", row(0)(0)
)(
  "Mixture 5", row(0)(0)
);

const MaterialPropertyMap1 nSf = material_property_map<rank1>
(
  "Mixture 1", row(4.4488e-3)(7.3753e-2)
)(
  "Mixture 2", row(5.5337e-3)(1.0581e-1)
)(
  "Mixture 3", row(7.0391e-3)(1.4964e-1)
)(
  "Mixture 4", row(0)(0)
)(
  "Mixture 5", row(0)(0)
);

const MaterialPropertyMap2 Ss = material_property_map<rank2>
(
  "Mixture 1", 
  matrix(
    row(0)(0)
  )(
    row(1.6893e-2)(0)
  )
)(
  "Mixture 2",
  matrix(
    row(0)(0)
  )(
    row(1.5912e-2)(0)
  )
)(
  "Mixture 3",
  matrix(
    row(0)(0)
  )(
    row(1.4888e-2)(0)
  )
)(
  "Mixture 4",
  matrix(
    row(0)(0)
  )(
    row(2.2264e-2)(0)
  )
)(
  "Mixture 5",
  matrix(
    row(0)(0)
  )(
    row(3.2262e-2)(0)
  )
);