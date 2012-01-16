#include "neutronics/material_properties.h"
using namespace Hermes::Hermes2D::Neutronics; 

// Reference eigenvalue.
const double REF_K_INF = 1.0;

//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "clanek.mesh";

// Boundary markers.
const std::string bdy_vacuum = "2";
const std::string bdy_symmetry = "1";

const RegionMaterialMap rm_map = region_material_map
  ("1", "Palivo")
  ("2", "Palivo")
  ("3", "Palivo")
  ("4", "Palivo")
  ("5", "Pokryti")
  ("6", "Voda");

const Hermes::vector<std::string> fission_regions = Hermes::vector<std::string>(
  "1", "2", "3", "4"
);

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap1 D = material_property_map<rank1>
(
  "Palivo", row(1.92E-1)(1.84)
)(
  "Pokryti", row(1.23E-2)(4.09E-2)
)(
  "Voda", row(9.22E-1)(6.15E-2)
);

const MaterialPropertyMap1 Sa = material_property_map<rank1>
(
  "Palivo", row(3.26E-2)(6.88E-1)
)(
  "Pokryti", row(3.49E-4)(1.05E-2)
)(
  "Voda", row(4.12E-4)(1.68E-2)
);

const MaterialPropertyMap1 Sf = material_property_map<rank1>
(
  "Palivo", row(1.31E-2)(5.70E-1)
)(
  "Pokryti", row(0)(0)
)(
  "Voda", row(0)(0)
);

const MaterialPropertyMap1 nu = material_property_map<rank1>
(
  "Palivo", row(2.47)(2.44)
)(
  "Pokryti", row(0)(0)
)(
  "Voda", row(0)(0)
);

const MaterialPropertyMap2 Ss = material_property_map<rank2>
(
  "Palivo", 
  matrix(
    row(2.87E-1)(1.97E-4)
  )(
    row(4.26E-4)(2.16E-1)
  )
)(
  "Pokryti",
  matrix(
    row(1.72E-1)(7.08E-5)
  )(
    row(2.13E-4)(8.89E-2)
  )  
)(
  "Voda",
  matrix(
    row(8.27E-1)(3.74E-4)
  )(
    row(4.19E-2)(2.93E+0)
  )    
);
