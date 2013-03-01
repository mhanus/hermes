#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain.mesh";

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const Hermes::vector<std::string> materials = HermesMultiArray<std::string>
  ("M1S")("M1")("Void");

const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "M1S", 6.400
)(
  "M1", 0.000
)(
  "Void", 0.000
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "M1S", row(0.2)
)(
  "M1", row(0.2)
)(
  "Void", row(0.0)  
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "M1S", page(matrix(row(0.19)))
)(
  "M1", page(matrix(row(0.19)))
)(
  "Void", page(matrix(row(0.00)))
);
