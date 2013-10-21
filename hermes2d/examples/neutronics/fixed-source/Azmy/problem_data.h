#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "../domain.mesh";

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const Hermes::vector<std::string> materials = HermesMultiArray<std::string>
  ("M1")("M2");

const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "M1", 1.000
)(
  "M2", 0.000
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "M1", row(1.0)
)(
  "M2", row(2.0)
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "M1", page(matrix(row(0.5)))
)(
  "M2", page(matrix(row(0.1)))
);
