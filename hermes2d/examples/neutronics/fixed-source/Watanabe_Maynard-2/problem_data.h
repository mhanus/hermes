#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "../domain.mesh";

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const Hermes::vector<std::string> materials = HermesMultiArray<std::string>
  ("I")("II")("III")("IV")("V")("Void");

const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "I", 0.000
)(
  "II", 0.000
)(
  "III", 0.000
)(
  "IV", 0.000
)(
  "V", 0.160
)(
  "Void", 0.000
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "I", row(1.0)
)(
  "II", row(2.0)
)(
  "III", row(1.0)
)(
  "IV", row(2.0)
)(
  "V", row(0.1)
)(
  "Void", row(0.0)  
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "I", page(matrix(row(0.9)))
)(
  "II", page(matrix(row(1.5)))
)(
  "III", page(matrix(row(0.9)))
)(
  "IV", page(matrix(row(1.2)))
)(
  "V", page(matrix(row(0.1)))
)(
  "Void", page(matrix(row(0.0)))
);

Hermes::vector<string> edit_regions = HermesMultiArray<string>
  ("III")("IV");
