#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain2.mesh";

const Hermes::vector<std::string> reflective_boundaries = HermesMultiArray<std::string>
  ("boxBL-bottom")("bottom")("boxBR-bottom")("boxBR-right")("right")("boxTR-right");
const Hermes::vector<std::string> vacuum_boundaries = HermesMultiArray<std::string>
  ("boxTR-top")("top")("boxTL-top")("boxTL-left")("left")("boxBL-left");

const RegionMaterialMap rm_map = region_material_map
  ("source", "source")
  ("boxBL", "strong_pure_scatterer")
  ("boxBR", "strong_pure_scatterer")
  ("boxTR", "detector")
  ("boxTL", "detector")
  ("obstacleTL", "strong_pure_absorber")
  ("obstacleTR", "weak_pure_absorber")
  ("medium", "vacuum");

Hermes::vector<string> edit_regions = HermesMultiArray<string>
  ("boxTR")("boxTL");
  
//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "source", 1.000
)(
  "strong_pure_scatterer", 0.000
)(
  "detector", 0.000
)(
  "strong_pure_absorber", 0.000
)(
  "weak_pure_absorber", 0.000
)(
  "vacuum", 0.000
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "source", row(0.1)
)(
  "strong_pure_scatterer", row(2.0)
)(
  "detector", row(1.0)
)(
  "strong_pure_absorber", row(2.0)
)(
  "weak_pure_absorber", row(0.5)
)(
  "vacuum", row(0.0)
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "source", page(matrix(row(0.05)))
)(
  "strong_pure_scatterer", page(matrix(row(2.0)))
)(
  "detector", page(matrix(row(0.2)))
)(
  "strong_pure_absorber", page(matrix(row(0.0)))
)(
  "weak_pure_absorber", page(matrix(row(0.0)))
)(
  "vacuum", page(matrix(row(0.0)))
);
