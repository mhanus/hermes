#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::MaterialProperties; 
using namespace Definitions; 

//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain.mesh";

const RegionMaterialMap rm_map = region_material_map
  ("pin 0", "Mixture 3")
  ("pin 1", "Mixture 2")
  ("pin 2", "Mixture 2")
  ("pin 3", "Mixture 2")
  ("pin 4", "Mixture 2")
  ("pin 5", "Mixture 2")
  ("pin 6", "Mixture 2")
  ("pin 7", "Mixture 2")
  ("pin 8", "Mixture 2")
  ("pin 9", "Mixture 2")
  ("cell 0", "Mixture 1")
  ("cell 1", "Mixture 1")
  ("cell 2", "Mixture 1")
  ("cell 3", "Mixture 1")
  ("cell 4", "Mixture 1")
  ("cell 5", "Mixture 1")
  ("cell 6", "Mixture 1")
  ("cell 7", "Mixture 1")
  ("cell 8", "Mixture 1")
  ("cell 9", "Mixture 1");

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const MaterialPropertyMap0 src = material_property_map<rank0>
(
  "Mixture 1", 1.000
)(
  "Mixture 2", 0.000
)(
  "Mixture 3", 0.000
);

const MaterialPropertyMap1 St = material_property_map<rank1>
(
  "Mixture 1", row(1.250)
)(
  "Mixture 2", row(0.625)
)(
  "Mixture 3", row(14.00)  
);

const MaterialPropertyMap3 Ssn = material_property_map<rank3>
(
  "Mixture 1", page(matrix(row(1.242)))
)(
  "Mixture 2", page(matrix(row(0.355)))
)(
  "Mixture 3", page(matrix(row(0.000)))
);

//////  Reference solutions.  /////////////////////////////////////////////////////////////////

const int n_pins = 10;
Hermes::vector<string> edit_regions = HermesMultiArray<string>
  ("pin 0")("cell 0")
  ("pin 1")("cell 1")
  ("pin 2")("cell 2")
  ("pin 3")("cell 3")
  ("pin 4")("cell 4")
  ("pin 5")("cell 5")
  ("pin 6")("cell 6")
  ("pin 7")("cell 7")
  ("pin 8")("cell 8")
  ("pin 9")("cell 9");
  
const double ref_integrated_fluxes[2*n_pins] = {
  2.3312E-01, 3.4114E+00, 
  1.0931E+01, 1.7956E+01, 
  1.2242E+01, 2.0084E+01, 
  1.2606E+01, 2.0660E+01, 
  1.1760E+01, 1.9264E+01, 
  2.4831E+01, 4.0680E+01, 
  2.5331E+01, 4.1491E+01, 
  1.2631E+01, 2.0688E+01, 
  2.5508E+01, 4.1771E+01, 
  1.2814E+01, 2.0982E+01
};
const double ref_integrated_absorption_rates[2*n_pins] = {
  3.2636E+00, 2.7291E-02,
  2.9513E+00, 1.4365E-01,
  3.3054E+00, 1.6067E-01,
  3.4037E+00, 1.6528E-01,
  3.1751E+00, 1.5411E-01,
  6.7043E+00, 3.2544E-01,
  6.8394E+00, 3.3193E-01,
  3.4103E+00, 1.6550E-01,
  6.8873E+00, 3.3417E-01,
  3.4598E+00, 1.6786E-01,
};
const double ref_regions_areas[2*n_pins] = {
  6.36173E-01, 9.26327E-01, 2.54469E+00, 3.70531E+00, 2.54469E+00,
  3.70531E+00, 2.54469E+00, 3.70531E+00, 2.54469E+00, 3.70531E+00,
  5.08938E+00, 7.41062E+00, 5.08938E+00, 7.41062E+00, 2.54469E+00,
  3.70531E+00, 5.08938E+00, 7.41062E+00, 2.54469E+00, 3.70531E+00
};