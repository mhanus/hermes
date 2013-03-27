#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain.mesh";

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const double scattering_ratio = 1-1e-8;
const double sigma_t = 1e4;
