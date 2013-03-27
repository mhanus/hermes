#include "hermes2d.h"
using namespace Hermes::Hermes2D;

#include "weakforms_neutronics.h"
using namespace Neutronics; 


//////  Geometric parameters.  /////////////////////////////////////////////////////////////////

// File with initial mesh specification.
const std::string mesh_file = "domain.mesh";

//////  Physical parameters.  /////////////////////////////////////////////////////////////////

const double extinction = 1e4;
const double thermalization = 1e-8;
