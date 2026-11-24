// include/physicsforge/physicsforge.h
//
// Umbrella header: includes the complete physicsforge C++ API.
//
//     #include <physicsforge/physicsforge.h>
//
// Prefer including the individual domain headers in large translation units to
// keep compile times down.

#ifndef PHYSICSFORGE_PHYSICSFORGE_H
#define PHYSICSFORGE_PHYSICSFORGE_H

#include "physicsforge/config.h"
#include "physicsforge/constants.h"
#include "physicsforge/version.h"

#include "physicsforge/math/linalg.h"
#include "physicsforge/math/vector3.h"

#include "physicsforge/mechanics/dynamics.h"
#include "physicsforge/mechanics/gravitation.h"
#include "physicsforge/mechanics/kinematics.h"
#include "physicsforge/mechanics/rocketry.h"
#include "physicsforge/mechanics/rotation.h"

#include "physicsforge/electromagnetism/circuits.h"
#include "physicsforge/electromagnetism/electrostatics.h"
#include "physicsforge/electromagnetism/maxwell.h"

#include "physicsforge/thermodynamics/heat.h"
#include "physicsforge/thermodynamics/heat_equation.h"

#include "physicsforge/modern/quantum.h"
#include "physicsforge/modern/relativity.h"

#include "physicsforge/numerics/calculus.h"
#include "physicsforge/numerics/integrate.h"
#include "physicsforge/numerics/optimize.h"
#include "physicsforge/numerics/roots.h"
#include "physicsforge/numerics/transforms.h"

#endif // PHYSICSFORGE_PHYSICSFORGE_H
