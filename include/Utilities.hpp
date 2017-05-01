//
// Created by matthias on 28.04.17.
//

#ifndef FIELDCALIBRATION_UTILITIES_H
#define FIELDCALIBRATION_UTILITIES_H

#endif //FIELDCALIBRATION_UTILITIES_H

#include "Laser.hpp"

// Split the laser track set into tracks that reached the expected exit point (within a configurable region) and others.
// First entry of the return vector is tracks that reach the exit point, second is the ones that do not reach it.
std::vector<Laser> ReachedExitPoint(const Laser&, float);
