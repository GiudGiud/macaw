//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VacuumParticleBC.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"

registerMooseObject("MaCawApp", VacuumParticleBC);

InputParameters
VacuumParticleBC::validParams()
{
  auto params = GeneralRayBC::validParams();
  params.addClassDescription("A RayBC that kills a particle on a boundary.");
  return params;
}

VacuumParticleBC::VacuumParticleBC(const InputParameters & params) : GeneralRayBC(params) {}

void
VacuumParticleBC::onBoundary(const unsigned int /* num_applying */)
{
  // After RayBCs are completed, ray->shouldContinue() is checked and this will kill the Ray
  currentRay()->setShouldContinue(false);

  // Contribute to leakage tally
#pragma omp atomic
  openmc::global_tally_leakage += currentRay()->auxData(1);
}
