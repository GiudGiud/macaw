//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// local includes
#include "CollisionKernel.h"
#include "MaCawUtils.h"

// openmc includes
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/material.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/random_lcg.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"

// Extra includes for debugging
#include "openmc/bank.h"

registerMooseObject("MaCawApp", CollisionKernel);

InputParameters
CollisionKernel::validParams()
{
  InputParameters params = IntegralRayKernelBase::validParams();

  params.addClassDescription("Computes distance to the next collision. "
                             "Samples the reaction. Computes trajectory changes. "
                             "Contributes to tallies.");

  params.addRequiredCoupledVar("temperature", "The temperature of the medium.");
  params.addRequiredParam<std::vector<SubdomainName>>("blocks",
                                                      "Blocks on which this kernel is defined.");
  params.addRequiredParam<std::vector<unsigned int>>("materials",
                                                     "OpenMC material ids for each block.");
  params.addParam<Real>("z_coord", 0, "The axial coordinate for 2D calculations");
  params.addParam<bool>("verbose", false, "Whether to print collision information");

  return params;
}

CollisionKernel::CollisionKernel(const InputParameters & params)
  : IntegralRayKernelBase(params),
    _T(coupledValue("temperature")),
    _is_2D(_mesh.dimension() == 2),
    _z_coord(getParam<Real>("z_coord")),
    _verbose(getParam<bool>("verbose"))
{
  // Check that the temperature variable is a constant monomial
  // TODO check properly
  if (getFieldVar("temperature", 0) && getFieldVar("temperature", 0)->order() != 0)
    paramError("temperature",
               "Only CONST MONOMIAL temperatures (current order: ",
               getFieldVar("temperature", 0)->order(),
               ") are supported.");

  // Check for 2D parameters
  if (_is_2D && !params.isParamSetByUser("z_coord"))
    mooseError("The z_coord parameter must be specified for 2D calculations");
  if (!_is_2D && params.isParamSetByUser("z_coord"))
    mooseError("The z_coord parameter must not be specified for 3D calculations");

  // Build map from subdomains to OpenMC materials
  const auto blocks = _mesh.getSubdomainIDs(getParam<std::vector<SubdomainName>>("blocks"));
  const auto materials = getParam<std::vector<unsigned int>>("materials");
  if (blocks.size() != materials.size())
    paramError("blocks", "The blocks parameter must be the same size as the materials parameter.");
  // materials are indexed into openmc::model:::materials directly, not a map by their id
  std::vector<unsigned int> mat_indexes(materials.size());
  for (unsigned int i = 0; i < materials.size(); i++)
    for (unsigned int j = 0; j < openmc::model::materials.size(); j++)
      if ((int)materials[i] == openmc::model::materials[j]->id_)
        mat_indexes[i] = j;
  for (unsigned int i = 0; i < blocks.size(); i++)
    _block_to_openmc_materials.insert(std::make_pair<int, int>(blocks[i], mat_indexes[i]));

  // Resize the neutrons objects used to call OpenMC routines and initialize the seeds
  // TODO optimization: create these neutrons in the study, re-use them everywhere
  _particles.resize(1); // libMesh::n_threads());
  // WE CANNOT DO THIS because the global simulation weight is modified
  // for (unsigned int i = 0; i < _particles.size(); i++)
  //   openmc::initialize_history(_particles[i], i + 1);
}

void
CollisionKernel::initialSetup()
{
  if (_verbose)
    _console << "Kernel initial setup" << std::endl;

  // Check that all materials do exist in OpenMC, otherwise it will crash at XS calculation
  for (auto & m : getParam<std::vector<unsigned int>>("materials"))
  {
    auto search = openmc::model::material_map.find(m);
    if (search == openmc::model::material_map.end())
      mooseError("Could not find material ", m, " in OpenMC materials.");
  }

  // Check that all blocks exist, as a sanity check only
  for (auto & b : _block_to_openmc_materials)
    if (!hasBlocks(b.first))
      mooseWarning("Could not find specified block ", b.first, " for the collision kernel.");
}

void
CollisionKernel::onSegment()
{
  // Use a fake neutron to compute the cross sections
  auto p = &_particles[0]; //_tid]; (kernel already threaded)

  if (_verbose)
    std::cout << "In kernel onSegment" << std::endl;

  // Reset the particle since we are reusing them
  // Do not clear the particle, we rely on it for keeping last_sqrtkbT and last_mat
  // p->clear();

  // Resize to account for filters specified in moose (only needed once)
  // TODO Move to inital setup once can show that all tallies have already been added
  // THIS CANNOT BE DONE CURRENTLY.
  // Note that the ray tracing cost is way higher anyway
  p->filter_matches().resize(openmc::model::tally_filters.size());

  // Set particle attributes
  p->sqrtkT() = std::sqrt(openmc::K_BOLTZMANN * 293.6); // _T[0]);   ////////////////
  const auto mat = _block_to_openmc_materials.find(_current_subdomain_id);
  if (mat == _block_to_openmc_materials.end())
    mooseError("Material not found in subdomain " + std::to_string(_current_subdomain_id));
  p->material() = mat->second;
  p->coord(p->n_coord() - 1).universe = _current_subdomain_id;
  p->coord(p->n_coord() - 1).cell = _current_elem->id();
  if (_is_2D)
  {
    const Real factor = sqrt(currentRay()->direction().norm_sq() /
                             (1 - currentRay()->auxData(9) * currentRay()->auxData(9)));
    p->u() = {currentRay()->direction()(0) / factor,
              currentRay()->direction()(1) / factor,
              currentRay()->auxData(9)};
    mooseAssert(MooseUtils::absoluteFuzzyEqual(p->u().norm(), 1),
                "Should be unit direction " +
                    Moose::stringify(Point(p->u()[0], p->u()[1], p->u()[2])));
  }
  else
    p->u() = {
        currentRay()->direction()(0), currentRay()->direction()(1), currentRay()->direction()(2)};
  p->E() = currentRay()->auxData(0);
  p->wgt() = currentRay()->auxData(1);
  p->n_progeny() = currentRay()->auxData(2);

  // Set the particle id at the right value for setting n_progeny in event_death
  p->id() = currentRay()->auxData(3);

  // To avoid an overflow in the n_progeny array for neutrons which would have
  // changed domain, we adjust the value, however this prevents us from using the
  // openmc sorting algorithm when using domain decomposition
  if (comm().size() > 1)
  {
    if (p->id() - 1 - openmc::simulation::work_index[comm().rank()] >=
        openmc::simulation::work_per_rank)
      p->id() = 1 + openmc::simulation::work_index[comm().rank()];
    else if (p->id() - 1 - openmc::simulation::work_index[comm().rank()] < 0)
      p->id() = 1 + openmc::simulation::work_index[comm().rank()];
  }

  // Previous tracking information
  // p->sqrtkT_last() = currentRay()->auxData(13);
  // p->material_last() = currentRay()->auxData(14);

  // Set the particle seed for consistent random number generation
  p->seeds(0) =
      (unsigned long long)(currentRay()->auxData(4)) + (long long)currentRay()->auxData(10);
  p->seeds(1) =
      (unsigned long long)(currentRay()->auxData(5)) + (long long)currentRay()->auxData(11);
  p->seeds(2) =
      (unsigned long long)(currentRay()->auxData(6)) + (long long)currentRay()->auxData(12);

  if (_verbose)
  {
    std::cout << "seeds " << p->seeds()[0] << " " << p->seeds()[1] << " " << p->seeds()[2] << " "
              << *p->current_seed() << std::endl;
    std::cout << "Skipping evaluation ?" << p->sqrtkT() << " " << p->sqrtkT_last() << " "
              << p->material() << " " << p->material_last() << std::endl;
  }

  // Set particle type
  p->type() = openmc::ParticleType(currentRay()->auxData(7));

  if (p->type() == openmc::ParticleType::electron || p->type() == openmc::ParticleType::positron ||
      p->type() == openmc::ParticleType::photon)
    mooseError("Unsupported particle");

  // Reset the OpenMC particle status
  p->alive() = true;

  // Set particle number of events
  p->n_event() = currentRay()->auxData(8);

  // Compute all cross sections
  p->event_calculate_xs();

  // Compute distance to next collision
  // p.event_advance();
  // TODO scores track-length tallies as well
  // TODO Can we use this??

  // Adjust for 2D
  const Real effective_segment_length =
      _is_2D ? _current_segment_length / cos(p->u()[2]) : _current_segment_length;

  const auto rand = openmc::prn(p->current_seed());
  const Real collision_distance = -std::log(rand) / p->macro_xs().total;
  const Real distance = std::min(collision_distance, effective_segment_length);

  // Contribute to the track-length k-eff estimator
  if (openmc::settings::run_mode == openmc::RunMode::EIGENVALUE &&
      p->type() == openmc::ParticleType::neutron)
    p->keff_tally_tracklength() += p->wgt() * distance * p->macro_xs().nu_fission;

  // Keep track of the particle seed for consistent random number generation
  Real round, diff;
  MaCaw::convertInt64(p->seeds(0), round, diff);
  currentRay()->auxData(4) = round;
  currentRay()->auxData(10) = diff;
  MaCaw::convertInt64(p->seeds(1), round, diff);
  currentRay()->auxData(5) = round;
  currentRay()->auxData(11) = diff;
  MaCaw::convertInt64(p->seeds(2), round, diff);
  currentRay()->auxData(6) = round;
  currentRay()->auxData(12) = diff;

  // Longer travel than the distance to the next intersection, ray tracing will take care of moving
  // the particle
  if (collision_distance > effective_segment_length)
  {
    // Previous tracking information
    // this was previously done in openmc::geometry::find_cell_inner
    // currentRay()->auxData(13) = p->sqrtkT();
    // currentRay()->auxData(14) = p->material();

    // p.event_cross_surface();
    // TODO Score track length tallies
    // TODO Score surface tallies
    return;
  }
  else
  {
    // Advance ray (really, move backwards)
    // if (_verbose)
    //   std::cout << currentRay()->currentPoint() << effective_segment_length << " "
    //             << collision_distance << " " << Point(p->u()[0], p->u()[1], p->u()[2]) <<
    //             std::endl;
    Point current_position = currentRay()->currentPoint() -
                             _current_segment_length * currentRay()->direction() +
                             collision_distance * Point(p->u()[0], p->u()[1], p->u()[2]);

    // Adapt for two dimensions
    if (_is_2D)
      current_position(2) = _z_coord;

    // Set the particle position to the collision location for banking sites
    p->r() = {current_position(0), current_position(1), current_position(2)};

    // Compute collision
    p->event_collide();

    mooseAssert(!std::isnan(p->u()[0]) && !std::isnan(p->u()[1]) && !std::isnan(p->u()[2]) &&
                    !std::isnan(p->E() && p->E() >= 0),
                "Invalid particle energy and direction after collision.");

    if (_verbose)
      _console << "Collision event " << int(p->event()) << " Energy " << currentRay()->auxData(0)
               << " -> " << p->E() << " block " << _current_elem->subdomain_id() << " material "
               << p->material() << " progeny " << p->n_progeny() << " ids " << p->id() << " / "
               << currentRay()->id() << std::endl;

    // Scattering direction is stored on neutron
    Point new_direction(p->u()[0], p->u()[1], p->u()[2]);

    // Adapt for two dimensions
    if (_is_2D)
    {
      currentRay()->auxData(9) = new_direction(2);
      new_direction(2) = 0;
    }

    if (p->event() == openmc::TallyEvent::SCATTER)
      changeRayStartDirection(current_position, new_direction);

    // Update Ray energy
    currentRay()->auxData(0) = p->E();

    // Keep track of weight (for implicit capture)
    currentRay()->auxData(1) = p->wgt();

    // Keep track of particle number of progeny
    currentRay()->auxData(2) = p->n_progeny();

    // Keep track of the particle seed for consistent random number generation
    Real round, diff;
    MaCaw::convertInt64(p->seeds(0), round, diff);
    currentRay()->auxData(4) = round;
    currentRay()->auxData(10) = diff;
    MaCaw::convertInt64(p->seeds(1), round, diff);
    currentRay()->auxData(5) = round;
    currentRay()->auxData(11) = diff;
    MaCaw::convertInt64(p->seeds(2), round, diff);
    currentRay()->auxData(6) = round;
    currentRay()->auxData(12) = diff;

    // Keep track of particle number of events
    currentRay()->auxData(8) = p->n_event();

    // Mark the ray as 'should not continue' if absorption was sampled
    // TODO Handle secondary particles
    // Need to create a new ray as soon as the particles are sampled, as the kernel can only
    // change direction and create rays in the same element
    // TODO Need to track ray ids for progeny when adding new rays
    if (p->event() == openmc::TallyEvent::KILL || p->event() == openmc::TallyEvent::ABSORB)
    {
      currentRay()->setShouldContinue(false);
      const auto current_progeny = openmc::simulation::progeny_per_particle[p->id() - 1];
      p->event_death();

      // Particles are sharing ids, can't overwrite the progeny number
      // FIXME: This isnt quite right, it's assigning more progeny to some banked fission sites
      if (comm().size() > 1)
        openmc::simulation::progeny_per_particle[p->id() - 1] += current_progeny;
    }

    // Previous tracking information
    // this was previously done in openmc::geometry::find_cell_inner
    // this is not useful since we are changing energy
    // currentRay()->auxData(13) = p->sqrtkT();
    // currentRay()->auxData(14) = p->material();

    // Check for unsupported event
    if (p->event() != openmc::TallyEvent::KILL && p->event() != openmc::TallyEvent::ABSORB &&
        p->event() != openmc::TallyEvent::SCATTER)
      mooseError("Unsupported event", int(p->event()));

    // Handle secondary particles immediately after they are sampled, as the particle
    // could change domain, and rays cannot be created in a different domain
    // p->alive() = false;
    // p->event_revive_from_secondary();
    // NOTE: event_revive_from_secondary is currently disabled in openmc
    if (true)
      return;

    // A secondary particle has been sampled
    // TODO: allow all particles, not just neutrons
    // TODO: handle multiple secondary particles
    if (p->n_event() == 0 && p->type() == openmc::ParticleType::neutron)
    {
      // Create a new Ray with starting information
      Point start(p->r()[0], p->r()[1], p->r()[2]);
      Point direction(p->u()[0], p->u()[1], p->u()[2]);
      if (_is_2D)
      {
        start(2) = _z_coord;
        currentRay()->auxData(9) = direction(2);
        // direction(2) = 0;  // cant do that, must keep real direction
      }
      std::shared_ptr<Ray> new_ray = acquireRay(start, direction);

      if (true)
        _console << "Sampled secondary particle type " << int(p->type()) << ", creating new ray "
                 << new_ray->id() << " at " << p->r() << " " << p->E() << "eV " << std::endl;

      // Store neutron information
      new_ray->auxData(0) = p->E();
      new_ray->auxData(1) = p->wgt();

      // Reset number of progeny particles
      new_ray->auxData(2) = 0;

      // Need to generate a new particle id
      // FIXME This is not deterministic
      // Particle id is used in event_death to track the number
      // of progeny per particle, then used for sorting the fission source
      openmc::Particle p2;
      new_ray->auxData(3) = p2.id();

      // Keep track of the particle seed for consistent random number generation
      // TODO conversion of uint64
      new_ray->auxData(4) = p->seeds(0);
      new_ray->auxData(5) = p->seeds(1);
      new_ray->auxData(6) = p->seeds(2);

      // Keep track of particle type
      new_ray->auxData(7) = int(p->type());

      moveRayToBuffer(new_ray);
    }
  }
}
