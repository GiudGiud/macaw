from math import log10

import matplotlib.pyplot as plt
import numpy as np
import openmc
import os

# Strong scaling study
# Increase the number of processors, keep the problem (neutrons, domain size) the same

mode = 'openmp'
mode = 'mpi'
scaling = 'strong'

num_procs = [1, 2, 4, 8, 16, 28, 32, 40, 50, 56, 64, 80, 112]
timings = [[], []]

# Some parameters for the study
n_batches = 100
n_active_batches = 90
n_neutrons = int(1e2)

grid_x = 10
grid_y = 10
grid_z = 10

for num_proc in num_procs:

    pin_pitch = 1
    if scaling == 'strong':
        box_size = 5
        neutrons = n_neutrons
    if scaling == 'weak':
        box_size = 1.25984 * np.sqrt(num_proc)
        neutrons = num_proc * n_neutrons

    ###############################################################################
    # Create materials for the problem

    fuel = openmc.Material(name='fuel')
    fuel.set_density('g/cc', 10.)
    fuel.add_element('U', 1., enrichment=5)
    fuel.add_element('O', 2.)

    water = openmc.Material(name='moderator')
    water.set_density('g/cc', 0.74)
    water.add_element('H', 5.0e-2)
    water.add_element('O', 2.4e-2)

    # Collect the materials together and export to XML
    materials = openmc.Materials([fuel, water])
    materials.export_to_xml()

    ###############################################################################
    # Define problem geometry

    # Create a region represented as the inside of a rectangular prism
    pitch = box_size
    pin = openmc.rectangular_prism(pin_pitch, pin_pitch)

    # get a few planes going
    x0 = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    x1 = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
    y0 = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    y1 = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')
    z0 = openmc.ZPlane(z0=-pitch/2, boundary_type='reflective')
    z1 = openmc.ZPlane(z0=+pitch/2, boundary_type='reflective')
    x0b = openmc.XPlane(x0=-pin_pitch/2)
    x1b = openmc.XPlane(x0=+pin_pitch/2)
    y0b = openmc.YPlane(y0=-pin_pitch/2)
    y1b = openmc.YPlane(y0=+pin_pitch/2)

    # Create cells, mapping materials to regions
    fuel = openmc.Cell(fill=fuel, region=+x0b&-x1b&+y0b&-y1b&+z0&-z1)
    water = openmc.Cell(fill=water, region=(+x0&-x1&+y0&-y1&+z0&-z1)&(-x0b|+x1b|-y0b|+y1b))

    # Create a geometry and export to XML
    geometry = openmc.Geometry([fuel, water])
    geometry.export_to_xml()

    ###############################################################################
    # Define problem settings

    # Indicate how many particles to run
    settings = openmc.Settings()
    settings.batches = n_batches
    settings.inactive = n_batches - n_active_batches
    settings.particles = neutrons

    # Create an initial uniform spatial source distribution over fissionable zones
    lower_left = (-pin_pitch/2, -pin_pitch/2, -pitch/2)
    upper_right = (pin_pitch/2, pin_pitch/2, pitch/2)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    settings.source = openmc.source.Source(space=uniform_dist)

    settings.export_to_xml()

    ###############################################################################
    # Define tallies

    lower_left = (-pitch/2, -pitch/2, -pitch/2)
    upper_right = (pitch/2, pitch/2, pitch/2)

    # Create a mesh that will be used for tallying
    mesh = openmc.RegularMesh()
    mesh.dimension = (grid_x, grid_y, grid_z)
    mesh.lower_left = (-pitch/2, -pitch/2, -pitch/2)
    mesh.upper_right = (pitch/2, pitch/2, pitch/2)

    # Create a mesh filter that can be used in a tally
    mesh_filter = openmc.MeshFilter(mesh)

    # Now use the mesh filter in a tally and indicate what scores are desired
    mesh_tally = openmc.Tally(name="Mesh tally")
    mesh_tally.filters = [mesh_filter]
    mesh_tally.scores = ['flux', 'fission', 'nu-fission']

    # Let's also create a tally to get the flux energy spectrum. We start by
    # creating an energy filter
    e_min, e_max = 1e-5, 20.0e6
    groups = 500
    energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
    energy_filter = openmc.EnergyFilter(energies)

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([mesh_tally])
    tallies.export_to_xml()

    ###############################################################################
    # Plot file
    plot1 = openmc.Plot()
    plot1.basis = 'xy'
    plot1.origin = (0, 0, 0)
    plot1.width = (pitch, pitch)
    plot1.pixels = (400, 400)

    plot2 = openmc.Plot()
    plot2.basis = 'xz'
    plot2.origin = (0, 0, 0)
    plot2.width = (pitch, pitch)
    plot2.pixels = (400, 400)

    plots = openmc.Plots([plot1, plot2])
    plots.export_to_xml()

    ###############################################################################
    # Run the code
    if mode == 'openmp':
        openmc.run("--threads", num_proc)
    elif mode == 'mpi':
        os.system("mpirun -n "+str(num_proc)+" openmc -s 1")
    else:
      print("Unrecognized run mode")

    # Get run time
    with openmc.StatePoint("statepoint.100.h5") as statepoint:
        timings[0].append(statepoint.runtime['active batches'])


# Plot results
plt.figure()
plt.loglog(num_procs, timings[0])
plt.title("OpenMC "+mode+"scaling")
plt.xlabel("Number of processors (-)")
plt.ylabel("Active batch time (s)")
plt.tight_layout()
plt.savefig("openmc_scaling_" + mode + "_time")

plt.figure()
plt.loglog(num_procs, np.array(num_procs) * n_active_batches * n_neutrons / timings[0])
plt.title("OpenMC scaling: "+mode)
plt.xlabel("Number of processors (-)")
plt.ylabel("Active neutrons per s (-)")
plt.tight_layout()
plt.savefig("openmc_scaling_" + mode + "_neutrons")
