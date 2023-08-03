# ==============================================================================
# GEOMETRY AND MESH
# ==============================================================================

# Full core has repeated 10010
# Quarter core has 10008 10009 10010
fuel_16 = 10010
fuel_24 = 10010
fuel_31 = 10010

[Mesh]
  [core]
    type = FileMeshGenerator
    file = 'save_mesh/2d_full_core.e'
    # file = 'quarter_core_2d.e'
  []
  [boundaries6]
    type = SideSetsFromNormalsGenerator
    input = core
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0'
    fixed_normal = true
    new_boundary = 'left right bottom top'
  []
  [shift_z]
    type = TransformGenerator
    input = 'boundaries6'
    transform = 'TRANSLATE'
    vector_value = '0 0 100'
  []

  # [Partitioner]
  #   type = HierarchicalGridPartitioner
  #   nx_nodes = 2
  #   ny_nodes = 5
  #   nz_nodes = 1

  #   # Does not work for 2D systems
  #   # 48 processors per node on sawtooth
  #   nx_procs = 1
  #   ny_procs = 1
  #   nz_procs = 1
  # []

  # Shared memory run
  # [Partitioner]
  #    type = GridPartitioner
  #    nx = 4
  #    ny = 1
  #    nz = 1
  #  []
[]

# ==============================================================================
# SET UP OPENMC SIMULATION IN MOOSE
# ==============================================================================

[Problem]
  solve = false
  kernel_coverage_check = false
  skip_nl_system_check = true
  verbose_multiapps = true
[]

# Main things we care about for the coupling
[AuxVariables]
  [temperature]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 300
  []
  [power]
    order = CONSTANT
    family = MONOMIAL
    block = '${fuel_16} ${fuel_24} ${fuel_31}'
  []
[]

[RayKernels/collision]
  type = CollisionKernel
  temperature = temperature
  # mesh block ids
  # blocks = "10000 10001 10004 10006 10008 10009 10010 10013 10014 10016"  # 2D core
  blocks = "10000 10001 2 3 ${fuel_16} ${fuel_24} ${fuel_31} 10013 4 10016 10004 10006 10014"  # 2D core - new mgs
  # blocks = "1 2 3 4 5 6 7 8 9"
  # openmc material id
  materials = "1      2 5 7 9      10     11    14  15   16     5     7  15"  # 2D core
  # materials = "14 17 1 7 9 10 11 16 15"
  # verbose = true
  z_coord = 100
[]

[RayBCs]
  # [reflect]
  #   type = ReflectRayBC
  #   boundary = '3 4' #'back front top right left bottom'
  # []
  # [vacuum]
  #   type = KillRayBC
  #   boundary = '1 2' #'bottom top right front'
  # []
  [vacuum]
    type = KillRayBC
    boundary = 'bottom left right top'
  []
[]

[UserObjects]
  [study]
    type = OpenMCStudy

    execute_on = TIMESTEP_END

    # Needed to cache trace information for RayTracingMeshOutput
    always_cache_traces = false
    segments_on_cache_traces = false

    # Needed to cache Ray data for RayTracingMeshOutput
    data_on_cache_traces = false
    aux_data_on_cache_traces = false

    # Parameters to make it work for now
    tolerate_failure = true
  []
[]

[Executioner]
  type = Transient
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]


# ==============================================================================
# TALLIES
# ==============================================================================

[UserObjects]
  # [tally]
  #   type = OpenMCTally
  #   id = 3
  #   particle_type = 'neutron'
  #   estimator = 'COLLISION'
  #   scores = 'flux scatter (n,fission) 16'
  #   filters = 'energy particle'
  #   energy_bins = '1e-5 1e3 2e7'
  #   # filter_ids = 3
  #   execute_on = 'initial'
  # []

  # [univtally]
  #   type = OpenMCTally
  #   id = 2
  #   particle_type = 'neutron'
  #   estimator = 'COLLISION'
  #   scores = 'kappa-fission'
  #   filters = 'universe'
  #   # filter_ids = 2
  #   execute_on = 'initial'
  # []

  [celltally]
    type = OpenMCTally
    id = 6
    particle_type = 'neutron'
    estimator = 'COLLISION'
    scores = 'kappa-fission'
    filters = 'cell'
    # filter_ids = 1
    execute_on = 'initial'
  []
[]

[AuxKernels]
  [cell_val]
    type = OpenMCTallyAux
    tally_id = 6
    execute_on = TIMESTEP_END
    variable = power
    granularity = 'cell'
    score = 'kappa-fission'
  []
[]

# Output on a pincell mesh
[VectorPostprocessors]
  [pin_powers]
    type = NearestPointIntegralVariablePostprocessor
    variable = 'power'
    block = '${fuel_16} ${fuel_24} ${fuel_31}'
    points_file = pin_file
    execute_on = 'TIMESTEP_END'
  []
[]

[MultiApps]
  [pin_mesh]
    type = TransientMultiApp
    execute_on = TIMESTEP_BEGIN
    input_files = 'output_pin.i'
  []
[]

[Transfers]
  [power_uo]
    type = MultiAppUserObjectTransfer
    to_multi_app = pin_mesh
    user_object = pin_powers
    variable = pin_power
    execute_on = TIMESTEP_END
  []
[]

# ==============================================================================
# SIMULATION PERFORMANCE STUDY
# ==============================================================================

# To look at domain decomposition
# [AuxVariables/domains]
#   order = CONSTANT
#   family = MONOMIAL
# []

# [AuxKernels]
#   [domains]
#     type = ProcessorIDAux
#     variable = domains
#   []
# []

[VectorPostprocessors]
  [mem]
    type = VectorMemoryUsage
    execute_on = 'INITIAL TIMESTEP_END'
    report_peak_value = true
    mem_units = megabytes
  []
[]
