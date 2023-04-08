# Large box size means less intersections, maximizing speed
box_size = 10
# Number of elements must be above number of MPI ranks to have
# decent partitioning

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
    xmin = ${fparse -box_size / 2}
    ymin = ${fparse -box_size / 2}
    zmin = ${fparse -box_size / 2}
    xmax = ${fparse box_size / 2}
    ymax = ${fparse box_size / 2}
    zmax = ${fparse box_size / 2}
  []
[]

[Problem]
  solve = false
  kernel_coverage_check = false
[]

##############################################################################
# Problem definition

# Main things we care about for the coupling
[AuxVariables/temperature]
  order = CONSTANT
  family = MONOMIAL
  initial_condition = 300
[]

[RayKernels/collision]
  type = CollisionKernel
  temperature = temperature
  blocks = "0 1 2"
  materials = "1 1 1" # openmc material id
  # verbose = true
[]

[RayBCs]
  [reflect]
    type = ReflectRayBC
    boundary = 'back front top right left bottom'
  []
[]

[UserObjects/study]
  type = OpenMCStudy

  execute_on = TIMESTEP_END

  # Needed to cache trace information for RayTracingMeshOutput
  always_cache_traces = false
  segments_on_cache_traces = false

  # Needed to cache Ray data for RayTracingMeshOutput
  data_on_cache_traces = false
  aux_data_on_cache_traces = false

  # Dont kill simulation on ray tracing failure
  tolerate_failure = true
[]

[Executioner]
  type = Transient
[]

##############################################################################
# Outputs

[Outputs]
  exodus = true
  csv = true
  [console]
    type = Console
    interval = 10
    hide = 'total_num_rays max_proc_m'
  []
[]

# To look at domain decomposition
[AuxVariables/domain]
[]

[AuxKernels]
  [domains]
    type = ProcessorIDAux
    variable = domain
  []
[]

##############################################################################
# Performance analysis

[Postprocessors]
  [total_mem]
    type = MemoryUsage
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [proc_mem]
    type = MemoryUsage
    value_type = "average"
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_proc_m]
    type = MemoryUsage
    value_type = "max_process"
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_time]
    type = PerfGraphData
    execute_on = 'INITIAL TIMESTEP_END'
    data_type = 'TOTAL'
    section_name = 'Root'
  []
  [run_time]
    type = ChangeOverTimePostprocessor
    postprocessor = total_time
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Neutrons per second to compare to openmc
  [num_rays]
    type = VectorPostprocessorComponent
    vectorpostprocessor = per_proc_ray_tracing
    index = 0
    vector_name = rays_traced
  []
  [total_num_rays]
    type = CumulativeValuePostprocessor
    postprocessor = num_rays
  []
  [neutrons_per_s_tot]
    type = ParsedPostprocessor
    pp_names = 'total_num_rays total_time'
    function = 'total_num_rays / total_time'
  []
  [neutrons_per_s_last]
    type = ParsedPostprocessor
    pp_names = 'num_rays run_time'
    function = 'num_rays / run_time'
  []
[]

[VectorPostprocessors/per_proc_ray_tracing]
  type = PerProcessorRayTracingResultsVectorPostprocessor
  execute_on = TIMESTEP_END
  study = study
  outputs = none
[]

##############################################################################
# Tallies

[UserObjects]
  [tally]
    type = OpenMCTally
    particle_type = 'neutron'
    estimator = 'COLLISION'
    scores = 'kappa-fission'
    filters = 'cell'
    execute_on = 'initial'
    id = 1
  []
[]

[AuxKernels]
  [cell_val]
    type = OpenMCTallyAux
    granularity = 'cell'
    score = 'kappa-fission'
    tally_id = 1
    execute_on = TIMESTEP_END
    variable = power
  []
[]

[AuxVariables]
  [power]
    order = CONSTANT
    family = MONOMIAL
  []
[]

# Examine deviation
[Postprocessors]
  [average_power]
    type = ElementAverageValue
    variable = power
  []
[]
