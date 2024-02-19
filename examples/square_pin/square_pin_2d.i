[Mesh]
  [gmg]
    type = CartesianMeshGenerator
    dim = 2
    ix = '1 1 1'
    iy = '1 1 1'
    dx = '4 2 4'
    dy = '4 2 4'
  []
  [center]
    type = TransformGenerator
    input = 'gmg'
    transform = 'TRANSLATE'
    vector_value = '-5 -5 0'
  []
  # [gmg]
  #   type = GeneratedMeshGenerator
  #   dim = 2
  #   nx = 5
  #   ny = 5
  #   xmin = -5
  #   ymin = -5
  #   xmax = 5
  #   ymax = 5
  # []
  [add_infinite_z_pin]
    input = center
    type = SubdomainBoundingBoxGenerator
    top_right = '1 1 100'
    bottom_left = '-1 -1 -100'
    block_id = 1
    block_name = 'pin'
  []
[]

[Problem]
  solve = false
  kernel_coverage_check = false
[]

# Main things we care about for the coupling
[AuxVariables/temperature]
  order = CONSTANT
  family = MONOMIAL
  # This is the temperature default for OpenMC
  initial_condition = 293.6
[]

[AuxVariables/power]
[]

[RayKernels/collision]
  type = CollisionKernel
  temperature = temperature
  blocks = "0 1"
  # moderator fuel
  materials = "20 19" # openmc material id
  # verbose = true
  z_coord = 0
[]

[RayBCs]
  [reflect]
    type = ReflectRayBC
    boundary = 'top right left bottom'
  []
[]

[UserObjects/study]
  type = OpenMCStudy

  execute_on = TIMESTEP_END

  verbose = false

  # Needed to cache trace information for RayTracingMeshOutput
  always_cache_traces = false
  segments_on_cache_traces = false

  # Needed to cache Ray data for RayTracingMeshOutput
  data_on_cache_traces = false
  aux_data_on_cache_traces = false

  tolerate_failure = true
  verify_rays = false
  verify_trace_intersections = false
[]

[Executioner]
  type = Transient
[]

[Outputs]
  exodus = true
  csv = true
  [console]
    type = Console
    interval = 10
    hide = 'total_num_rays max_proc_m'
  []
  perf_graph = true
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

# To measure performance
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
  [neutrons_per_s]
    type = ParsedPostprocessor
    pp_names = 'total_num_rays total_time'
    function = 'total_num_rays / total_time'
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
