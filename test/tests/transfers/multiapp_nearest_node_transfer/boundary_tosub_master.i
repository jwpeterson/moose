# The mesh for the master app consists of two disconnected 2D parts
# arranged vertically, the top part consists of 4 elements in a 2x2
# configuration while the bottom part is a single 2D element.  The
# subdomains (blocks) and boundaries are numbered as follows:
#
#          1
#  -----------------
# |    block 1      |
#  -----------------
#          2
#
#          3
#  -----------------
# |    block 2      |
#  -----------------
#          4
#
# and we solve the diffusion equation on each piece with Dirichlet BCs
# on each of the numbered boundaries. The true solution is therefore
# u = a*y + b in each piece, for appropriate constants a, b. The
# spatial coordinates of the two pieces are:
#
#            -----------------* (8,0)
#           |    block 1      |
#    (0,-1) *-----------------
#
#            -----------------* (8,-2)
#           |    block 2      |
#    (0,-3) *-----------------
#
#           o--------o--------o  <= 1D SubApp @ y=-4 with two elements
#
# There is a single 1D SubApp set up at the position y=-4 (as shown)
# which has 4 different values transferred to it, one corresponding to
# each of the boundaries 1-4. The actual values are controlled by
# the DirichletBCs, so we should get:
# from_master_1 = 4
# from_master_2 = 3
# from_master_3 = 2
# from_master_4 = 1
[Mesh]
  file = 2blk.e
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./pid]
    order = constant
    family = monomial
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxKernels]
  [./pid]
    type = ProcessorIDAux
    variable = pid
  [../]
[]

[BCs]
  [./top_1]
    type = DirichletBC
    variable = u
    boundary = '1'
    value = 4
  [../]
  [./top_2]
    type = DirichletBC
    variable = u
    boundary = '2'
    value = 3
  [../]

  [./bot_3]
    type = DirichletBC
    variable = u
    boundary = '3'
    value = 2
  [../]
  [./bot_4]
    type = DirichletBC
    variable = u
    boundary = '4'
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = MooseTestApp
    execute_on = timestep_end
    positions = '0 -4 0'
    input_files = boundary_tosub_sub.i
  [../]
[]

[Transfers]
  [./to_sub_1]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = u
    source_boundary = 1
    variable = from_master_1
  [../]
  [./to_sub_2]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = u
    source_boundary = 2
    variable = from_master_2
  [../]
  [./to_sub_3]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = u
    source_boundary = 3
    variable = from_master_3
  [../]
  [./to_sub_4]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = u
    source_boundary = 4
    variable = from_master_4
  [../]
[]
