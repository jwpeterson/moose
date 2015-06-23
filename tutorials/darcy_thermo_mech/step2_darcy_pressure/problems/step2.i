[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 0.304 # Length of test chamber
  ymax = 0.0257 # Test chamber radius
[]

[Variables]
  [./pressure]
    # Scaling this example up helps with convergence issues due to the small permeability
    # scaling = 1e6
  [../]
[]

[Kernels]
  [./darcy_pressure]
    type = DarcyPressure
    variable = pressure
    # permeability = 0.8451e-9 # (m^2) 1mm balls.  From paper
    # The DarcyPressure Kernel requires both permeability and viscosity parameters.
    permeability = 1
    viscosity = 1
  [../]
[]

[BCs]
  [./inlet]
    type = DirichletBC
    # type = PresetBC
    variable = pressure
    boundary = left
    # value = 4000 # (Pa) From Figure 2 from paper.  First data point for 1mm balls.
    value = 1
  [../]
  [./outlet]
    type = DirichletBC
    # type = PresetBC
    variable = pressure
    boundary = right
    value = 0 # (Pa) Gives the correct pressure drop from Figure 2 for 1mm balls
  [../]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = X
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  l_tol = 1.e-10
[]

[Outputs]
  output_initial = true
  exodus = true
  print_perf_log = true
  print_linear_residuals = true
[]
