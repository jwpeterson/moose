# This input file tests several different things:
# .) The axisymmetric (RZ) form of the governing equations.
# .) An open boundary.
# .) Integrating the pressure by parts.
# .) Natural boundary condition at the outlet.
[GlobalParams]
  # mu=.5e-2  # Re=100
  mu=1e-3     # Re=500
  # mu=1      # Re=1/2
  rho=1
  gravity = '0 0 0'

  # Stabilization parameters
  supg = true
  pspg = true
  alpha = 1
[]

[Mesh]
  # file = '2d_cone.msh'
  file = 'cone_linear.e'
  # file = 'cone_quadratic.e'
  # This version of the quadratic mesh happens to not have
  # a Tri6 with all three vertices on the boundary, but this does not seem to have
  # any effect on the simulation... the only thing that matters is apparently the
  # Reynolds number.
  # file = 'cone_quadratic_qtri.e'
[]

[Problem]
  coord_type = RZ
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = Newton
  [../]
[]

[Executioner]
  type = Transient
  # dt = 0.005

  [./TimeStepper]
    dt = .005
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    growth_factor = 1.2
    optimal_iterations = 5
  [../]
  trans_ss_check = true
  ss_check_tol = 1e-10

  dtmin = 0.001
  num_steps = 1000
  l_max_its = 300

  # Note: The Steady executioner can be used for this problem, if you
  # drop the INSMomentumTimeDerivative kernels and use the following
  # direct solver options.
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'

  # Block Jacobi works well for this problem, as does "-pc_type asm
  # -pc_asm_overlap 2", but an overlap of 1 does not work for some
  # reason?
  # petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_levels'
  # petsc_options_value = 'bjacobi  ilu          4'

  # petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  # petsc_options_value = 'asm      2               ilu          3'

  # Set the linear tolerance dynamically based on Eisenstat-Walker formula. This is
  # only relevant when not using a direct solver. It generally requires more nonlinear
  # steps, so it may not be the best approach when an expensive preconditioner is being
  # used.
  petsc_options = '-snes_ksp_ew'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-14
  nl_max_its = 20
[]

[Outputs]
  csv = true
  console = true
  [./out]
    type = Exodus
  [../]
[]

[Variables]
  [./vel_x]
    # Velocity in radial (r) direction
    # order = SECOND
  [../]
  [./vel_y]
    # Velocity in axial (z) direction
    # order = SECOND
  [../]
  [./p]
  [../]
[]

[BCs]
  [./u_in]
    type = DirichletBC
    boundary = bottom
    variable = vel_x
    value = 0
  [../]
  [./v_in]
    type = FunctionDirichletBC
    boundary = bottom
    variable = vel_y
    function = 'inlet_func'
  [../]
  [./u_axis_and_walls]
    type = DirichletBC
    boundary = 'left right'
    variable = vel_x
    value = 0
  [../]
  [./v_no_slip]
    type = DirichletBC
    boundary = 'right'
    variable = vel_y
    value = 0
  [../]
[]


[Kernels]
  [./x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  [../]
  [./y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  [../]
  [./mass]
    type = INSMassRZ
    variable = p
    u = vel_x
    v = vel_y
    p = p
  [../]
  [./x_momentum_space]
    type = INSMomentumLaplaceFormRZ
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumLaplaceFormRZ
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '${GlobalParams/rho} ${GlobalParams/mu}'
  [../]
[]

[Functions]
  [./inlet_func]
    type = ParsedFunction
    value = '-4 * x^2 + 1'
  [../]
[]

[Postprocessors]
  [./flow_in]
    type = VolumetricFlowRate
    vel_x = vel_x
    vel_y = vel_y
    boundary = 'bottom'
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
  [./flow_out]
    type = VolumetricFlowRate
    vel_x = vel_x
    vel_y = vel_y
    boundary = 'top'
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
[]
