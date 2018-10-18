# This is a modification of Yaqi's original problem that simply solves du/dt = t
# without involving stateful AuxVariables/Postprocessors.
[Mesh]
  type = GeneratedMesh
  dim  = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx   = 2
  ny   = 2
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./forcing_fn]
    type = ParsedFunction
    value = 't'
  [../]

  [./exact_fn]
    type = ParsedFunction
    value = '0.5*t*t'
  [../]
[]

[Kernels]
  [./ie]
    type = TimeDerivative
    variable = u
  [../]

  [./cv]
    type = BodyForce
    variable = u
    function = forcing_fn
  [../]
[]

[Postprocessors]
  # Estimate spatial norm of error at fixed time, ||e||_{L2}
  [./l2_err]
    type = ElementL2Error
    variable = u
    function = exact_fn
  [../]
  # Estimate spacetime norm ||e||_{L2, \infty}
  [./max_l2_err]
    type = TimeExtremeValue
    value_type = max
    postprocessor = l2_err
  [../]
  # Estimate spacetime norm ||e||_{L2, L1}
  [./cumulative_l2_err]
    type = CumulativeValuePostprocessor
    postprocessor = l2_err
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time   = 1.0
  # dt = 1.0
  # dt = 0.5
  # dt = 0.25
  # dt = 0.125
  # dt = 0.0625
  dt = 0.03125
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12

  # Use to set different TimeIntegrators
  [./TimeIntegrator]
    # type = ImplicitEuler
    # type = CrankNicolson
    type = LStableDirk2
  [../]
[]

[Outputs]
  csv = true
[]
