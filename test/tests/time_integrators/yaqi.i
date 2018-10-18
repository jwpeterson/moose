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

[AuxVariables]
  [./v]
  [../]
  [./c]
    initial_condition = 1
  [../]
[]

[Kernels]
  [./ie]
    type = TimeDerivative
    variable = u
  [../]

  [./cv]
    type = CoupledForce
    variable = u
    v = v
  [../]
[]

[AuxKernels]
  [./vt]
    type = VariableTimeIntegrationAux
    variable = v
    variable_to_integrate = c
  [../]
[]

[Postprocessors]
  [./intu]
    type = ElementIntegralVariablePostprocessor
    variable = u
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time   = 1.0
  dt = 0.1
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
[]

[Outputs]
  csv = true
[]
