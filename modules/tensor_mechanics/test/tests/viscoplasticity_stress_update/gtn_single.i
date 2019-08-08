[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = 0.002
  ymax = 0.002
[]

[Modules/TensorMechanics/Master/All]
  strain = FINITE
  add_variables = true
  base_name = 'asdf'
  generate_output = 'strain_xx strain_yy strain_xy hydrostatic_stress vonmises_stress'
[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x = '0 0.1'
    y = '0 1e-5'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e10
    poissons_ratio = 0.3
    base_name = 'asdf'
  [../]
  [./strain]
    type = ComputeMultiplePorousInelasticStress
    inelastic_models = gtn
    initial_porosity = 0.1
    outputs = all
    base_name = 'asdf'
  [../]

  [./gtn]
    type = ViscoplasticityStressUpdate
    total_strain_base_name = 'asdf'
    coefficient = 'coef'
    power = 3
    base_name = 'gtn'
    viscoplasticity_model = GTN
    outputs = all
    relative_tolerance = 1e-11
  [../]
  [./coef]
    type = ParsedMaterial
    f_name = coef
    function = '1e-18 * exp(-4e4 / 1.987 / 1200)'
  [../]
[]

[BCs]
  [./no_disp_x]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./no_disp_y]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./pull_disp_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = pull
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.01
  end_time = 0.12
[]

[Postprocessors]
  [./disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  [../]
  [./disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./avg_hydro]
    type = ElementAverageValue
    variable = asdf_hydrostatic_stress
  [../]
  [./avg_vonmises]
    type = ElementAverageValue
    variable = asdf_vonmises_stress
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./num_lin]
    type = NumLinearIterations
    outputs = console
  [../]
  [./num_nonlin]
    type = NumNonlinearIterations
    outputs = console
  [../]
  [./eff_creep_strain]
    type = ElementAverageValue
    variable = gtn_effective_viscoplasticity
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = porosity
  [../]
[]

[Outputs]
  csv = true
[]
