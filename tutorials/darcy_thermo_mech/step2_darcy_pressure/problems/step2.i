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
    # scaling = 1e4
  [../]
[]

[Kernels]
  [./darcy_pressure]
    type = DarcyPressure
    variable = pressure
    permeability = 0.8451e-9 # (m^2) 1mm balls.  From paper
  [../]
[]

[BCs]
  [./inlet]
    type = DirichletBC
    variable = pressure
    boundary = left
    value = 4000 # (Pa) From Figure 2 from paper.  First data point for 1mm balls.
  [../]
  [./outlet]
    type = DirichletBC
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

  # The linear and nonlinear residuals agree perfectly if we use
  # solve_type = NEWTON for this example. So is the issue something
  # with computing the action of the Jacobian slightly incorrectly
  # with PJFNK, perhaps due to finite differencing parameters?
  # Variable scaling is also unnecessary if we use solve_type = NEWTON.
  solve_type = PJFNK
  # solve_type = NEWTON

  # PETSc options that control the finite differencing algorithm and
  # -mat_mffd_type <(null)>: Matrix free type (one of) ds wp (MatMFFDSetType)
  # -mat_mffd_err <1.49012e-08>: set sqrt relative error in function (MatMFFDSetFunctionError)
  petsc_options_iname = '-pc_type -pc_hypre_type -mat_mffd_err'
  petsc_options_value = 'hypre boomeramg 1.0'

  # Convergence for -mat_mffd_type wp is the same as if we don't specify:
  # 0 Nonlinear |R| = 1.326650e+04
  #       0 Linear |R| = 1.326650e+04
  #       1 Linear |R| = 2.944220e-05
  #       2 Linear |R| = 1.038678e-06
  #       3 Linear |R| = 3.625419e-08
  #       4 Linear |R| = 2.946512e-09
  # 1 Nonlinear |R| = 1.864784e+00
  #       0 Linear |R| = 1.864784e+00
  #       1 Linear |R| = 4.125705e-09
  #       2 Linear |R| = 1.545564e-10
  #       3 Linear |R| = 5.822819e-12
  #       4 Linear |R| = 4.933235e-13
  # 2 Nonlinear |R| = 7.953179e-07

  # Convergence while setting -mat_mffd_err 1.e-5
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944635e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625930e-08
  #      4 Linear |R| = 2.946928e-09
  # 1 Nonlinear |R| = 4.481955e-03
  #      0 Linear |R| = 4.481955e-03
  #      1 Linear |R| = 3.859683e-10
  #      2 Linear |R| = 1.935739e-11
  #      3 Linear |R| = 9.429797e-13
  #      4 Linear |R| = 6.275555e-14
  #      5 Linear |R| = 3.090261e-15
  # 2 Nonlinear |R| = 4.826324e-10

  # Convergence while setting -mat_mffd_err 1.e-4
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944634e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625928e-08
  #      4 Linear |R| = 2.946927e-09
  # 1 Nonlinear |R| = 3.926530e-04
  #      0 Linear |R| = 3.926530e-04
  #      1 Linear |R| = 3.849709e-10
  #      2 Linear |R| = 1.909964e-11
  #      3 Linear |R| = 9.714352e-13
  #      4 Linear |R| = 6.189219e-14
  #      5 Linear |R| = 2.975686e-15
  #      6 Linear |R| = 1.944565e-16
  # 2 Nonlinear |R| = 1.176417e-10

  # Convergence while setting -mat_mffd_err 1.e-3
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944634e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625929e-08
  #      4 Linear |R| = 2.946927e-09
  # 1 Nonlinear |R| = 2.684445e-06

  # Convergence while setting -mat_mffd_err 1.e-2
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944634e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625929e-08
  #      4 Linear |R| = 2.946927e-09
  # 1 Nonlinear |R| = 2.684445e-06

  # Convergence while setting -mat_mffd_err 1.e-1
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944634e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625929e-08
  #      4 Linear |R| = 2.946926e-09
  # 1 Nonlinear |R| = 2.403328e-07

  # Convergence while setting -mat_mffd_err 1.0 (best one yet!)
  # 0 Nonlinear |R| = 1.326650e+04
  #      0 Linear |R| = 1.326650e+04
  #      1 Linear |R| = 2.944634e-05
  #      2 Linear |R| = 1.038824e-06
  #      3 Linear |R| = 3.625929e-08
  #      4 Linear |R| = 2.946926e-09
  # 1 Nonlinear |R| = 4.542393e-08

  l_tol = 1.e-12
[]

[Outputs]
  output_initial = true
  exodus = true
  print_perf_log = true
  print_linear_residuals = true
[]
