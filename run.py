from solvers import *
from problems import *
from meshes import *

params = dict(
    stationary = True,
    solver_class = Nonlinear,
    problem_class = ChannelFlow,#ShiftedVortex, 
    level = LogLevel.WARNING,
    mesh = ChannelSymm(32),#UnitSquareMesh(16, 16, 'crossed'), #
    linear_solver = 'mumps',
    preconditioner = 'none'
)

# print('Test convection')
# run(params, convection = ['divergence', 'advective', 'skewsymmetric', 'rotation'])
# print('Test solvers')
# run(params, solver_class = solvers_list)
# print('Test problems')
# run(params, problem_class = problems_list)

run(params)

# bicgstab, cg, gmres, minres, richardson, tfqmr
# mumps, petsc, superlu, umfpack
# amg, hypre_amg, hypre_euclid, hypre_parasails, icc, ilu, jacobi, none, petsc_amg, sor  
# icc, jacobi, none, 