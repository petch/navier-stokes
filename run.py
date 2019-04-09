from solvers import *
from problems import *
from meshes import *

params = dict(
    stationary = True,
    solver = Nonlinear,
    problem = PoiseuilleFlow,#ChannelFlow,ShiftedVortex,
    level = LogLevel.WARNING,
    mesh = Square(16),#ChannelThin(32),
    linear_solver = 'mumps',
    preconditioner = 'none',
)

# run(params, convection = ['divergence', 'advective', 'skewsymmetric', 'rotation'])
# run(params, solver = solvers_list)
# run(params, problem = problems_list)

run(params)

# bicgstab, cg, gmres, minres, richardson, tfqmr
# mumps, petsc, superlu, umfpack
# amg, hypre_amg, hypre_euclid, hypre_parasails, icc, ilu, jacobi, none, petsc_amg, sor  
# icc, jacobi, none, 