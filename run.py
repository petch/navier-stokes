from solvers import *
from problems import *
from meshes import *

run(dict(
        stationary = True,
        solver_class = Nonlinear,
        problem_class = ShiftedVortex,
        level = LogLevel.WARNING,
        mesh = UnitSquareMesh(16, 16, 'crossed'), #ChannelSymm(32),
        linear_solver = 'mumps',
        preconditioner = 'none',
        elements = ('CG', 2, 'DG', 0),
    ),
    # elements = [('CG', 2, 'DG', 0)]
    # convection = ['divergence', 'advective', 'skewsymmetric', 'rotation']
    # solver_class = [Nonlinear, Linear, IPCS, Duoglas],
)
# bicgstab, cg, gmres, minres, richardson, tfqmr
# mumps, petsc, superlu, umfpack
# amg, hypre_amg, hypre_euclid, hypre_parasails, icc, ilu, jacobi, none, petsc_amg, sor  
# icc, jacobi, none, 