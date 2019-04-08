__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.coupled import *
from meshes import Bound

class Nonlinear(Coupled):
    def __init__(self, problem, **params):
        super(Nonlinear, self).__init__(problem, **params)
        self.solver_parameters = {'newton_solver': self.solver_parameters}

    def forms(self):
        problem = self.problem
        tau, rho, mu, f, n, ds = problem.tau, problem.rho, problem.mu, problem.f, problem.n, problem.ds
        u_, p_, v, q, un, pn = self.u_, self.p_, self.v, self.q, self.un, self.pn
        
        F = rho*dot(problem.convection(u_, u_), v)*dx \
          + inner(problem.stress(u_, p_), sym(nabla_grad(v)))*dx \
          - dot(f, v)*dx \
          + dot(div(u_), q)*dx
        #   - dot(dot(problem.stress(u_, p_), n), v)*ds \
        if not self.stationary:
            F += rho/tau*dot(u_ - un, v)*dx
        for i, g in problem.un:
            F -= dot(problem.stress_n(u_, p_, n, g), v)*ds(int(i))
        if len(problem.pd) == 0:
            ce, c_, r = problem.ce, self.c_, self.r
            F += (p_ - ce)*r*dx + c_*q*dx

        return F, 0
