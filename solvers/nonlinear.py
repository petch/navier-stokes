__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.coupled import *
from meshes import Bound

class Nonlinear(Coupled):
    def __init__(self, problem_data, **params):
        super(Nonlinear, self).__init__(problem_data, **params)
        self.solver_parameters = {'newton_solver': self.solver_parameters}

    def forms(self):
        problem = self.problem
        tau, rho, mu, f, n, ds = problem.tau, problem.rho, problem.mu, problem.f, problem.n, problem.ds
        u_, p_, v, q, u0, p0 = self.u_, self.p_, self.v, self.q, self.u0, self.p0
        
        F = rho*dot(problem.convection(u_, u_), v)*dx \
          + inner(problem.stress(u_, p_), sym(nabla_grad(v)))*dx \
          - dot(f, v)*dx \
          + dot(div(u_), q)*dx

        for i, dun, dut in problem.un:
            F -= dot(problem.stress_n(u_, p_, n, dun, dut), v)*ds(int(i))

        if problem.constrained():
            ce, c_, r = problem.ce, self.c_, self.r
            F += (p_ - ce)*r*dx + c_*q*dx

        if not self.stationary:
            F += rho/tau*dot(u_ - u0, v)*dx

        return F, 0
