__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.coupled import *

class Linear(Coupled):
    def forms(self):
        problem = self.problem
        rho, tau, mu, f, n, ds = problem.rho, problem.tau, problem.mu, problem.f, problem.n, problem.ds
        u, p, v, q, u0, p0 = self.u, self.p, self.v, self.q, self.u0, self.p0

        F = rho/tau*dot(u - u0, v)*dx \
          + rho*dot(problem.convection(u, u0), v)*dx \
          + inner(problem.stress(u, p), sym(nabla_grad(v)))*dx \
          - dot(f, v)*dx \
          + dot(div(u), q)*dx

        for i, dun, dut in problem.un:
            F -= dot(problem.stress_n(u, p, n, dun, dut), v)*ds(int(i))

        if problem.constrained():
            ce, c, r = problem.ce, self.c, self.r
            F += (p - ce)*r*dx + c*q*dx

        return lhs(F), rhs(F)
        