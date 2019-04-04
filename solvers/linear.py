__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.coupled import *

class Linear(Coupled):
    def forms(self):
        problem = self.problem
        rho, tau, mu, f, n = problem.rho, problem.tau, problem.mu, problem.f, problem.n
        u, p, c, v, q, r, un, pn, cn = self.u, self.p, self.c, self.v, self.q, self.r, self.un, self.pn, self.cn

        F = rho/tau*dot(u - un, v)*dx \
          + rho*dot(problem.convection(u, un), v)*dx \
          + inner(problem.stress(u, p), sym(nabla_grad(v)))*dx \
          - dot(dot(problem.stress(u, p), n), v)*ds \
          - dot(f, v)*dx \
          + dot(div(u), q)*dx \
          + p*r*dx + c*q*dx

        return lhs(F), rhs(F)
        