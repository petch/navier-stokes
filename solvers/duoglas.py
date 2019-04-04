__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.decoupled import *

# Duoglas-Rachford
class Duoglas(Decoupled):
    def forms(self):
        problem = self.problem
        rho, tau, mu, f, n = problem.rho, problem.tau, problem.mu, problem.f, problem.n
        u, p, v, q, u_, p_, un, pn = self.u, self.p, self.v, self.q, self.u_, self.p_, self.un, self.pn 

        F1 = rho/tau*dot(u - un, v)*dx \
           + rho*dot(problem.convection(un, u), v)*dx \
           + inner(problem.stress(u, pn), sym(nabla_grad(v)))*dx \
           - dot(dot(problem.stress(u, pn), n), v)*ds \
           - dot(f, v)*dx
        F2 = dot(nabla_grad(p - pn), nabla_grad(q))*dx + rho/tau*div(u_)*q*dx
        F3 = dot(u - u_, v)*dx + tau/rho*dot(nabla_grad(p_ - pn), v)*dx

        return lhs(F1), rhs(F1), lhs(F2), rhs(F2), lhs(F3), rhs(F3)
