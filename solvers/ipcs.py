__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.decoupled import *

class IPCS(Decoupled):
    def forms(self):
        problem = self.problem
        rho, tau, mu, f, n, ds = problem.rho, problem.tau, problem.mu, problem.f, problem.n, problem.ds
        u, p, v, q, u_, p_, u0, p0 = self.u, self.p, self.v, self.q, self.u_, self.p_, self.u0, self.p0
        U = (u + u0)/2

        F1 = rho/tau*dot(u - u0, v)*dx \
           + rho*dot(problem.convection(u0, u0), v)*dx \
           + inner(problem.stress(U, p0), sym(nabla_grad(v)))*dx \
           - dot(f, v)*dx

        for i, dun, dut in problem.un:
            F1 -= dot(problem.stress_n(U, p0, n, dun, dut), v)*ds(int(i))

        F2 = dot(nabla_grad(p - p0), nabla_grad(q))*dx + rho/tau*div(u_)*q*dx
        
        if problem.constrained():
            ce, c, r = problem.ce, self.c, self.r
            F2 += (p - ce)*r*dx + c*q*dx
        
        F3 = dot(u - u_, v)*dx + tau/rho*dot(nabla_grad(p_ - p0), v)*dx

        return lhs(F1), rhs(F1), lhs(F2), rhs(F2), lhs(F3), rhs(F3)
        
