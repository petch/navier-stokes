from .decoupled import *

# Duoglas-Rachford
class Duoglas(Decoupled):
    def forms(self):
        problem = self.problem
        rho, dt, mu, f = problem.rho, problem.dt, problem.mu, problem.f
        u, p, v, q, u_, p_, u0, p0 = self.u, self.p, self.v, self.q, self.u_, self.p_, self.u0, self.p0 

        F1 = rho/dt*dot(u - u0, v)*dx \
           + rho*dot(problem.convection(u0, u), v)*dx \
           + inner(problem.stress(u, p0), sym(nabla_grad(v)))*dx \
           + problem.neumann_u(u, p0, v) \
           - dot(f, v)*dx

        F2 = dot(nabla_grad(p - p0), nabla_grad(q))*dx + rho/dt*div(u_)*q*dx

        if problem.constrained():
            ce, c, r = problem.ce, self.c, self.r
            F2 += (p - ce)*r*dx + c*q*dx

        F3 = dot(u - u_, v)*dx + dt/rho*dot(nabla_grad(p_ - p0), v)*dx

        return lhs(F1), rhs(F1), lhs(F2), rhs(F2), lhs(F3), rhs(F3)
