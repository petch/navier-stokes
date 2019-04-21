from .coupled import *

class Linear(Coupled):
    def forms(self):
        problem = self.problem
        rho, dt, mu, f, n, ds = problem.rho, problem.dt, problem.mu, problem.f, problem.n, problem.ds
        u, p, v, q, u0, p0 = self.u, self.p, self.v, self.q, self.u0, self.p0

        F = rho/dt*dot(u - u0, v)*dx \
          + rho*dot(problem.convection(u, u0), v)*dx \
          + inner(problem.stress(u, p), sym(nabla_grad(v)))*dx \
          + problem.neumann_u(u, p, v) \
          - dot(f, v)*dx \
          + dot(div(u), q)*dx

        if problem.constrained():
            ce, c, r = problem.ce, self.c, self.r
            F += (p - ce)*r*dx + c*q*dx

        return lhs(F), rhs(F)
        