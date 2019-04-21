from .coupled import *

class Nonlinear(Coupled):
    def __init__(self, problem_data, **params):
        super(Nonlinear, self).__init__(problem_data, **params)
        self.solver_parameters = {'newton_solver': self.solver_parameters}

    def forms(self):
        problem = self.problem
        dt, rho, mu, f = problem.dt, problem.rho, problem.mu, problem.f
        u_, p_, v, q, u0, p0 = self.u_, self.p_, self.v, self.q, self.u0, self.p0
        
        F = rho*dot(problem.convection(u_, u_), v)*dx \
          + inner(problem.stress(u_, p_), sym(nabla_grad(v)))*dx \
          + problem.neumann_u(u_, p_, v) \
          - dot(f, v)*dx \
          + dot(div(u_), q)*dx

        if problem.constrained():
            ce, c_, r = problem.ce, self.c_, self.r
            F += (p_ - ce)*r*dx + c_*q*dx

        if not self.problem.stationary:
            F += rho/dt*dot(u_ - u0, v)*dx

        return F, 0
