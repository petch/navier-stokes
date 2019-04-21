from .basesolver import *

class Coupled(BaseSolver):
    def solve_unconstrained(self):
        problem = self.problem

        self.W = FunctionSpace(problem.mesh, MixedElement([self.ufe, self.pfe]))
        self.w = TrialFunction(self.W)
        self.z = TestFunction(self.W)
        self.w_ = Function(self.W)
        self.w0 = Function(self.W)
        self.dw = Function(self.W)
        self.ew = Function(self.W)

        self.V,  self.Q  = self.W.sub(0), self.W.sub(1)
        self.u,  self.p  = split(self.w)
        self.v,  self.q  = split(self.z)
        self.u_, self.p_ = split(self.w_)
        self.u0, self.p0 = split(self.w0)
        self.du, self.dp = self.dw.split()
        self.eu, self.ep = self.ew.split(True)

        lhs, rhs = self.forms()
        self.u_, self.p_ = self.w_.split() # only after self.forms()
        bcs = problem.dirichlet_u(self.V) + problem.dirichlet_p(self.Q)
        
        we = Function(self.W)
        ue = project(problem.ue, self.V.collapse())
        pe = project(problem.pe, self.Q.collapse())
        assign(we, [ue, pe])
        if not self.problem.stationary:
            assign(self.w_, we)
            assign(self.w0, we)
        
        def step():
            solve(lhs == rhs, self.w_, bcs, solver_parameters=self.solver_parameters)
            self.dw.assign(self.w_ - self.w0)
            self.w0.assign(self.w_)

            if not self.problem.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
                assign(we, [ue, pe, ce])
            self.ew.assign(we - self.w_)
            self.eu.assign(self.ew.cpp_object().sub(0))
            self.ep.assign(self.ew.cpp_object().sub(1))
            
        self.simulate(step)

    def solve_constrained(self):
        problem = self.problem

        self.W = FunctionSpace(problem.mesh, MixedElement([self.ufe, self.pfe, self.rfe]))
        self.w = TrialFunction(self.W)
        self.z = TestFunction(self.W)
        self.w_ = Function(self.W)
        self.w0 = Function(self.W)
        self.dw = Function(self.W)
        self.ew = Function(self.W)

        self.V,  self.Q,  self.R  = self.W.sub(0), self.W.sub(1), self.W.sub(2)
        self.u,  self.p,  self.c  = split(self.w)
        self.v,  self.q,  self.r  = split(self.z)
        self.u_, self.p_, self.c_ = split(self.w_)
        self.u0, self.p0, self.c0 = split(self.w0)
        self.du, self.dp, self.dc = self.dw.split()
        self.eu, self.ep, self.ec = self.ew.split(True)

        lhs, rhs = self.forms()
        self.u_, self.p_, self.c_ = self.w_.split() # only after self.forms()
        bcs = problem.dirichlet_u(self.V) + problem.dirichlet_p(self.Q)

        we = Function(self.W)
        ue = project(problem.ue, self.V.collapse())
        pe = project(problem.pe, self.Q.collapse())
        ce = project(problem.ce, self.R.collapse())
        assign(we, [ue, pe, ce])
        if not self.problem.stationary:
            assign(self.w_, we)
            assign(self.w0, we)
        
        def step():
            solve(lhs == rhs, self.w_, bcs, solver_parameters=self.solver_parameters)
            self.dw.assign(self.w_ - self.w0)
            self.w0.assign(self.w_)

            if not self.problem.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
                problem.ce.t = problem.t
                ce.assign(problem.ce)
                assign(we, [ue, pe, ce])
            self.ew.assign(we - self.w_)
            self.eu.assign(self.ew.cpp_object().sub(0))
            self.ep.assign(self.ew.cpp_object().sub(1))
            
        self.simulate(step)