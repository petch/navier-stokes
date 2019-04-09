__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.basesolver import *

class Decoupled(BaseSolver):
    def solve_unconstrained(self):
        problem = self.problem
        
        self.V = FunctionSpace(problem.mesh, self.ufe)
        self.Q = FunctionSpace(problem.mesh, self.pfe)
        self.u = TrialFunction(self.V)
        self.p = TrialFunction(self.Q)
        self.v = TestFunction(self.V)
        self.q = TestFunction(self.Q)
        self.u_ = Function(self.V)
        self.p_ = Function(self.Q)
        self.u0 = Function(self.V)
        self.p0 = Function(self.Q)
        self.du = Function(self.V)
        self.dp = Function(self.Q)
        self.eu = Function(self.V)
        self.ep = Function(self.Q)

        a1, L1, a2, L2, a3, L3 = self.forms()
        bcu = problem.dirichlet_u(self.V) 
        bcp = problem.dirichlet_p(self.Q)

        ue = project(problem.ue, self.V)
        pe = project(problem.pe, self.Q)
        if not self.stationary:
            assign(self.u_, ue)
            assign(self.p_, pe)
            assign(self.u0, ue)
            assign(self.p0, pe)
        
        def step():
            solve(a1 == L1, self.u_, bcu, solver_parameters=self.solver_parameters)
            solve(a2 == L2, self.p_, bcp, solver_parameters=self.solver_parameters)
            solve(a3 == L3, self.u_,      solver_parameters=self.solver_parameters)

            self.du.assign(self.u_ - self.u0)
            self.dp.assign(self.p_ - self.p0)
            self.u0.assign(self.u_)
            self.p0.assign(self.p_)

            if not self.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
            self.eu.assign(ue - self.u_)
            self.ep.assign(pe - self.p_)

        self.simulate(step)

    def solve_constrained(self):
        problem = self.problem
        
        self.V = FunctionSpace(problem.mesh, self.ufe)
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)
        self.u_ = Function(self.V)
        self.u0 = Function(self.V)
        self.du = Function(self.V)
        self.eu = Function(self.V)

        self.W = FunctionSpace(problem.mesh, MixedElement([self.pfe, self.rfe]))
        self.w = TrialFunction(self.W)
        self.z = TestFunction(self.W)
        self.w_ = Function(self.W)
        self.w0 = Function(self.W)
        self.dw = Function(self.W)
        self.ew = Function(self.W)

        self.Q,  self.R  = self.W.sub(0), self.W.sub(1)
        self.p,  self.c  = split(self.w)
        self.q,  self.r  = split(self.z)
        self.p_, self.c_ = split(self.w_)
        self.p0, self.c0 = split(self.w0)
        self.dp, self.dc = self.dw.split()
        self.ep, self.ec = self.ew.split(True)

        a1, L1, a2, L2, a3, L3 = self.forms()
        self.p_, self.c_ = self.w_.split() # only after self.forms()
        bcu = problem.dirichlet_u(self.V)
        bcp = problem.dirichlet_p(self.Q)

        ue = project(problem.ue, self.V)
        we = Function(self.W)
        pe = project(problem.pe, self.Q.collapse())
        ce = project(problem.ce, self.R.collapse())
        assign(we, [pe, ce])
        if not self.stationary:
            assign(self.u_, ue)
            assign(self.u0, ue)
            assign(self.w_, we)
            assign(self.w0, we)

        def step():
            solve(a1 == L1, self.u_, bcu, solver_parameters=self.solver_parameters)
            solve(a2 == L2, self.w_, bcp, solver_parameters=self.solver_parameters)
            solve(a3 == L3, self.u_,      solver_parameters=self.solver_parameters)

            self.du.assign(self.u_ - self.u0)
            self.dw.assign(self.w_ - self.w0)
            self.u0.assign(self.u_)
            self.w0.assign(self.w_)

            if not self.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
                problem.ce.t = problem.t
                ce.assign(problem.ce)
                assign(we, [pe, ce])
            self.eu.assign(ue - self.u_)
            self.ew.assign(we - self.w_)
            self.ep.assign(self.ew.cpp_object().sub(0))

        self.simulate(step)
