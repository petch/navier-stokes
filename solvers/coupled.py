__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.basesolver import *

class Coupled(BaseSolver):
    def solve(self):
        if len(self.problem.pd) > 0:
            self.solve_unconstrained()
        else:
            self.solve_constrained()

    def solve_unconstrained(self):
        problem = self.problem

        me = MixedElement([self.ufe, self.pfe])
        self.W = FunctionSpace(problem.mesh, me)
        self.V = self.W.sub(0).collapse()
        self.Q = self.W.sub(1).collapse()
        self.w = TrialFunction(self.W)
        self.z = TestFunction(self.W)
        self.w_ = Function(self.W)
        self.wn = Function(self.W)
        self.dw = Function(self.W)
        self.e  = Function(self.W)
        self.u,  self.p  = split(self.w)
        self.v,  self.q  = split(self.z)
        self.u_, self.p_ = split(self.w_)
        self.un, self.pn = split(self.wn)
        self.du, self.dp = self.dw.split()
        self.eu, self.ep = self.e.split(True)

        bcs = []
        for bc in problem.ud:
            space = self.W.sub(0) if len(bc) < 3 else self.W.sub(0).sub(bc[2])
            bcs.append(DirichletBC(space, bc[1], problem.bounds, bc[0]))
        for bc in problem.pd:
            bcs.append(DirichletBC(self.W.sub(1), bc[1], problem.bounds, bc[0]))

        lhs, rhs = self.forms()
        self.u_, self.p_ = self.w_.split() # only after self.forms()
        we = Function(self.W)
        ue = project(problem.ue, self.V)
        pe = project(problem.pe, self.Q)
        assign(we, [ue, pe])
        if not self.stationary:
            assign(self.w_, we)
            assign(self.wn, we)
        
        def step():
            solve(lhs == rhs, self.w_, bcs, solver_parameters=self.solver_parameters)
            self.dw.assign(self.w_ - self.wn)
            self.wn.assign(self.w_)

            if not self.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
                assign(we, [ue, pe, ce])
            self.e.assign(we - self.w_)
            self.eu.assign(self.e.cpp_object().sub(0))
            self.ep.assign(self.e.cpp_object().sub(1))
            
        self.simulate(step)

    def solve_constrained(self):
        problem = self.problem

        me = MixedElement([self.ufe, self.pfe, self.rfe])
        self.W = FunctionSpace(problem.mesh, me)
        self.V = self.W.sub(0).collapse()
        self.Q = self.W.sub(1).collapse()
        self.R = self.W.sub(2).collapse()
        self.w = TrialFunction(self.W)
        self.z = TestFunction(self.W)
        self.w_ = Function(self.W)
        self.wn = Function(self.W)
        self.dw = Function(self.W)
        self.e  = Function(self.W)
        self.u,  self.p,  self.c  = split(self.w)
        self.v,  self.q,  self.r  = split(self.z)
        self.u_, self.p_, self.c_ = split(self.w_)
        self.un, self.pn, self.cn = split(self.wn)
        self.du, self.dp, self.dc = self.dw.split()
        self.eu, self.ep, self.ec = self.e.split(True)

        bcs = []
        for bc in problem.ud:
            if len(bc) < 3:
                bcs.append(DirichletBC(self.W.sub(0), bc[1], problem.bounds, bc[0]))
            else:
                bcs.append(DirichletBC(self.W.sub(0).sub(bc[2]), bc[1], problem.bounds, bc[0]))
        for bc in problem.pd:
            bcs.append(DirichletBC(self.W.sub(1), bc[1], problem.bounds, bc[0]))

        lhs, rhs = self.forms()
        self.u_, self.p_, self.c_ = self.w_.split() # only after self.forms()
        we = Function(self.W)
        ue = project(problem.ue, self.V)
        pe = project(problem.pe, self.Q)
        ce = project(problem.ce, self.R)
        assign(we, [ue, pe, ce])
        if not self.stationary:
            assign(self.w_, we)
            assign(self.wn, we)
        
        def step():
            solve(lhs == rhs, self.w_, bcs, solver_parameters=self.solver_parameters)
            self.dw.assign(self.w_ - self.wn)
            self.wn.assign(self.w_)

            if not self.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
                problem.ce.t = problem.t
                ce.assign(problem.ce)
                assign(we, [ue, pe, ce])
            self.e.assign(we - self.w_)
            self.eu.assign(self.e.cpp_object().sub(0))
            self.ep.assign(self.e.cpp_object().sub(1))
            
        self.simulate(step)