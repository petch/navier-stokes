__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from solvers.basesolver import *

class Decoupled(BaseSolver):
    def solve(self):
        problem = self.problem
        
        self.V = FunctionSpace(problem.mesh, self.ufe)
        self.Q = FunctionSpace(problem.mesh, self.pfe)
        self.R = FunctionSpace(problem.mesh, self.rfe)
        self.u = TrialFunction(self.V)
        self.p = TrialFunction(self.Q)
        self.v = TestFunction(self.V)
        self.q = TestFunction(self.Q)
        self.u_ = Function(self.V)
        self.p_ = Function(self.Q)
        self.un = Function(self.V)
        self.pn = Function(self.Q)
        self.du = Function(self.V)
        self.dp = Function(self.Q)
        self.eu = Function(self.V)
        self.ep = Function(self.Q)

        bcu = []
        for bc in problem.ud:
            bcu.append(DirichletBC(self.V, bc[0], bc[1]))
        bcp = []
        for bc in problem.pd:
            bcp.append(DirichletBC(self.Q, bc[0], bc[1]))

        a1, L1, a2, L2, a3, L3 = self.forms()
        ue = project(problem.ue, self.V)
        pe = project(problem.pe, self.Q)
        ce = project(problem.ce, self.R)
        if not self.stationary:
            assign(self.u_, ue)
            assign(self.p_, pe)
            assign(self.un, ue)
            assign(self.pn, pe)
        
        def step():
            solve(a1 == L1, self.u_, bcu, solver_parameters=solver_parameters)
            solve(a2 == L2, self.p_, bcp, solver_parameters=solver_parameters)
            solve(a3 == L3, self.u_,      solver_parameters=solver_parameters)

            self.p_.vector()[:] -= assemble(self.p_*dx) - assemble(ce*dx)

            self.du.assign(self.u_ - self.un)
            self.dp.assign(self.p_ - self.pn)
            self.un.assign(self.u_)
            self.pn.assign(self.p_)

            if not self.stationary:
                problem.ue.t = problem.t
                ue.assign(problem.ue)
                problem.pe.t = problem.t
                pe.assign(problem.pe)
            self.eu.assign(ue - self.u_)
            self.ep.assign(pe - self.p_)

        self.simulate(step)
