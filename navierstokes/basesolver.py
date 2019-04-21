from dolfin import *
from enum import Enum
import matplotlib.pyplot as plt
import sys

class LinearSolver(Enum):
    default = 0
    mumps, petsc, superlu, umfpack = 1, 2, 3, 4
    bicgstab, cg, gmres, minres, richardson, tfqmr = 5, 6, 7, 8, 9, 10

class Preconditioner(Enum):
    default = 0
    none, icc, jacobi, ilu, sor = 1, 2, 3, 4, 5
    amg, hypre_amg, petsc_amg = 6, 7, 8
    hypre_euclid, hypre_parasails = 9, 10

class BaseSolver(object):
    defaults = dict(
        linear_solver = LinearSolver.default,
        preconditioner = Preconditioner.default,
    )

    def solve_constrained(self):
        raise NotImplementedError

    def solve_unconstrained(self):
        raise NotImplementedError

    def forms(self):
        raise NotImplementedError

    def __init__(self, problem, **params):
        self.problem = problem
        self.params = {**self.__class__.defaults,  **params}
        self.__dict__.update(self.params)

        uf, uo, pf, po = 'CG', 2, 'CG', 1
        self.ufe = VectorElement(uf, self.problem.mesh.ufl_cell(), uo)
        self.pfe = FiniteElement(pf, self.problem.mesh.ufl_cell(), po)
        self.rfe = FiniteElement('Real', self.problem.mesh.ufl_cell(), 0)

        self.solver_parameters = { 'linear_solver': self.linear_solver.name, 'preconditioner': self.preconditioner.name}
        self.path = problem.path

    def solve(self):
        if self.problem.constrained():
            self.solve_constrained()
        else:
            self.solve_unconstrained()

    def rename(self):
        self.u_.rename('u', 'Velocity')
        self.p_.rename('p', 'Pressure')
        self.eu.rename('eu', 'Error of velocity')
        self.ep.rename('ep', 'Error of pressure')
        self.du.rename('du', 'Difference of velocity time layers')
        self.dp.rename('dp', 'Difference of pressure time layers')

    def simulate(self, step):
        problem = self.problem

        fu  = XDMFFile(f'{self.path}/u.xdmf')
        fp  = XDMFFile(f'{self.path}/p.xdmf')
        feu = XDMFFile(f'{self.path}/eu.xdmf')
        fep = XDMFFile(f'{self.path}/ep.xdmf')

        self.rename()
        t = float(problem.t)
        fu.write(self.u_, t)
        fp.write(self.p_, t)
        feu.write(self.eu, t)
        fep.write(self.ep, t)

        timer = Timer()
        compute_time = 0.0
        while problem.hasnext():
            t, dt = problem.next()
            timer.start()
            step()
            timer.stop()
            compute_time += timer.elapsed()[0]

            fu.write(self.u_, t)
            fp.write(self.p_, t)
            feu.write(self.eu, t)
            fep.write(self.ep, t)

            if self.monitor is not None and self.monitor(self):
                break
        print(f'# compute time: {compute_time}')

        plt.figure('u');  plt.colorbar(plot(self.u_)); plt.savefig(f'{self.path}/u.png',  bbox_inches='tight'); plt.clf()
        plt.figure('p');  plt.colorbar(plot(self.p_)); plt.savefig(f'{self.path}/p.png',  bbox_inches='tight'); plt.clf()
        plt.figure('eu'); plt.colorbar(plot(self.eu)); plt.savefig(f'{self.path}/eu.png', bbox_inches='tight'); plt.clf()
        plt.figure('ep'); plt.colorbar(plot(self.ep)); plt.savefig(f'{self.path}/ep.png', bbox_inches='tight'); plt.clf()
