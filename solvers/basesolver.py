__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from dolfin import *
from mpi4py import MPI
import sys
import os
import matplotlib.pyplot as plt

class Logger(object):
    def __init__(self, path, level):
        set_log_level(level)
        self.rank = MPI.COMM_WORLD.Get_rank()
        if (self.rank != 0):
            return
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.out = sys.stdout
        self.log = open(path, "w")
    def write(self, message):
        if (self.rank != 0):
            return
        self.out.write(message)
        self.log.write(message)
    def flush(self):
        if (self.rank != 0):
            return
        self.out.flush()
        self.log.flush()

class BaseSolver(object):
    def __init__(self, problem_data, **params):
        self.problem = problem_data

        self.stationary = params.get('stationary', True)
        self.dudtmin    = params.get('dudtmin', 1e-9)
        self.dudtmax    = params.get('dudtmax', 1e+5)

        uf, uo, pf, po = params.get('elements', ('CG', 2, 'CG', 1))
        self.ufe = VectorElement(uf, self.problem.mesh.ufl_cell(), uo)
        self.pfe = FiniteElement(pf, self.problem.mesh.ufl_cell(), po)
        self.rfe = FiniteElement('Real', self.problem.mesh.ufl_cell(), 0)

        self.linear_solver     = params.get('linear_solver', 'default')
        self.preconditioner    = params.get('preconditioner', 'default')
        self.solver_parameters = { 'linear_solver': self.linear_solver, 'preconditioner': self.preconditioner}

        self.level = params.get('level', LogLevel.WARNING)
        self.title = params.get('title', 'Default')
        self.path  = params.get('path', f'results/{self.problem.__class__.__name__}/{self.title}/')

    def solve(self):
        if self.problem.constrained():
            self.solve_constrained()
        else:
            self.solve_unconstrained()

    def solve_constrained(self):
        raise NotImplementedError

    def solve_unconstrained(self):
        raise NotImplementedError

    def forms(self):
        raise NotImplementedError

    def simulate(self, step):
        sys.stdout = Logger(f'{self.path}!out.txt', self.level)
        print(f'# Solving {self.title}')
        problem = self.problem

        fu  = XDMFFile(f'{self.path}u.xdmf')
        fp  = XDMFFile(f'{self.path}p.xdmf')
        feu = XDMFFile(f'{self.path}eu.xdmf')
        fep = XDMFFile(f'{self.path}ep.xdmf')

        timer = Timer()
        compute_time = 0.0
        print('time\tdudt\tdpdt\teumax\tepmax\teul2\tepl2\teuhd0\teph10')
        while problem.t < problem.T - problem.dt/2:
            problem.t += problem.dt

            timer.start()
            step()
            timer.stop()
            compute_time += timer.elapsed()[0]

            fu.write(self.u_, problem.t)
            fp.write(self.p_, problem.t)
            feu.write(self.eu, problem.t)
            fep.write(self.ep, problem.t)

            dudt  = norm(self.du)/problem.dt
            dpdt  = norm(self.dp)/problem.dt
            eumax = norm(self.eu.vector(), 'linf')
            epmax = norm(self.ep.vector(), 'linf')
            eul2  = norm(self.eu)
            epl2  = norm(self.ep)
            euhd0 = norm(self.eu, 'hdiv0') if self.ufe.degree() > 0 else 0
            eph10 = norm(self.ep, 'h10') if self.pfe.degree() > 0 else 0
            print('{:.3f}\t{:.10f}\t{:.10f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}'.format(problem.t, dudt, dpdt, eumax, epmax, eul2, epl2, euhd0, eph10))
            if self.stationary and (dudt < self.dudtmin or dudt > self.dudtmax):
                break
        print(f'# compute time: {compute_time}')

        plt.figure('mesh'); plot(problem.mesh); plt.savefig(f'{self.path}mesh.png', bbox_inches='tight'); plt.clf()
        plt.figure('u');  plt.colorbar(plot(self.u_)); plt.savefig(f'{self.path}u.png',  bbox_inches='tight'); plt.clf()
        plt.figure('p');  plt.colorbar(plot(self.p_)); plt.savefig(f'{self.path}p.png',  bbox_inches='tight'); plt.clf()
        plt.figure('eu'); plt.colorbar(plot(self.eu)); plt.savefig(f'{self.path}eu.png', bbox_inches='tight'); plt.clf()
        plt.figure('ep'); plt.colorbar(plot(self.ep)); plt.savefig(f'{self.path}ep.png', bbox_inches='tight'); plt.clf()

        sys.stdout = sys.__stdout__
