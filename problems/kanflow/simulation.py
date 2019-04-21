from navierstokes import *
from .domain import *
from .problem import *

class Simulation(BaseSimulation):
    defaults = {**BaseSimulation.defaults, **dict(
        solver = Nonlinear,
        linear_solver = LinearSolver.mumps,
    )}

    def prepare(self):
        self.Domain = Domain
        self.Problem = Problem
        self.runs = [
            dict(cellscale = [2.0, 1.0, 0.5, 0.25]),
        ]

if __name__ == '__main__':
    Simulation().run()
