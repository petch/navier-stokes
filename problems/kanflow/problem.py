from navierstokes import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        wall_u = (0.0, 0.0),
        inflow_u = ('0', '-sin(pi*(x[0]*x[0]*x[0]-3*x[0]*x[0]+3*x[0]))'),
        outflow_g = (0.0, 0.0),
    )}

    def conditions(self):
        self.ud = [
            ([4], self.inflow_u),
            ([2, 5], self.wall_u),
        ]
        self.un = [
            ([3], self.outflow_g),
        ]

if __name__ == '__main__':
    Nonlinear(Problem()).solve()