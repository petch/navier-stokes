from .domain import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        ue = ('0.5/mu*x[1]*(1-x[1])', '0'),
        pe = '1.0-x[0]',
        wall_u = (0.0, 0.0),
        inflow_p = 1.0,
        inflow_u = (0.0, 0.0),
        outflow_p = 0.0,
        outflow_u = (0.0, 0.0),
    )}

    def conditions(self):
        self.ud = [
            ([2, 4], self.wall_u),
        ]
        self.pd = [
            ([5], self.inflow_p),
            ([3], self.outflow_p),
        ]
        self.un = [
            ([5], self.inflow_u),
            ([3], self.outflow_u),
        ]

if __name__ == '__main__':
    Domain().generate()
    Nonlinear(Problem()).solve()