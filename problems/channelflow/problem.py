from navierstokes import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        h         = 0.2,
        wall_u    = (0.0, 0.0),
        bottom_u  = 0.0,
        outflow_g = (0.0, 0.0),
        outflow_p = 0.0,
        inflow_u  = ('(h-x[1])*(h+x[1])/h/h', '0'),
    )}

    def conditions(self):
        self.ud = [
            ([4, 5, 6], self.wall_u),
            ([7], self.inflow_u),
            ([2], self.bottom_u, 1),
        ]
        self.un = [
            ([3], self.outflow_g)
        ]
        self.pd = [
            ([3], self.outflow_p)
        ]

if __name__ == '__main__':
    Nonlinear(Problem()).solve()