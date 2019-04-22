from .domain import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        stationary = False,
        ue = ('(1+mu*t)*2*x[0]*x[0]*(x[0]-1)*(x[0]-1)*x[1]*(2*x[1]-1)*(x[1]-1)', '-(1+mu*t)*2*x[0]*(2*x[0]-1)*(x[0]-1)*x[1]*x[1]*(x[1]-1)*(x[1]-1)'),
        pe = 'rho*x[1]',
        ce = 'rho*0.5',
    )}

    def conditions(self):
        self.ud = [
            ([2, 3, 4, 5], self.ue),
        ]

if __name__ == '__main__':
    Domain.generate()
    Nonlinear(Problem()).solve()