from .domain import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        stationary = False,
        mu = 0.01,
        ue = ('-cos(pi*x[0])*sin(pi*x[1])*exp(-2*pi*pi*mu*t)', 'sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*mu*t)'),
        pe = '-0.25*rho*(cos(2*pi*x[0])+cos(2*pi*x[1]))*exp(-4*pi*pi*mu*t)',
        ce = '0',
    )}

    def conditions(self):
        self.ud = [
            ([2, 3, 4, 5], self.ue),
        ]

if __name__ == '__main__':
    Domain().generate()
    Nonlinear(Problem()).solve()