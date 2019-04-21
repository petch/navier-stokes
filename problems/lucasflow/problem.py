from navierstokes import *

class Problem(BaseProblem):
    defaults = {**BaseProblem.defaults, **dict(
        ue = ('-cos(pi*x[0])/pi', '-x[1]*sin(pi*x[0])'),
        pe = '0',
        ce = '0',
    )}

    def conditions(self):
        self.ud = [
            ([2, 3, 4, 5], self.ue),
        ]

if __name__ == '__main__':
    Nonlinear(Problem()).solve()