__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class KanFlow(BaseProblem):
    def __init__(self, **params):
        p = dict(
            f  = Constant((0, 0)),
            ud = [
                (Expression(('0', '-sin(pi*(x[0]*x[0]*x[0]-3*x[0]*x[0]+3*x[0]))'), degree=2), 'near(x[1], 1)'),
                (Constant((0, 0)), 'near(x[0], 0) || near(x[1], 0)')
            ]
        )
        super(KanFlow, self).__init__(**{**p, **params})
