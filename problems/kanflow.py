__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class KanFlow(BaseProblem):
    def __init__(self, **params):
        super(KanFlow, self).__init__(**params)
        self.ud = [
            (Bound.TOP, Expression(('0', '-sin(pi*(x[0]*x[0]*x[0]-3*x[0]*x[0]+3*x[0]))'), degree=2)),
            (Bound.LEFT, Constant((0, 0))),
            (Bound.BOTTOM, Constant((0, 0))),
        ]
        self.un = [
            (Bound.RIGHT, Constant((0, 0)), None)
        ]
