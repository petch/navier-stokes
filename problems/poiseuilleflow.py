__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class PoiseuilleFlow(BaseProblem):
    def __init__(self, **params):
        p = dict(
            ue = ('0.5/mu*x[1]*(1-x[1])', '0'),
            pe = '1.0-x[0]',
            pa = '0.5',
            ud = [(Constant((0, 0)), 'near(x[1], 0) || near(x[1], 1)')],
            pd = [(Constant(1), 'near(x[0], 0)'), (Constant(0), 'near(x[0], 1)')]
        )
        super(PoiseuilleFlow, self).__init__(**{**p, **params})
