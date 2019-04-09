__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class LucasFlow(BaseProblem):
    def __init__(self, **params):
        p = dict(
            ue = ('-cos(pi*x[0])/pi', '-x[1]*sin(pi*x[0])'),
            pe = '0',
            ce = '0',
        )
        super(LucasFlow, self).__init__(**{**p, **params})
        self.ud = [
            (Bound.LEFT,   self.ue),
            (Bound.RIGHT,  self.ue),
            (Bound.TOP,    self.ue),
            (Bound.BOTTOM, self.ue),
        ]
