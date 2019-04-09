__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *
from meshes import Bound

class ChannelFlow(BaseProblem):
    def __init__(self, **params):
        super(ChannelFlow, self).__init__(**params)
        self.ud = [
            (Bound.TOP, Constant((0, 0))),
            (Bound.OTHER, Constant((0, 0))),
            (Bound.LEFT, Expression(('(0.2-x[1])*(0.2+x[1])/0.2/0.2', '0'), degree=2)),
            (Bound.BOTTOM, Constant(0), 1),
        ]
        self.un = [
            (Bound.RIGHT, Constant((0, 0)), None)
        ]
        self.pd = [
            (Bound.RIGHT, Constant(0))
        ]