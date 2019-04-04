__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class ChannelFlow(BaseProblem):
    def __init__(self, **params):
        p = dict(
            f  = Constant((0, 0)),
        )
        super(ChannelFlow, self).__init__(**{**p, **params})
        mesh = self.mesh
        self.ud = [
            (Constant((0, 0)), mesh.walls),
            (Expression(mesh.profile, degree=2), mesh.inflow),
        ]
        self.ud1 = [(Constant(0), mesh.symm)]
