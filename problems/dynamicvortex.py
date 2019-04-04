__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class DynamicVortex(BaseProblem):
    def __init__(self, **params):
        p = dict(
            ue = ('(1+mu*t)*2*x[0]*x[0]*(x[0]-1)*(x[0]-1)*x[1]*(2*x[1]-1)*(x[1]-1)', '-(1+mu*t)*2*x[0]*(2*x[0]-1)*(x[0]-1)*x[1]*x[1]*(x[1]-1)*(x[1]-1)'),
            pe = 'rho*x[1]',
            ce = 'rho*0.5',
        )
        super(DynamicVortex, self).__init__(**{**p, **params})
