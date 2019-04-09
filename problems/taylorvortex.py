__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class TaylorVortex(BaseProblem):
    def __init__(self, **params):
        p = dict(
            ue = ('-cos(pi*x[0])*sin(pi*x[1])*exp(-2*pi*pi*mu*t)', 'sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*mu*t)'),
            pe = '-0.25*rho*(cos(2*pi*x[0])+cos(2*pi*x[1]))*exp(-4*pi*pi*mu*t)',
            ce = '0',
        )
        super(TaylorVortex, self).__init__(**{**p, **params})
