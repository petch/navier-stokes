__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from problems.baseproblem import *

class ShiftedVortex(BaseProblem):
    def __init__(self, **params):
        p = dict(
            mu = Constant(0.01),
            dt = 0.01,
            T  = 100.0,
            ue = ('40*2*exp(x[0])*pow(x[0]-1,2)*pow(x[0],2)*x[1]*(x[1]-1)*(2*x[1]-1)', '-40*exp(x[0])*(x[0]-1)*x[0]*(x[0]*(3+x[0])-2)*pow(x[1]-1,2)*pow(x[1],2)'),
            pe = '10*(-424+156*exp(1)+(x[1]*x[1]-x[1])*(-456+exp(x[0])*(456+pow(x[0],2)*(228-5*(x[1]*x[1]-x[1]))+2*x[0]*(-228+(x[1]*x[1]-x[1]))+2*pow(x[0],3)*(-36+(x[1]*x[1]-x[1]))+pow(x[0],4)*(12+(x[1]*x[1]-x[1])))))',
        )
        super(ShiftedVortex, self).__init__(**{**p, **params})
