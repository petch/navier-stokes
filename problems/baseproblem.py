__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from dolfin import *

class BaseProblem(object):
    def __init__(self, **params):
        self.convection = getattr(self, params.get('convection', 'divergence'))

        self.t  = params.get('t', 0.0)
        self.T  = params.get('T', 10.0)
        self.dt = params.get('dt', 1.0)

        self.mesh = params.get('mesh', UnitSquareMesh(16, 16, 'crossed'))
        self.n    = FacetNormal(self.mesh)
        if hasattr(self.mesh, 'bounds'):
            self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.mesh.bounds)
        else:
            self.ds = ds

        self.tau = Constant(self.dt)
        self.rho = params.get('rho', Constant(1))
        self.mu  = params.get('mu',  Constant(1))

        self.ue = Expression(params.get('ue', ('0', '0')), degree=2, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        self.pe = Expression(params.get('pe', '0'), degree=1, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        self.ce = Expression(params.get('ce', '0'), degree=1, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        
        if self.ue is None or self.pe is None:
            self.f = params.get('f', Constant((0, 0)))
        else:
            self.f = self.source(self.ue, self.pe)

        self.ud = params.get('ud', [(self.ue, 'on_boundary')])
        self.pd = params.get('pd', [])

    def divergence(self, u1, u2):
        return div(outer(u1, u2))

    def advective(self, u1, u2):
        return dot(u1, nabla_grad(u2))

    def skewsymmetric(self, u1, u2):
        return 0.5*self.divergence(u1, u2) + 0.5*self.advective(u1, u2)

    def rotation(self, u1, u2):
        return curl(u2)*as_vector((-u1[1], u1[0])) + 0.5*nabla_grad(dot(u1, u2))

    def stress(self, u, p):
        return 2*self.mu*sym(nabla_grad(u)) - p*Identity(len(u))

    def source(self, u, p):
        return self.rho*self.convection(u, u) - div(self.stress(u, p))
