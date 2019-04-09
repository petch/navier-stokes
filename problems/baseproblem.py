__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

from dolfin import *
from meshes import *

class BaseProblem(object):
    def __init__(self, **params):
        self.convection = getattr(self, params.get('convection', 'divergence'))

        self.t  = params.get('t', 0.0)
        self.T  = params.get('T', 10.0)
        self.dt = params.get('dt', 1.0)

        self.mesh    = params.get('mesh', Square(16))
        self.domains = self.mesh.domains
        self.bounds  = self.mesh.bounds

        self.dx = Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.bounds)
        self.n  = FacetNormal(self.mesh)

        self.tau = Constant(self.dt)
        self.rho = params.get('rho', Constant(1))
        self.mu  = params.get('mu',  Constant(1))

        self.ue = Expression(params.get('ue', ('0', '0')), degree=2, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        self.pe = Expression(params.get('pe', '0'), degree=1, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        self.ce = Expression(params.get('ce', '0'), degree=1, domain=self.mesh, rho=self.rho, mu=self.mu, t=self.t)
        
        self.f = params.get('f', self.source(self.ue, self.pe))

        self.ud = params.get('ud', [])
        self.pd = params.get('pd', [])
        self.un = params.get('un', [])

    def dirichlet_u(self, V):
        bcs = []
        for bc in self.ud:
            space = V if len(bc) < 3 else V.sub(bc[2])
            bcs.append(DirichletBC(space, bc[1], self.bounds, bc[0]))
        return bcs

    def dirichlet_p(self, Q):
        bcs = []
        for bc in self.pd:
            bcs.append(DirichletBC(Q, bc[1], self.bounds, bc[0]))
        return bcs

    def constrained(self):
        return len(self.pd) == 0

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

    def stress_n(self, u, p, n, dun=None, dut=None):
        sn = -dot(p*Identity(len(u)), n)
        if dun is None:
            sn += dot(self.mu*nabla_grad(u).T, n)
        else:
            sn += self.mu*dun
        if dut is None:
            sn += dot(self.mu*nabla_grad(u), n)
        else:
            sn += self.mu*dut
        return sn

    def source(self, u, p):
        return self.rho*self.convection(u, u) - div(self.stress(u, p))
