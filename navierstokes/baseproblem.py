from dolfin import *
from enum import Enum
import inspect
import os

class Convection(Enum):
    divergence    = 0
    advective     = 1
    skewsymmetric = 2
    rotation      = 3

class BaseProblem(object):
    defaults = dict(
        convection = Convection.divergence,
        stationary = True,
        t  = 0.0,
        T  = 10.0,
        dt  = 1.0,
        rho = 1.0,
        mu  = 1.0,
        ue = ('0', '0'),
        pe = '0',
        ce = '0',
        ud = [],
        pd = [],
        un = [],
    )

    def __init__(self, title='default', **params):
        self.title = title
        self.params = {**self.__class__.defaults,  **params}
        self.__dict__.update(self.params)

        self.convection = getattr(self, self.convection.name)

        self.path = f'{os.path.dirname(inspect.getfile(self.__class__))}/results/{self.title}'
        self.mesh = Mesh(f'{self.path}/mesh.xml')
        self.subdomains = MeshFunction('size_t', self.mesh, f'{self.path}/subdomains.xml') 
        self.boundaries = MeshFunction('size_t', self.mesh, f'{self.path}/boundaries.xml') 

        self.dx = Measure('dx', domain=self.mesh, subdomain_data=self.subdomains)
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = FacetNormal(self.mesh)

        self.constants = dict()
        for k, v in self.params.items():
            if isinstance(v, float) or isinstance(v, tuple) and isinstance(v[0], float):
                self.constants[k] = Constant(v)
        self.__dict__.update(self.constants)

        self.expressions = dict()
        for k, v in self.params.items():
            if isinstance(v, str) or isinstance(v, tuple) and isinstance(v[0], str):
                self.expressions[k] = Expression(v, degree=2, domain=self.mesh, **self.constants)
        self.__dict__.update(self.expressions)

        if 'f' not in self.params:
            self.f = self.source(self.ue, self.pe)

        self.conditions()

    def hasnext(self):
        return float(self.t) < float(self.T) - float(self.dt)/2

    def next(self):
        self.t.assign(float(self.t) + float(self.dt))
        return float(self.t), float(self.dt)

    def dirichlet_u(self, V):
        bcs = []
        for bc in self.ud:
            space = V if len(bc) < 3 else V.sub(bc[2])
            for i in bc[0]:
                bcs.append(DirichletBC(space, bc[1], self.boundaries, i))
        return bcs

    def dirichlet_p(self, Q):
        bcs = []
        for bc in self.pd:
            for i in bc[0]:
                bcs.append(DirichletBC(Q, bc[1], self.boundaries, i))
        return bcs

    def neumann_u(self, u, p, v):
        F = dot(Constant((0, 0)), v)*self.ds
        for indices, dun in self.un:
            for i in indices:
                F -= dot(self.stress_n(u, p, self.n, dun, None), v)*self.ds(i)
        return F

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
