from fenics import *
from mshr import *
from enum import Enum

class Bound(Enum):
    LEFT   = 1
    RIGHT  = 2
    BOTTOM = 3
    TOP    = 4
    FRONT  = 5
    BACK   = 6
    OTHER  = 7
    def __int__(self):
        return self.value

def Square(N=16, diagonal='crossed'):
    mesh = UnitSquareMesh(N, N, diagonal)
    mesh.domains = MeshFunction('size_t', mesh, mesh.topology().dim())
    mesh.bounds  = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    CompiledSubDomain('near(x[0], 0)').mark(mesh.bounds, Bound.LEFT)
    CompiledSubDomain('near(x[0], 1)').mark(mesh.bounds, Bound.RIGHT)
    CompiledSubDomain('near(x[1], 0)').mark(mesh.bounds, Bound.BOTTOM)
    CompiledSubDomain('near(x[1], 1)').mark(mesh.bounds, Bound.TOP)
    return mesh

def ChannelCylinder(N=16, H=0.41, L=2.2, R=0.05):
    geo  = Rectangle(Point(0, 0), Point(L, H)) - Circle(Point(H/2, H/2), R)
    mesh = generate_mesh(geo, N)
    mesh.domains = MeshFunction('size_t', mesh, mesh.topology().dim())
    mesh.bounds  = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    CompiledSubDomain(f'on_boundary').mark(mesh.bounds, Bound.OTHER)
    CompiledSubDomain(f'near(x[0], {0})').mark(mesh.bounds, Bound.LEFT)
    CompiledSubDomain(f'near(x[0], {L})').mark(mesh.bounds, Bound.RIGHT)
    CompiledSubDomain(f'near(x[1], {0})').mark(mesh.bounds, Bound.BOTTOM)
    CompiledSubDomain(f'near(x[1], {H})').mark(mesh.bounds, Bound.TOP)
    return mesh

def ChannelThin(N=16, H=0.4, L=1.0, h=0.2, l=0.4):
    geo  = Polygon([Point(0, 0), Point(L, 0), Point(L, H), Point(l + H - h, H), Point(l, h), Point(0, h)])
    mesh = generate_mesh(geo, N)
    mesh.domains = MeshFunction('size_t', mesh, mesh.topology().dim())
    mesh.bounds  = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    CompiledSubDomain(f'on_boundary').mark(mesh.bounds, Bound.OTHER)
    CompiledSubDomain(f'near(x[0], {0})').mark(mesh.bounds, Bound.LEFT)
    CompiledSubDomain(f'near(x[0], {L})').mark(mesh.bounds, Bound.RIGHT)
    CompiledSubDomain(f'near(x[1], {0})').mark(mesh.bounds, Bound.BOTTOM)
    CompiledSubDomain(f'near(x[1], {H})').mark(mesh.bounds, Bound.TOP)
    return mesh