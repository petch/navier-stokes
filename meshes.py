from fenics import *
from mshr import *

def Square(N=16, diagonal='crossed'):
    mesh = UnitSquareMesh(N, N, diagonal)
    mesh.name = f'Square{N}'
    return mesh

def Channel(N=32, H=1.0, L=1.0):
    domain = Rectangle(Point(0, 0), Point(L, H))
    mesh = generate_mesh(domain, N)
    mesh.name     = f'Channel{N}'
    mesh.inflow   = f'near(x[0], 0)'
    mesh.outflow  = f'near(x[0], {L})'
    mesh.walls    = f'near(x[1], 0) || near(x[1], {H})'
    mesh.obstacle = None
    mesh.profile  = (f'4.0*x[1]*({H} - x[1])/{H}/{H}', '0')
    return mesh

def ChannelCylinder(N=32, H=0.41, L=2.2, R=0.05):
    domain = Rectangle(Point(0, 0), Point(L, H))
    domain -= Circle(Point(H/2, H/2), R)
    mesh = generate_mesh(domain, N)
    mesh.name     = f'ChannelCylinder{N}'
    mesh.inflow   = f'near(x[0], 0)'
    mesh.outflow  = f'near(x[0], {L})'
    mesh.walls    = f'near(x[1], 0) || near(x[1], {H})'
    mesh.obstacle = f'on_boundary && x[0]>{H/2-R*1.1} && x[0]<{H/2+R*1.1} && x[1]>{H/2-R*1.1} && x[1]<{H/2+R*1.1}'
    mesh.profile  = (f'4.0*x[1]*({H} - x[1])/{H}/{H}', '0')
    return mesh

def ChannelThin(N, H=0.4, L=1.0, h=0.1, l=0.3):
    domain = Polygon([Point(0.0, H-h), Point(l, H-h), Point(l+H-h, 0.0), Point(L, 0.0), Point(L, H), Point(0.0, H)])
    mesh = generate_mesh(domain, N)
    mesh.name     = f'ChannelThin{N}'
    mesh.inflow   = f'near(x[0], 0)'
    mesh.outflow  = f'near(x[0], {L})'
    mesh.walls    = f'on_boundary && (!near(x[0], 0.0) || near(x[1], {H-h}) || near(x[1], {H})) && (!near(x[0], {L}) || near(x[1], 0.0) || near(x[1], {H}))'
    mesh.obstacle = None
    mesh.profile  = (f'4.0*(x[1] - {H-h})*({H} - x[1])/{h}/{h}', '0')
    return mesh

def ChannelSymm(N, H=0.4, L=1.0, h=0.2, l=0.4):
    domain = Polygon([Point(0, 0), Point(L, 0), Point(L, H), Point(l + H - h, H), Point(l, h), Point(0, h)])
    mesh = generate_mesh(domain, N)
    mesh.name     = f'ChannelSymm{N}'
    mesh.inflow   = f'near(x[0], 0)'
    mesh.outflow  = f'near(x[0], {L})'
    mesh.walls    = f'on_boundary && !near(x[1], 0) && (!near(x[0], 0) || near(x[1], {h})) && (!near(x[0], {L}) || near(x[1], {H}))'
    mesh.symm     = f'near(x[1], 0)'
    mesh.profile  = (f'(x[1] + {h})*({h} - x[1])/{h}/{h}', '0')
    mesh.bounds = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    CompiledSubDomain(mesh.inflow).mark(mesh.bounds, 1)
    CompiledSubDomain(mesh.outflow).mark(mesh.bounds, 2)
    CompiledSubDomain(mesh.walls).mark(mesh.bounds, 3)
    CompiledSubDomain(mesh.symm).mark(mesh.bounds, 4)
    return mesh