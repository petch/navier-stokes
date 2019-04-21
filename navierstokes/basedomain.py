from dolfin import *
from enum import Enum
import matplotlib.pyplot as plt
import subprocess
import pygmsh
import inspect
import os

class Dimension(Enum):
    two = 2
    three = 3

class BaseDomain(object):
    defaults = dict(
        dim       = Dimension.two,
        cellsize  = 0.05,
        cellscale = 1.0,
    )

    def geometry(self):
        raise NotImplementedError

    def __init__(self, title='Default', **params):
        self.title = title
        self.params = {**self.__class__.defaults, **params}
        self.__dict__.update(self.params)
        self.path = f'{os.path.dirname(inspect.getfile(self.__class__))}/results/{self.title}'
        os.makedirs(self.path, exist_ok=True)

    def execute(self, command):
        print(f"Running: {command}")
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            line = process.stdout.readline()
            if not line:
                break
            print(line.decode("utf-8"), end="")
        process.communicate()
        print(f"Exited with code {process.returncode}.")

    def generate(self):
        self.geom = pygmsh.opencascade.Geometry()
        self.shapes = []
        self.geometry()
        
        for shape in self.shapes:
            if hasattr(shape, 'surface'):
                self.geom.add_physical(shape.surface)
            if hasattr(shape, 'plane_surface'):
                self.geom.add_physical(shape.plane_surface)
            if hasattr(shape, 'lines'):
                for l in shape.lines:
                    self.geom.add_physical(l)

        open(f'{self.path}/domain.geo', 'w').write(self.geom.get_code())

        self.execute(f"gmsh -{self.dim.value} -clmin {self.cellsize} -clmax {self.cellsize} -clscale {self.cellscale} -o {self.path}/mesh.msh {self.path}/domain.geo")
        self.execute(f"dolfin-convert {self.path}/mesh.msh {self.path}/mesh.xml")
        os.rename(f"{self.path}/mesh_physical_region.xml", f"{self.path}/subdomains.xml")
        os.rename(f"{self.path}/mesh_facet_region.xml", f"{self.path}/boundaries.xml")

        self.mesh       = Mesh(f'{self.path}/mesh.xml')
        self.subdomains = MeshFunction('size_t', self.mesh, f'{self.path}/subdomains.xml')
        self.boundaries = MeshFunction('size_t', self.mesh, f'{self.path}/boundaries.xml')

        XDMFFile(f'{self.path}/mesh.xdmf').write(self.mesh)
        XDMFFile(f'{self.path}/subdomains.xdmf').write(self.subdomains)
        XDMFFile(f'{self.path}/boundaries.xdmf').write(self.boundaries)

        self.plot('mesh', self.mesh)
        self.plot('subdomains', self.subdomains)

    def plot(self, name, target):
        plt.figure(name)
        plot(target)
        plt.savefig(f'{self.path}/{name}.png', bbox_inches='tight')
        plt.clf()
