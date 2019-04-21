from dolfin import *
from mpi4py import MPI
import sys
import os

class Logger(object):
    def __init__(self, path, level):
        set_log_level(level)
        self.rank = MPI.COMM_WORLD.Get_rank()
        if (self.rank != 0):
            return
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.out = sys.stdout
        self.log = open(path, "w")
    def write(self, message):
        if (self.rank != 0):
            return
        self.out.write(message)
        self.log.write(message)
    def flush(self):
        if (self.rank != 0):
            return
        self.out.flush()
        self.log.flush()
