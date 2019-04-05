__author__ = "Petr Zakharov <zapetch@gmail.com>"
__date__ = "2019-03-29"
__copyright__ = "Copyright (C) 2019 " + __author__
__license__  = "CC0"

import inspect
from solvers.nonlinear import *
from solvers.linear import *
from solvers.ipcs import *
from solvers.duoglas import *

solvers_list = [Nonlinear, Linear, IPCS, Duoglas]

def run(params, **vary):
    if len(vary) == 0:
        problem = params['problem_class'](**params)
        solver  = params['solver_class'](problem, **params)
        solver.solve()
        return
    key, values = vary.popitem()
    title = params.get('title', '')
    p = params.copy()
    for i in range(len(values)):
        v = values[i]
        if isinstance(v, int) or isinstance(v, float):
            t = f'{v}'
        elif isinstance(v, str):
            t = f'{v.title()}'
        elif inspect.isclass(v):
            t = v.__name__
        else:
            t = f'{key.title()}{i}'
        p['title'] = f'{title}{t}'
        p[key] = v
        run(p, **vary)

