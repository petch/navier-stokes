from .basesolver import *
from .logger import *
from pprint import pformat

class BaseSimulation(object):
    defaults = dict(
        level = LogLevel.WARNING,
        dudtmin = 1e-9,
        dudtmax = 1e+5,
    )

    def __init__(self):
        self.params = self.__class__.defaults
        self.params['monitor'] = self.monitor
        self.__dict__.update(self.params)

    def prepare(self):
        raise NotImplementedError
        
    def run(self):
        self.prepare()
        for vary in self.runs:
            self.vary_params(self.params, vary)

    def vary_params(self, params, vary):
        if len(vary) == 0:
            self.solve(params)
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
                t = f'{i}'
            p['title'] = f'{title}{key.title()}{t}'
            p[key] = v
            self.vary_params(p, vary)

    def solve(self, params):
        domain = self.Domain(**params)
        sys.stdout = Logger(f'{domain.path}/!out.txt', self.level)
        self.norms = open(f'{domain.path}/norms.txt', 'w')
        self.norms.write('t\tdudt\tdpdt\teumax\tepmax\teul2\tepl2\teuhd0\teph10\n')
        open(f'{domain.path}/params.txt', 'w').write(pformat(params))
        domain.generate()
        problem = self.Problem(**params)
        params['solver'](problem, **params).solve()
        sys.stdout = sys.__stdout__

    def monitor(self, solver):
        t  = float(solver.problem.t)
        dt = float(solver.problem.dt)
        dudt  = norm(solver.du)/dt
        dpdt  = norm(solver.dp)/dt
        eumax = norm(solver.eu.vector(), 'linf')
        epmax = norm(solver.ep.vector(), 'linf')
        eul2  = norm(solver.eu)
        epl2  = norm(solver.ep)
        euhd0 = norm(solver.eu, 'hdiv0')
        eph10 = norm(solver.ep, 'h10')
        norms_text = '{:.3f}\t{:.10f}\t{:.10f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}'.format(t, dudt, dpdt, eumax, epmax, eul2, epl2, euhd0, eph10)
        self.norms.write(norms_text + '\n')
        print(norms_text)
        return solver.problem.stationary and (dudt < self.dudtmin or dudt > self.dudtmax)
