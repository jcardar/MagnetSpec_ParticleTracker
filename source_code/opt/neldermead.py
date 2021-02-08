import numpy
import minitest

from . import OptimizationProblem, Goal

class NelderMead:
    """An interface to the Nelder-Mead optimization algorithm."""
    def __init__(self, opt_problem):
        self.problem = opt_problem
        self.simplex = None
        self.values = None
        self.params = {
                'reflect': 1,
                'expand': 2,
                'contract': 0.5,
                'shrink': 0.5,
                }
    def init_simplex(self, selection):
        """Initializes the simplex from a selection of sample points.

        The argument must be a sequence of at least n + 1 points, where n is the
        dimensionality of the parameter space.
        """
        n = numpy.prod(self.problem.shape)
        if len(selection) < n + 1:
            raise ValueError('not enough samples to initialize simplex')
        self.simplex = numpy.array(selection)[:n+1]
        self.values = []
    def seed(self, centre, scale, seed=None):
        """Initializes the simplex from a seed.

        Initializes the simplex from an rng seed, a centre point, and a scale.
        """
        n = numpy.prod(self.problem.shape)
        rng = numpy.random.RandomState(seed)
        centre = numpy.asarray(centre)
        scale = numpy.asarray(scale)

        selection = []
        if self.problem.validate(centre):
            selection.append(centre)

        for idx in numpy.ndindex(*self.problem.shape):
            val = numpy.array(centre)
            if rng.randint(2):
                val[idx] += scale[idx]
            else:
                val[idx] -= scale[idx]
            if self.problem.validate(val):
                selection.append(val)

        while len(selection) < n + 1:
            diff = rng.normal(0, 1, self.problem.shape)
            val = centre + scale * diff
            if self.problem.validate(val):
                selection.append(val)

        self.init_simplex(selection)

    def init_values(self):
        if self.simplex is None:
            raise RuntimeError('simplex not initialized')
        if type(self.values) is not list:
            raise RuntimeError('values already initialized')
        # Don't evaluate samples already evaluated
        for i in range(len(self.values), len(self.simplex)):
            val = self.problem.evaluate_min(self.simplex[i])
            if val is None:
                raise ValueError('invalid point in initial simplex')
            self.values.append(val)
        self.values = numpy.array(self.values)
        self.best_idx = numpy.argmin(self.values)

    def step(self):
        """Runs a single step of the Nelder-Mead algorithm."""
        if self.simplex is None:
            raise RuntimeError('simplex not initialized')
        if type(self.values) is list:
            self.init_values()
        order = numpy.argsort(self.values)
        # order[0] is the index of the best (minimum) point
        centroid = numpy.mean(self.simplex[order[:-1]], axis=0)
        worst = self.simplex[order[-1]]

        # Reflect
        for z in (1, 0.5):
            reflect = centroid + z * self.params['reflect'] * (centroid - worst)
            reflect_val = self.problem.evaluate_min(reflect)
            if reflect_val is None:
                # Try again with less reflection
                continue
            if self.values[order[0]] <= reflect_val < self.values[order[-2]]:
                self.simplex[order[-1]] = reflect
                self.values[order[-1]] = reflect_val
                return
            break
        
        # Expand
        if reflect_val is not None:
            if reflect_val < self.values[order[0]]:
                expand = centroid + self.params['expand'] * (reflect - centroid)
                expand_val = self.problem.evaluate_min(expand)
                if expand_val is not None and expand_val < reflect_val:
                    self.simplex[order[-1]] = expand
                    self.values[order[-1]] = expand_val
                else:
                    self.simplex[order[-1]] = reflect
                    self.values[order[-1]] = reflect_val
                self.best_idx = order[-1]
                return

        contract = centroid + self.params['contract'] * (worst - centroid)
        contract_val = self.problem.evaluate_min(contract)
        if contract_val is not None and contract_val < self.values[order[-1]]:
            self.simplex[order[-1]] = contract
            self.values[order[-1]] = contract_val
            if contract_val < self.values[order[-1]]:
                self.best_idx = order[-1]
        else:
            self.shrink(order)

    def shrink(self, order):
        # This is the shrink step of the Nelder-Mead algorithm
        # Shrink every point but the best
        best = self.simplex[order[0]]
        for i in order[1:]:
            point = best + self.params['shrink'] * (self.simplex[i] - best)
            value = self.problem.evaluate_min(point)
            if value is not None:
                self.simplex[i] = point
                self.values[i] = value
        self.best_idx = numpy.argmin(self.values)

    def optimize(self, steps):
        """Runs the optimization algorithm for a given number of steps."""
        for i in range(steps):
            self.step()
        return self.simplex[self.best_idx]
    
    def best_point(self):
        """The best point found so far."""
        if self.simplex is None:
            raise RuntimeError('simplex not initialized')
        if self.values is None:
            self.init_values()
        return numpy.array(self.simplex[self.best_idx])

    def best_val(self):
        """The optimal value found so far."""
        if self.simplex is None:
            raise RuntimeError('simplex not initialized')
        if self.values is None:
            self.init_values()
        val = self.values[self.best_idx]
        return val if self.problem.goal == Goal.minimize else -val

def rosenbrock(x, y):
    a = 1
    b = 100
    return (a - x)**2 + b * (y - x**2)**2

@minitest.test
def test_banana():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.minimize, lbound=(-2, -1), ubound=(2, 3))
    nm = NelderMead(problem)
    nm.seed((0.0, 2.0), (0.1, 0.2), 45)
    r = nm.optimize(100)
    assert numpy.all(r == nm.best_point())
    assert problem.evaluate(r) == nm.best_val()
    assert problem.evaluate(r) < 1

@minitest.test
def test_multidim():
    problem = OptimizationProblem((2, 3, 4), lambda x: numpy.sum(x**2),
            Goal.minimize)
    nm = NelderMead(problem)
    nm.seed(numpy.ones((2, 3, 4)), numpy.ones((2, 3, 4)), 100)
    r = nm.optimize(200)
    assert problem.evaluate(r) < 10


@minitest.test
def test_exception_safe():
    eval_counter = 0
    def evaluate(params):
        nonlocal eval_counter
        eval_counter += 1
        if eval_counter == 25:
            raise RuntimeError
        return rosenbrock(*params)

    problem = OptimizationProblem((2,), evaluate,
            Goal.minimize, lbound=(-2, -1), ubound=(2, 3))

    nm = NelderMead(problem)
    nm.seed((0.0, 2.0), (0.1, 0.2), 45)
    try:
        nm.optimize(50)
    except RuntimeError:
        pass
    else:
        assert False

    r = nm.optimize(50)

    assert problem.evaluate(r) < 10

@minitest.test
def test_exception_safe_2():
    eval_counter = 0
    def evaluate(params):
        nonlocal eval_counter
        eval_counter += 1
        if eval_counter == 3:
            raise RuntimeError
        return rosenbrock(*params)

    problem = OptimizationProblem((2,), evaluate,
            Goal.minimize, lbound=(-2, -1), ubound=(2, 3))

    nm = NelderMead(problem)
    nm.seed((0.0, 2.0), (0.1, 0.2), 45)
    try:
        nm.step()
    except RuntimeError:
        pass
    else:
        assert False

    nm.step()

    # This checks that the first two points in the simplex were not evaluated twice
    assert eval_counter <= 6

    r = nm.optimize(100)

    assert problem.evaluate(r) < 10

@minitest.runner
def test(interactive):
    pass
