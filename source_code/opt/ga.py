import random
import numpy
import minitest

from . import OptimizationProblem, Goal

class BasicGaDriver:
    """A straightforward but flexible genetic algorithm.

    This is a naive implementation of a genetic algorithm. In order to determine
    the size of any mutation, it relies on a mutation_size array, with the same
    shape as the data. Its behaviour can be further configured through a number
    of keys:

        population: the number of members of each generation
        selection: the number of members selected to form the next generation
        elitism: the number of members that enter the next generation unchanged

        mutation_size: a global multiplier for the mutation size
        mutation_count: the number of array elements that will be mutated
        mutation_probability: the probability that any given array element will
            be mutated. Independent of mutation_count

        crossover_rate: the fraction of each genome taken from the second parent

    Useful attributes:
        generation_number: the number of generations of the GA evaluated
        population: the members of the current generation
        selection: the members selected to form the next generation, from best
            to worst
        selection_fitness: the fitness values corresponding to the members of
            the selection list
    """
    def __init__(self, opt_problem, mutation_size):
        self.problem = opt_problem
        mutation_size = numpy.array(mutation_size)
        if mutation_size.shape != opt_problem.shape:
            raise ValueError('mutation_size has an incorrect shape - {}, '
                'should be {}'.format(mutation_size.shape, opt_problem.shape))
        self.mutation_size = mutation_size
        self._fitness = None
        self.fitness = [] ####
        self.generation_number = 0
        self.population = []
        self.selection = []
        self.selection_fitness = []
        self.random = random.Random()
        self.nprandom = numpy.random.RandomState()
        self.params = {
                'population': 10,
                'selection': 3,
                'elitism': 5,
                'mutation_size': 1.0,
                'mutation_count': 0,
                'mutation_probability': 1.0,
                'crossover_rate': 0.5,
                }
    def __getitem__(self, key):
        return self.params[key]
    def __setitem__(self, key, value):
        if key in self.params:
            self.params[key] = value
        else:
            raise KeyError('{} is not a valid parameter name'.format(key))
    def seed(self, seed):
        """Seeds the algorithm's internal random number generator."""
        self.random.seed(seed)
        self.nprandom.seed(seed)
    def uniform(self):
        """Initialises the first generation across all parameter space.

        Uses a uniform distribution across all parameter space to initialise the
        first generation. Only works if both lower and upper bounds are set.
        """
        lbound, ubound = self.problem.bounds()
        if lbound is None:
            raise RuntimeError('no lower bound set')
        if ubound is None:
            raise RuntimeError('no upper bound set')

        self._fitness = None
        self.generation_number = 0
        self.population.clear()
        self.selection.clear()
        self.selection_fitness.clear()
        self.fitness.clear() ###

        while len(self.population) < self['population']:
            value = self.nprandom.uniform(lbound, ubound)
            if self.problem.validate(value):
                self.population.append(value)
    def populate(self, point_zero, scale=2.0):
        """Initialises the first generation from a central point.

        Population members are generated by applying mutations to the supplied
        point. The various mutation-related parameters are ignored. The optional
        scale argument can change the size of these initial mutations.
        """
        point_zero = numpy.asarray(point_zero)
        if point_zero.shape != self.problem.shape:
            raise ValueError('point_zero has an incorrect shape - {}, '
                'should be {}'.format(point_zero.shape, self.problem.shape))
        self._fitness = None
        self.generation_number = 0
        self.population.clear()
        self.selection.clear()
        self.selection_fitness.clear()
        self.fitness = [] ####
        if self['elitism'] and self.problem.validate(point_zero):
            self.population.append(point_zero)
        while len(self.population) < self['population']:
            mutation = self.nprandom.normal(0.0, scale, self.problem.shape)
            mutation *= self.mutation_size
            value = point_zero + mutation
            value = self.problem.clip(value)
            if self.problem.validate(value):
                self.population.append(value)

    def crossover(self, ma, pa):
        """The crossover operator."""
        choices = self.nprandom.binomial(1, self['crossover_rate'], ma.shape)
        choices = choices == 1
        child = numpy.array(ma)
        child[choices] = pa[choices]
        return child

    def mutate(self, value):
        """The mutation operator."""
        if self['mutation_size'] <= 0:
            return self.problem.clip(value)
        mutations = self.nprandom.normal(0, self['mutation_size'], value.shape)
        mutations *= self.mutation_size
        mut_idxs = numpy.zeros(value.shape)
        if self['mutation_count']:
            indices = self.random.sample(range(len(mut_idxs.flat)),
                    k=self['mutation_count'])
            mut_idxs.flat[indices] = 1
        if self['mutation_probability']:
            mut_idxs += self.nprandom.binomial(1, self['mutation_probability'],
                    mut_idxs.shape)
        mut_idxs = mut_idxs > 0
        value[mut_idxs] += mutations[mut_idxs]
        return self.problem.clip(value)

    def generate(self):
        """Executes one generation of the genetic algorithm."""
        if self._fitness is None:
            self._begin_generation()
        self._resume_generation()

    def _begin_generation(self):
        # If an exception occurs during evaluation, the generate function may
        # be called to resume evaluation of the current generation. This is
        # indicated by self._fitness not being None. This function is only
        # executed at the start of a generation
        if self.generation_number > 0:
            # Produce the next generation by crossover and mutation
            self.population.clear()
            for i in range(self['elitism']):
                self.population.append(self.selection[i])
            while len(self.population) < self['population']:
                # Choose two parents, or possibly the same parent twice
                ma = self.random.choice(self.selection)
                pa = self.random.choice(self.selection)
                child = self.crossover(ma, pa)
                child = self.mutate(child)
                if self.problem.validate(child):
                    self.population.append(child)
        if not self.population:
            raise RuntimeError('call populate() before generate()')
        self.generation_number += 1
        self.selection.clear()
        self.selection_fitness.clear()
        self._fitness = []
        self.fitness = []

    def _resume_generation(self):
        # This function is called to resume and complete the evaluation of a
        # generation. The number of population members already evaluated is
        # indicated by the length of self._fitness
        start_from = len(self._fitness)
        for x in self.population[start_from:]:
            self._fitness.append(self.problem.evaluate(x))

        if self.problem.goal == Goal.maximize:
            rev = True
        elif self.problem.goal == Goal.minimize:
            rev = False
        else:
            raise ValueError('unknown Goal: {}'.format(self.problem.goal))

        fitness = sorted(enumerate(self._fitness), key=lambda x: x[1], reverse=rev)
        selection = fitness[:self.params['selection']]
        self.selection_fitness = []
        self.selection = []
        for (i, f) in selection:
            self.selection_fitness.append(f)
            self.selection.append(self.population[i])
        self._fitness = None


        ####
        for (i, f) in fitness:
            self.fitness.append(f)

class DifferentialGaDriver(BasicGaDriver):
    """A genetic algorithm that uses differential mutation.

    In this algorithm, the mutations are derived by taking the differences
    between randomly-selected members of the population. This ensures that the
    mutations become smaller as the algosrithm approaches convergence.

    This algorithm is inspired by the differential evolution algorithm that
    appears in scipy.
    """
    def __init__(self, opt_problem):
        super().__init__(opt_problem, numpy.zeros(opt_problem.shape))
        del self.mutation_size

    def mutate(self, value):
        """The mutation operator."""
        if self['mutation_size'] <= 0:
            return self.problem.clip(value)
        a, b = self.random.sample(self.selection, k=2)
        mutations = (a - b) * self.nprandom.normal(0, self['mutation_size'])
        mut_idxs = numpy.zeros(value.shape)
        if self['mutation_count']:
            indices = self.random.sample(range(len(mut_idxs.flat)),
                    k=self['mutation_count'])
            mut_idxs.flat[indices] = 1
        if self['mutation_probability']:
            mut_idxs += self.nprandom.binomial(1, self['mutation_probability'],
                    mut_idxs.shape)
        mut_idxs = mut_idxs > 0
        value[mut_idxs] += mutations[mut_idxs]
        return self.problem.clip(value)

    def populate(self, point_zero, scale):
        """Initialises the first generation from a central point.

        Population members are generated by applying mutations to the supplied
        point. The sizes of these mutations are determined by the scale
        argument.
        """
        point_zero = numpy.asarray(point_zero)
        if point_zero.shape != self.problem.shape:
            raise ValueError('point_zero has an incorrect shape - {}, '
                'should be {}'.format(point_zero.shape, self.problem.shape))
        scale = numpy.asarray(scale)
        if scale.shape != self.problem.shape:
            raise ValueError('scale has an incorrect shape - {}, '
                'should be {}'.format(scale.shape, self.problem.shape))

        self.generation_number = 0
        self.population.clear()
        self.selection.clear()
        self.selection_fitness.clear()
        if self['elitism'] and self.problem.validate(point_zero):
            self.population.append(point_zero)
        while len(self.population) < self['population']:
            mutation = self.nprandom.normal(0.0, 1.0, self.problem.shape)
            mutation *= scale
            value = point_zero + mutation
            value = self.problem.clip(value)
            if self.problem.validate(value):
                self.population.append(value)

def rosenbrock(x, y):
    a = 1
    b = 100
    return (a - x)**2 + b * (y - x**2)**2

@minitest.test
def test_banana_min():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.minimize, lbound=(-2, -1), ubound=(2, 3))
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga.seed(45)
    ga.populate((0, 2), scale=0.5)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 10

@minitest.test
def test_banana_max():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.maximize, lbound=(-2, -1), ubound=(2, 3))
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga.seed(45)
    ga.populate((0, 2), scale=0.5)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] > 700

@minitest.test
def test_banana_elite():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.minimize, lbound=(-2, -1), ubound=(2, 3))
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga['elitism'] = 1
    ga.seed(45)
    ga.populate((0, 2), scale=0.5)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 10

@minitest.test
def test_banana_valid():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.maximize, lbound=(-2, -1), ubound=(2, 3))
    problem.set_validator(lambda x: x[1] > x[0]**2)
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga.seed(45)
    ga.populate((0, 2), scale=4.0)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 1500

@minitest.test
def test_banana_uniform():
    problem = OptimizationProblem((2,), lambda x: rosenbrock(x[0], x[1]),
            Goal.maximize, lbound=(-2, -1), ubound=(2, 3))
    problem.set_validator(lambda x: x[1] > x[0]**2)
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga['population'] = 10
    ga.seed(45)
    ga.uniform()
    assert len(ga.population) == 10
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 1500

@minitest.test
def test_half_bounded():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize, lbound=[0]*20)
    ga = BasicGaDriver(problem, [1.0] * 20)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_size'] = 0.1
    ga.seed(45)
    ga.populate([1] * 20)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 5

@minitest.test
def test_unbounded():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = BasicGaDriver(problem, [1.0] * 20)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_size'] = 0.1
    ga.seed(45)
    ga.populate([0] * 20)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 15

@minitest.test
def test_mutation():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = BasicGaDriver(problem, [1.0] * 20)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_probability'] = 0.1
    ga['mutation_count'] = 2
    ga.seed(45)
    ga.populate([0] * 20)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 10

@minitest.test
def test_mutation_zero_1():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = BasicGaDriver(problem, [1.0] * 20)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_probability'] = 0.0
    ga['mutation_count'] = 0
    ga['crossover_rate'] = 0
    ga.seed(45)
    ga.populate([0] * 20)
    ga.generate()
    optimal = ga.selection_fitness[0]
    for i in range(3):
        ga.generate()
        assert ga.selection_fitness[0] == optimal

@minitest.test
def test_mutation_zero_2():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = BasicGaDriver(problem, [1.0] * 20)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_size'] = 0
    ga['crossover_rate'] = 0
    ga.seed(45)
    ga.populate([0] * 20)
    ga.generate()
    optimal = ga.selection_fitness[0]
    for i in range(3):
        ga.generate()
        assert ga.selection_fitness[0] == optimal

@minitest.test
def test_multidim():
    problem = OptimizationProblem((2, 3, 4), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = BasicGaDriver(problem, numpy.ones((2, 3, 4)))
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_size'] = 0.2
    ga.seed(45)
    ga.populate(numpy.zeros((2, 3, 4)))
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 15

@minitest.test
def test_differential():
    problem = OptimizationProblem((20,), lambda x: numpy.sum(x**2),
            Goal.minimize)
    ga = DifferentialGaDriver(problem)
    ga['population'] = 20
    ga['selection'] = 5
    ga['mutation_size'] = 1.0
    ga.seed(45)
    ga.populate([0] * 20, [1] * 20)
    for i in range(10):
        ga.generate()
    assert ga.selection_fitness[0] < 5

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
    ga = BasicGaDriver(problem, (0.25, 0.5))
    ga['population'] = 10
    ga.seed(45)
    ga.populate((0, 2), scale=0.5)
    for i in range(2):
        ga.generate()
    try:
        ga.generate()
    except RuntimeError:
        pass
    else:
        assert False

    for i in range(8):
        ga.generate()

    assert ga.generation_number == 10
    assert eval_counter == 101
    assert ga.selection_fitness[0] < 10

@minitest.runner
def test(interactive):
    pass
