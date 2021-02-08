from enum import Enum
import numpy

class Goal(Enum):
    maximize = 0
    minimize = 1

class OptimizationProblem:
    """The specification for an optimization problem

    This is the description of an optimization problem, consisting of the
    dimensionality of the parameter space (in the form of a numpy-style array
    shape), optional upper and/or lower bounds, a goal function, and a
    specification of whether to maximize or minimize the goal function. It is
    also possible to specify a validator, which can be used to reject some
    points early.

    An instance generally created by the user, then passed to one or more
    drivers. The drivers then use the public API to access the problem
    description and perform function evaluations.

    The solver state is stored in the driver, rather than in here, so it needs
    to be moved over manually when switching drivers.
    """
    def __init__(self, shape, goal_func, goal_type, *,
            lbound=None, ubound=None, validator=None):
        self.shape = shape
        self.goal_func = goal_func
        self.goal = goal_type

        self.set_bounds(lbound, ubound)
        self.set_validator(validator)
        self.set_callback(None)

        self.fevals = 0

    def set_bounds(self, lbound, ubound):
        """Sets the upper and lower bounds.
        
        These must be array-like with the correct shape. Specifying None for
        either the upper or lower bound will remove that bound, if it already
        exists.
        """
        if lbound is not None:
            lbound = numpy.asarray(lbound)
            if lbound.shape != self.shape:
                raise ValueError('lbound has an incorrect shape - {}, should '
                        'be {}'.format(lbound.shape, self.shape))

        if ubound is not None:
            ubound = numpy.asarray(ubound)
            if ubound.shape != self.shape:
                raise ValueError('ubound has an incorrect shape - {}, should '
                        'be {}'.format(ubound.shape, self.shape))

        self.lbound = lbound
        self.ubound = ubound

    def set_validator(self, validator):
        """Sets the validator function.

        The validator function provides a user-programmable way of checking
        points for validity, before attempting to evaluate them. It provides a
        way of communicating to the optimization algorithm that a set of points
        is outside the domain of the goal function. The validator function will
        be passed the point, which will be inside any bounds set by set_bounds,
        and must return True for valid points and False for invalid ones.
        Passing None to set_validator removes any existing validator (all points
        are considered valid)
        """
        if validator is None:
            self.validator = lambda _: True
        else:
            self.validator = validator
    
    def set_callback(self, callback):
        """Sets the evaluation callback.

        The evaluation callback is called every time the goal function is
        evaluated. It is passed the point at which the function was evaluated
        and the function value. Call set_callback(None) to remove any existing
        callback.
        """
        self.callback = callback

    def bounds(self):
        """Returns the upper and lower bounds, as a length-2 tuple."""
        return (self.lbound, self.ubound)

    def clip(self, value):
        """Clips a given value to the upper and lower bounds."""
        if self.lbound is not None:
            value = numpy.maximum(value, self.lbound)
        if self.ubound is not None:
            value = numpy.minimum(value, self.ubound)
        return value

    def __check_shape(self, value):
        if value.shape != self.shape:
            raise ValueError('incorrect shape - {}, should be {}'.format(
                value.shape, self.shape))

    def validate(self, value):
        """Runs the validation function on the given point.
        
        Return True for a valid point (inside any bounds, and validator returned
        True) and False otherwise.
        """
        self.__check_shape(value)
        if self.lbound is not None:
            if numpy.any(value < self.lbound):
                return False
        if self.ubound is not None:
            if numpy.any(value > self.ubound):
                return False

        return self.validator(value)

    def evaluate(self, value):
        """Evaluates the goal function at a point.

        If the supplied point is not valid (validate(point) returns False) this
        function returns None. Otherwise, it evaluates the goal function and
        returns that.
        """
        if not self.validate(value):
            return None
        val = self.goal_func(value)
        self.fevals += 1
        if self.callback is not None:
            self.callback(value, val)
        return val

    def evaluate_min(self, value):
        """Changes the problem to a minimization one.

        Both minimization and maximization problems are supported, and in
        principle the driver must identify which kind of problem is being solved
        and act appropriately. This function allows every problem to be treated
        as a minimization one, simply by changing the sign of the goal function
        for maximization problems.

        If other problem types are added in the future, use of this function
        with those problems are likely to raise an exception.
        """
        val = self.evaluate(value)
        if val is not None:
            if self.goal == Goal.minimize:
                return val
            elif self.goal == Goal.maximize:
                return -val
            else:
                raise ValueError('unknown Goal: {}'.format(self.goal))

    def reset_feval_counter(self):
        """Resets the function evaluation counter."""
        self.fevals = 0
