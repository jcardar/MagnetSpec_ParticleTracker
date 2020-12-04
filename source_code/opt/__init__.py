"""A framework for solving optimization problems.
"""

from .core import OptimizationProblem, Goal

def test(interactive):
    from . import ga, neldermead
    ga.test(interactive)
    neldermead.test(interactive)
