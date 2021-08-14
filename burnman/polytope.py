# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import cdd
import numpy as np
from sympy import Matrix, Rational
from fractions import Fraction
from scipy.spatial import Delaunay
from scipy.special import comb
from copy import copy

from .reductions import row_reduce
from .material import cached_property


def independent_row_indices(array):
    """

    """
    m = Matrix(array.shape[0], array.shape[1],
               lambda i, j: Rational(array[i, j]).limit_denominator(1000))
    _, pivots, swaps = row_reduce(m, iszerofunc=lambda x: x.is_zero,
                                  simpfunc=lambda x: Rational(x).limit_denominator(1000))
    indices = np.array(range(len(array)))
    for swap in np.array(swaps):
        indices[swap] = indices[swap[::-1]]
    return indices[:len(pivots)]


def generate_complete_basis(incomplete_basis, complete_basis):
    """

    """
    a = np.concatenate((incomplete_basis, complete_basis))
    return a[independent_row_indices(a)]


def complete_basis(basis):
    """
    # Creates a full basis by filling remaining rows with
    # rows of the identity matrix with row indices not
    # in the column pivot list of the basis RREF
    """

    n, m = basis.shape
    if n < m:
        pivots = list(Matrix(basis).rref()[1])
        return np.concatenate((basis,
                               np.identity(m)[[i for i in range(m)
                                               if i not in pivots], :]),
                              axis=0)
    else:
        return basis


class SimplexGrid(object):
    def __init__(self, vertices, points_per_edge):

        assert vertices >= 2, 'need at least two vertices'
        assert points_per_edge >= 2, 'need at least 2 points per edge'

        self.vertices = vertices
        self.points_per_edge = points_per_edge

    def generate(self, generate_type='list'):
        """
        Generate the grid points of the simplex in lexicographic order.
        Returns
        -------
        generator of ndarrays (int, ndim=1)
            Grid points of the simplex.
        """

        if generate_type == 'list':
            x = [0]*self.vertices
        elif generate_type == 'array':
            x = np.zeros(self.vertices, dtype=int)
        else:
            raise Exception('generate should create either lists or '
                            'arrays of indices')

        x[self.vertices-1] = self.points_per_edge-1

        h = self.vertices
        while True:
            yield copy(x)

            h -= 1
            if h == 0:
                return

            val = x[h]
            x[h] = 0
            x[self.vertices-1] = val - 1
            x[h-1] += 1
            if val != 1:
                h = self.vertices

    def grid(self, generate_type='list'):
        """

        """
        if generate_type == 'list':
            return list(self.generate(generate_type))
        else:
            return np.array(list(self.generate(generate_type)))

    def n_points(self):
        """

        """
        return comb(self.vertices+self.points_per_edge-2,
                    self.vertices-1, exact=True)


class MaterialPolytope(object):
    """

    """

    def __init__(self, equalities,
                 inequalities,
                 number_type='fraction',
                 return_fractions=False,
                 independent_endmember_occupancies=None):
        self.return_fractions = return_fractions
        self.equality_matrix = equalities[:, 1:]
        self.equality_vector = -equalities[:, 0]

        self.polytope_matrix = cdd.Matrix(equalities, linear=True,
                                          number_type=number_type)
        self.polytope_matrix.rep_type = cdd.RepType.INEQUALITY
        self.polytope_matrix.extend(inequalities, linear=False)
        self.polytope = cdd.Polyhedron(self.polytope_matrix)

        if independent_endmember_occupancies is not None:
            self.independent_endmember_occupancies = independent_endmember_occupancies

    def set_return_type(self, return_fractions=False):
        """

        """
        try:
            del self.__dict__['endmember_occupancies']
        except KeyError:
            pass
        self.return_fractions = return_fractions

    @cached_property
    def raw_vertices(self):
        """

        """
        return self.polytope.get_generators()[:]

    @cached_property
    def site_occupancy_limits(self):
        """

        """
        return np.array(self.polytope.get_inequalities(), dtype=float)

    @cached_property
    def n_endmembers(self):
        """

        """
        return len(self.raw_vertices)

    @cached_property
    def endmember_occupancies(self):
        """

        """
        if self.return_fractions:
            if self.polytope.number_type == 'fraction':
                v = np.array([[Fraction(value) for value in v]
                              for v in self.raw_vertices])
            else:
                v = np.array([[Rational(value).limit_denominator(1000000)
                               for value in v]
                              for v in self.raw_vertices])
        else:
            v = np.array([[float(value) for value in v]
                          for v in self.raw_vertices])

        if len(v.shape) == 1:
            raise ValueError("The combined equality and positivity "
                             "constraints result in a null polytope.")

        return v[:, 1:] / v[:, 0, np.newaxis]

    @cached_property
    def independent_endmember_occupancies(self):
        """

        """
        arr = self.endmember_occupancies
        return arr[independent_row_indices(arr)]

    def independent_endmember_proportions(self, endmember_occupancies):
        """

        """
        ind = self.independent_endmember_occupancies
        return np.array(Matrix(ind.T).pinv_solve(Matrix(endmember_occupancies.T)).T)

    @cached_property
    def endmembers_as_independent_endmember_proportions(self):
        """

        """
        ind = self.independent_endmember_occupancies

        sol = np.linalg.lstsq(np.array(ind.T).astype(float),
                              np.array(self.endmember_occupancies.T).astype(
                                  float),
                              rcond=0)[0].round(decimals=12).T
        return sol

    def _decompose_polytope_into_endmember_simplices(self, vertices):
        """

        """
        # Delaunay triangulation only works in dimensions > 1
        # and we remove the nullspace (sum(fractions) = 1)
        if len(vertices) > 2:
            nulls = np.repeat(vertices[:, -1],
                              vertices.shape[1]).reshape(vertices.shape)
            tri = Delaunay((vertices - nulls)[:, :-1])
            return tri.simplices
        else:
            return [[0, 1]]

    @cached_property
    def independent_endmember_polytope(self):
        """
        The polytope involves the first n-1 independent endmembers.
        The last endmember proportion makes the sum equal to one.
        """
        arr = self.endmembers_as_independent_endmember_proportions
        arr = np.hstack((np.ones((len(arr), 1)), arr[:, :-1]))
        M = cdd.Matrix(arr, number_type='fraction')
        M.rep_type = cdd.RepType.GENERATOR
        return cdd.Polyhedron(M)

    @cached_property
    def independent_endmember_limits(self):
        """

        """
        return np.array(self.independent_endmember_polytope.get_inequalities(),
                        dtype=float)

    def subpolytope_from_independent_endmember_limits(self, limits):
        """

        """
        modified_limits = self.independent_endmember_polytope.get_inequalities().copy()
        modified_limits.extend(limits, linear=False)
        return cdd.Polyhedron(modified_limits)

    def subpolytope_from_site_occupancy_limits(self, limits):
        """

        """
        modified_limits = self.polytope_matrix.copy()
        modified_limits.extend(limits, linear=False)
        return cdd.Polyhedron(modified_limits)

    def grid(self, points_per_edge=2, unique_sorted=True,
             grid_type='independent endmember proportions', limits=None):
        """
        """
        if limits is None:
            if grid_type == 'independent endmember proportions':
                f_occ = self.endmembers_as_independent_endmember_proportions / \
                    (points_per_edge-1)
            elif grid_type == 'site occupancies':
                f_occ = self.endmember_occupancies/(points_per_edge-1)
            else:
                raise Exception(
                    'grid type not recognised. Should be one of independent endmember proportions or site occupancies')

            simplices = self._decompose_polytope_into_endmember_simplices(
                vertices=self.endmembers_as_independent_endmember_proportions)
        else:
            if grid_type == 'independent endmember proportions':
                ppns = np.array(self.subpolytope_from_independent_endmember_limits(
                    limits).get_generators()[:])[:, 1:]
                last_ppn = np.array([1. - sum(p)
                                    for p in ppns]).reshape((len(ppns), 1))
                vertices_as_independent_endmember_proportions = np.hstack(
                    (ppns, last_ppn))
                f_occ = vertices_as_independent_endmember_proportions / \
                    (points_per_edge-1)

            elif grid_type == 'site occupancies':
                occ = np.array(self.subpolytope_from_site_occupancy_limits(
                    limits).get_generators()[:])[:, 1:]
                f_occ = occ/(points_per_edge-1)

                ind = self.independent_endmember_occupancies

                vertices_as_independent_endmember_proportions = np.linalg.lstsq(np.array(ind.T).astype(float),
                                                                                np.array(occ.T).astype(
                                                                                    float),
                                                                                rcond=None)[0].round(decimals=12).T
            else:
                raise Exception('grid_type not recognised. '
                                'Should be one of '
                                'independent endmember proportions '
                                'or site occupancies')

            simplices = self._decompose_polytope_into_endmember_simplices(
                vertices=vertices_as_independent_endmember_proportions)

        n_ind = f_occ.shape[1]
        n_simplices = len(simplices)
        dim = len(simplices[0])
        simplex_grid = SimplexGrid(dim, points_per_edge)
        grid = simplex_grid.grid('array')
        points_per_simplex = simplex_grid.n_points()
        n_points = n_simplices*points_per_simplex

        points = np.empty((n_points, n_ind))
        idx = 0
        for i in range(0, n_simplices):
            points[idx:idx+points_per_simplex] = grid.dot(f_occ[simplices[i]])
            idx += points_per_simplex

        if unique_sorted:
            points = np.unique(points, axis=0)
        return points
