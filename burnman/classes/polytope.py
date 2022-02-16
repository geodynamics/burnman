# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import importlib
import numpy as np
from sympy import Rational
from fractions import Fraction
from scipy.spatial import Delaunay
from scipy.special import comb
from copy import copy

from .material import cached_property

from ..utils.math import independent_row_indices

try:
    cdd = importlib.import_module('cdd')
except ImportError as err:
    print(f'Warning: {err}. '
          'For full functionality of BurnMan, please install pycddlib.')


class SimplexGrid(object):
    """
    A class that creates objects that can efficiently generate a set of points
    that grid a simplex with a user-defined number of vertices. The class
    contains both a generator method and a grid method. It also contains
    an n_points attribute that returns the number of points in the gridded
    simplex.

    This class is available as :class:`burnman.polytope.SimplexGrid`.
    """

    def __init__(self, vertices, points_per_edge):
        """
        Initialize SimplexGrid object with the desired number of vertices
        and points per edge.
        """
        assert vertices >= 2, 'need at least two vertices'
        assert points_per_edge >= 2, 'need at least 2 points per edge'

        self.vertices = vertices
        self.points_per_edge = points_per_edge

    def generate(self, generate_type='list'):
        """
        Generates the grid points of the simplex in lexicographic order.

        Parameters
        ----------
        generate_type : 'list' or 'array'
            Determines whether the generator returns lists or arrays
            corresponding to each point in the simplex grid.

        Returns
        -------
        generator of lists or ndarrays (int, ndim=1)
            Grid points of the simplex.
        """

        if generate_type == 'list':
            x = [0]*self.vertices
        elif generate_type == 'array':
            x = np.zeros(self.vertices, dtype=int)
        else:
            raise Exception('generate_type must be of type list or array.')

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
        Returns either a list or a numpy array
        corresponding the the points in the simplex grid, depending on
        whether the user chooses 'list' (default) or 'array' as
        the generate_type parameter.
        """
        if generate_type == 'list':
            return list(self.generate(generate_type))
        elif generate_type == 'array':
            return np.array(list(self.generate(generate_type)))
        else:
            raise Exception('generate_type must be of type list or array.')

    def n_points(self):
        """
        The number of points corresponding to the number of vertices and
        points per edge chosen by the user.
        """
        return comb(self.vertices+self.points_per_edge-2,
                    self.vertices-1, exact=True)


class MaterialPolytope(object):
    """
    A class that can be instantiated to create pycddlib polytope objects.
    These objects can be interrogated to provide the vertices satisfying the
    input constraints.

    This class is available as :class:`burnman.polytope.MaterialPolytope`.
    """

    def __init__(self, equalities,
                 inequalities,
                 number_type='fraction',
                 return_fractions=False,
                 independent_endmember_occupancies=None):
        """
        Initialization function for the MaterialPolytope class.
        Declares basis attributes of the class.

        Parameters
        ----------
        equalities: 2D numpy array
            A numpy array containing all the equalities of the polytope.
            Each row should evaluate to 0.
        inequalities: 2D numpy array
            A numpy array containing all the inequalities of the polytope.
            Each row should evaluate to <= 0.
        number_type: 'fraction' or 'float' (default is 'fraction')
            Whether pycddlib should read the input arrays as
            fractions or floats.
        return_fractions : boolean (default is False)
            Whether the generated polytope object should return fractions or
            floats.
        independent_endmember_occupancies : 2D numpy array (or None)
            If specified, this array provides the independent endmember set
            against which the dependent endmembers are defined.
        """
        self.set_return_type(return_fractions)
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
        Sets the return_type for the polytope object. Also deletes the cached
        endmember_occupancies property.

        Parameters
        ----------
        return_fractions : boolean (default is False)
            Whether the generated polytope object should return fractions or
            floats.
        """
        try:
            del self.__dict__['endmember_occupancies']
        except KeyError:
            pass
        self.return_fractions = return_fractions

    @cached_property
    def raw_vertices(self):
        """
        Returns a list of the vertices of the polytope without any
        postprocessing. See also endmember_occupancies.
        """
        return self.polytope.get_generators()[:]

    @cached_property
    def limits(self):
        """
        Return the limits of the polytope (the set of bounding inequalities).
        """
        return np.array(self.polytope.get_inequalities(), dtype=float)

    @cached_property
    def n_endmembers(self):
        """
        Return the number of endmembers
        (the number of vertices of the polytope).
        """
        return len(self.raw_vertices)

    @cached_property
    def endmember_occupancies(self):
        """
        Return the endmember occupancies
        (a processed list of all of the vertex locations).
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
        Return an independent set of endmember occupancies
        (a linearly-independent set of vertex locations)
        """
        arr = self.endmember_occupancies
        return arr[independent_row_indices(arr)]

    @cached_property
    def endmembers_as_independent_endmember_amounts(self):
        """
        Return a list of all the endmembers as a linear sum of
        the independent endmembers.
        """
        ind = self.independent_endmember_occupancies

        sol = np.linalg.lstsq(np.array(ind.T).astype(float),
                              np.array(self.endmember_occupancies.T).astype(
                                  float),
                              rcond=0)[0].round(decimals=12).T
        return sol

    def _decompose_vertices_into_simplices(self, vertices):
        """
        Decomposes a set of vertices into simplices by Delaunay triangulation.
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
        Returns the polytope expressed in terms of proportions of the
        independent endmembers. The polytope involves the first
        n-1 independent endmembers. The last endmember proportion makes
        the sum equal to one.
        """
        arr = self.endmembers_as_independent_endmember_amounts
        arr = np.hstack((np.ones((len(arr), 1)), arr[:, :-1]))
        M = cdd.Matrix(arr, number_type='fraction')
        M.rep_type = cdd.RepType.GENERATOR
        return cdd.Polyhedron(M)

    @cached_property
    def independent_endmember_limits(self):
        """
        Gets the limits of the polytope as a function of the independent
        endmembers.
        """
        return np.array(self.independent_endmember_polytope.get_inequalities(),
                        dtype=float)

    def subpolytope_from_independent_endmember_limits(self, limits):
        """
        Returns a smaller polytope by applying additional limits to the amounts
        of the independent endmembers.
        """
        modified_limits = self.independent_endmember_polytope.get_inequalities().copy()
        modified_limits.extend(limits, linear=False)
        return cdd.Polyhedron(modified_limits)

    def subpolytope_from_site_occupancy_limits(self, limits):
        """
        Returns a smaller polytope by applying additional limits to the
        individual site occupancies.
        """
        modified_limits = self.polytope_matrix.copy()
        modified_limits.extend(limits, linear=False)
        return cdd.Polyhedron(modified_limits)

    def grid(self, points_per_edge=2, unique_sorted=True,
             grid_type='independent endmember proportions', limits=None):
        """
        Create a grid of points which span the polytope.

        Parameters
        ----------
        points_per_edge : integer (default is 2)
            Number of points per edge of the polytope.
        unique_sorted : boolean (default is True)
            The gridding is done by splitting the polytope into
            a set of simplices. This means that points will be duplicated along
            vertices, faces etc. If unique_sorted is True, this function
            will sort and make the points unique. This is an expensive
            operation for large polytopes, and may not always be necessary.
        grid_type : 'independent endmember proportions' (default) or 'site occupancies'
            Whether to grid the polytope in terms of
            independent endmember proportions or site occupancies.
        limits : 2D numpy array
            Additional inequalities restricting the gridded area of the polytope.
        Returns
        -------
        points : 2D numpy array
            A list of points gridding the polytope.
        """
        if limits is None:
            if grid_type == 'independent endmember proportions':
                f_occ = (self.endmembers_as_independent_endmember_amounts
                         / (points_per_edge - 1))
            elif grid_type == 'site occupancies':
                f_occ = self.endmember_occupancies/(points_per_edge-1)
            else:
                raise Exception('grid type not recognised. Should be one of '
                                'independent endmember proportions '
                                'or site occupancies')

            simplices = self._decompose_vertices_into_simplices(
                self.endmembers_as_independent_endmember_amounts)
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

            simplices = self._decompose_vertices_into_simplices(
                vertices_as_independent_endmember_proportions)

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
