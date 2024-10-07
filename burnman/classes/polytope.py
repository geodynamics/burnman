# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import importlib
import numpy as np
from fractions import Fraction
from scipy.spatial import Delaunay
from scipy.special import comb
from copy import copy
import cdd as cdd_float

from .material import cached_property

from ..utils.math import independent_row_indices


try:
    cdd_fraction = importlib.import_module("cdd.gmp")
    cdd_gmp_loaded = True
except ImportError:
    cdd_fraction = importlib.import_module("cdd")
    cdd_gmp_loaded = False


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
        assert vertices >= 2, "need at least two vertices"
        assert points_per_edge >= 2, "need at least 2 points per edge"

        self.vertices = vertices
        self.points_per_edge = points_per_edge

    def generate(self, generate_type="list"):
        """
        Generates the grid points of the simplex in lexicographic order.

        :param generate_type: Determines whether the generator returns
            lists or arrays corresponding to each point in the simplex grid.
            Valid options are 'list' or 'array'.
        :type generate_type: str

        :returns: Grid points of the simplex.
        :rtype: generator of lists or ndarrays (int, ndim=1)
        """

        if generate_type == "list":
            x = [0] * self.vertices
        elif generate_type == "array":
            x = np.zeros(self.vertices, dtype=int)
        else:
            raise Exception("generate_type must be of type list or array.")

        x[self.vertices - 1] = self.points_per_edge - 1

        h = self.vertices
        while True:
            yield copy(x)

            h -= 1
            if h == 0:
                return

            val = x[h]
            x[h] = 0
            x[self.vertices - 1] = val - 1
            x[h - 1] += 1
            if val != 1:
                h = self.vertices

    def grid(self, generate_type="list"):
        """
        Returns either a list or a numpy array
        corresponding the the points in the simplex grid, depending on
        whether the user chooses 'list' (default) or 'array' as
        the generate_type parameter.
        """
        if generate_type == "list":
            return list(self.generate(generate_type))
        elif generate_type == "array":
            return np.array(list(self.generate(generate_type)))
        else:
            raise Exception("generate_type must be of type list or array.")

    def n_points(self):
        """
        The number of points corresponding to the number of vertices and
        points per edge chosen by the user.
        """
        return comb(
            self.vertices + self.points_per_edge - 2, self.vertices - 1, exact=True
        )


class MaterialPolytope(object):
    """
    A class for creating and manipulating polytope objects using pycddlib.

    This class is available as :class:`burnman.polytope.MaterialPolytope`.
    """

    def __init__(
        self,
        equalities,
        inequalities,
        return_fractions=False,
        independent_endmember_occupancies=None,
    ):
        """
        Initialization function for the MaterialPolytope class.
        Declares basis attributes of the class.

        :param equalities: A numpy array containing all the
            equalities defining the polytope. Each row should evaluate to 0.
        :type equalities: numpy.array (2D) of floats or Fractions
        :param inequalities: A numpy array containing all the inequalities
            defining the polytope. Each row should evaluate to <= 0.
        :type inequalities: numpy.array (2D) of the same type as equalities
        :param return_fractions: Whether or not to return fractions.
        :type return_fractions: bool
        :param independent_endmember_occupancies: If specified, this array
            provides the independent endmember set against which the
            dependent endmembers are defined.
        :type independent_endmember_occupancies: numpy.array (2D) or None
        """
        if equalities.dtype != inequalities.dtype:
            raise Exception(
                f"The equalities and inequalities arrays should have the same type ({equalities.dtype} != {inequalities.dtype})."
            )

        self.set_return_type(return_fractions)
        self.equality_matrix = equalities[:, 1:]
        self.equality_vector = -equalities[:, 0]

        if equalities.dtype == Fraction and cdd_gmp_loaded is True:
            self.number_type = Fraction
            self.cdd = cdd_fraction
        elif equalities.dtype == float or cdd_gmp_loaded is False:
            self.number_type = float
            self.cdd = cdd_float
        else:
            raise Exception(
                "equalities should be arrays of either floats or Fractions."
            )

        self.polytope_matrix = self.cdd.matrix_from_array(equalities)
        self.polytope_matrix.lin_set = set(range(len(equalities)))
        self.polytope_matrix.rep_type = self.cdd.RepType.INEQUALITY

        self.cdd.matrix_append_to(
            self.polytope_matrix, self.cdd.matrix_from_array(inequalities)
        )
        self.polytope = self.cdd.polyhedron_from_matrix(self.polytope_matrix)

        if independent_endmember_occupancies is not None:
            self.independent_endmember_occupancies = independent_endmember_occupancies

    def set_return_type(self, return_fractions=False):
        """
        Sets the return_type for the polytope object. Also deletes the cached
        endmember_occupancies property.

        :param return_fractions: Choose whether the generated polytope object
            should return fractions or floats.
        :type return_fractions: bool
        """
        try:
            del self.__dict__["endmember_occupancies"]
        except KeyError:
            pass
        self.return_fractions = return_fractions

    @cached_property
    def raw_vertices(self):
        """
        Returns a list of the vertices of the polytope without any
        postprocessing. See also endmember_occupancies.
        """
        return self.cdd.copy_generators(self.polytope).array

    @cached_property
    def limits(self):
        """
        Return the limits of the polytope (the set of bounding inequalities).
        """
        return np.array(self.cdd.copy_inequalities(self.polytope).array, dtype=float)

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
        if self.number_type == float and self.return_fractions:
            v = np.array(
                [
                    [Fraction(value).limit_denominator(1000000) for value in v]
                    for v in self.raw_vertices
                ]
            )
        elif self.number_type == Fraction and not self.return_fractions:
            v = np.array(self.raw_vertices).astype(float)
        else:
            v = np.array(self.raw_vertices)

        if len(v.shape) == 1:
            raise ValueError(
                "The combined equality and positivity "
                "constraints result in a null polytope."
            )
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

        sol = (
            np.linalg.lstsq(
                np.array(ind.T).astype(float),
                np.array(self.endmember_occupancies.T).astype(float),
                rcond=0,
            )[0]
            .round(decimals=12)
            .T
        )
        return sol

    def _decompose_vertices_into_simplices(self, vertices):
        """
        Decomposes a set of vertices into simplices by Delaunay triangulation.
        """
        # Delaunay triangulation only works in dimensions > 1
        # and we remove the nullspace (sum(fractions) = 1)
        if len(vertices) > 2:
            nulls = np.repeat(vertices[:, -1], vertices.shape[1]).reshape(
                vertices.shape
            )
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
        arr = [[Fraction(value).limit_denominator(1000000) for value in v] for v in arr]
        M = cdd_fraction.matrix_from_array(arr)
        M.rep_type = cdd_fraction.RepType.GENERATOR
        return cdd_fraction.polyhedron_from_matrix(M)

    @cached_property
    def independent_endmember_limits(self):
        """
        Gets the limits of the polytope as a function of the independent
        endmembers.
        """
        ind_poly = self.independent_endmember_polytope
        inequalities = cdd_fraction.copy_inequalities(ind_poly).array
        return np.array(inequalities, dtype=float)

    def subpolytope_from_independent_endmember_limits(self, limits):
        """
        Returns a smaller polytope by applying additional limits to the amounts
        of the independent endmembers.
        """
        ind_poly = self.independent_endmember_polytope
        modified_limits = cdd_fraction.copy_inequalities(ind_poly)
        cdd_fraction.matrix_append_to(
            modified_limits, cdd_fraction.matrix_from_array(limits)
        )
        return cdd_fraction.polyhedron_from_matrix(modified_limits)

    def subpolytope_from_site_occupancy_limits(self, limits):
        """
        Returns a smaller polytope by applying additional limits to the
        individual site occupancies.
        """
        modified_limits = copy(self.polytope_matrix)
        self.cdd.matrix_append_to(modified_limits, self.cdd.matrix_from_array(limits))
        return self.cdd.polyhedron_from_matrix(modified_limits)

    def grid(
        self,
        points_per_edge=2,
        unique_sorted=True,
        grid_type="independent endmember proportions",
        limits=None,
    ):
        """
        Create a grid of points which span the polytope.

        :param points_per_edge: Number of points per edge of the polytope.
        :type points_per_edge: int
        :param unique_sorted: The gridding is done by splitting the polytope
            into a set of simplices. This means that points will be duplicated
            along vertices, faces etc. If unique_sorted is True, this function
            will sort and make the points unique. This is an expensive
            operation for large polytopes, and may not always be necessary.
        :type unique_sorted: bool
        :param grid_type: Whether to grid the polytope in terms of
            independent endmember proportions or site occupancies.
            Choices are 'independent endmember proportions' or
            'site occupancies'
        :type grid_type: str
        :param limits: Additional inequalities restricting the
            gridded area of the polytope.
        :type limits: numpy.array (2D)

        :returns: A list of points gridding the polytope.
        :rtype: numpy.array (2D)
        """
        if limits is None:
            if grid_type == "independent endmember proportions":
                f_occ = self.endmembers_as_independent_endmember_amounts / (
                    points_per_edge - 1
                )
            elif grid_type == "site occupancies":
                f_occ = self.endmember_occupancies / (points_per_edge - 1)
            else:
                raise Exception(
                    "grid type not recognised. Should be one of "
                    "independent endmember proportions "
                    "or site occupancies"
                )

            simplices = self._decompose_vertices_into_simplices(
                self.endmembers_as_independent_endmember_amounts
            )
        else:
            if grid_type == "independent endmember proportions":
                plims = self.subpolytope_from_site_occupancy_limits(limits)
                ppns = np.array(self.cdd.copy_generators(plims).array)[:, 1:]
                last_ppn = np.array([1.0 - sum(p) for p in ppns]).reshape(
                    (len(ppns), 1)
                )
                vertices_as_independent_endmember_proportions = np.hstack(
                    (ppns, last_ppn)
                )
                f_occ = vertices_as_independent_endmember_proportions / (
                    points_per_edge - 1
                )

            elif grid_type == "site occupancies":
                plims = self.subpolytope_from_site_occupancy_limits(limits)
                occ = np.array(self.cdd.copy_generators(plims).array)[:, 1:]
                f_occ = occ / (points_per_edge - 1)

                ind = self.independent_endmember_occupancies

                vertices_as_independent_endmember_proportions = (
                    np.linalg.lstsq(
                        np.array(ind.T).astype(float),
                        np.array(occ.T).astype(float),
                        rcond=None,
                    )[0]
                    .round(decimals=12)
                    .T
                )
            else:
                raise Exception(
                    "grid_type not recognised. "
                    "Should be one of "
                    "independent endmember proportions "
                    "or site occupancies"
                )

            simplices = self._decompose_vertices_into_simplices(
                vertices_as_independent_endmember_proportions
            )

        n_ind = f_occ.shape[1]
        n_simplices = len(simplices)
        dim = len(simplices[0])
        simplex_grid = SimplexGrid(dim, points_per_edge)
        grid = simplex_grid.grid("array")
        points_per_simplex = simplex_grid.n_points()
        n_points = n_simplices * points_per_simplex

        points = np.empty((n_points, n_ind))
        idx = 0
        for i in range(0, n_simplices):
            points[idx : idx + points_per_simplex] = grid.dot(f_occ[simplices[i]])
            idx += points_per_simplex

        if unique_sorted:
            points = np.unique(points, axis=0)
        return points
