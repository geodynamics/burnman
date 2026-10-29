# This file is part of BurnMan, released under GNU GPL v2 or later.

import numpy as np
from scipy.stats import t
from util import BurnManTest
from burnman.optimize.nonlinear_fitting import (
    NonLinearModel,
    confidence_prediction_bands,
)


class LineModel(NonLinearModel):
    def __init__(self, scale=1.0):
        x = np.linspace(-2.0, 3.0, 20)
        y = 2.5 * x - 1.3 + np.random.RandomState(7).normal(scale=0.3, size=len(x))
        self.data = np.column_stack([x * scale, y])
        design = np.column_stack([x, np.ones(len(x))])
        self.params = np.linalg.solve(design.T @ design, design.T @ y) / [scale, 1.0]
        self.dof = len(x) - 2
        residual = y - design @ (self.params * [scale, 1.0])
        self.output_variance = residual @ residual / self.dof
        transform = np.diag([1.0 / scale, 1.0])
        self.pcov = (
            transform
            @ np.linalg.inv(design.T @ design)
            @ transform
            * self.output_variance
        )
        self.delta_params = np.array([1.0e-5 / scale, 1.0e-5])
        self.noise_variance = self.output_variance / (1.0 + self.params[0] ** 2)
        self.flags = [None] * len(x)

    def get_params(self):
        return self.params

    def set_params(self, p):
        self.params[:] = p

    def function(self, x, flag=None):
        return np.array([x[0], self.params[0] * x[0] + self.params[1]])

    def normal(self, x, flag=None):
        n = np.array([self.params[0], -1.0])
        return n / np.linalg.norm(n)


class test_confidence_prediction_bands(BurnManTest):
    def test_closed_form_regression_limits(self):
        model = LineModel()
        query = np.array([[-1, 999], [0, -999], [2, 0]])
        actual = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        X = np.column_stack([query[:, 0], np.ones(len(query))])
        mean = X @ model.params
        variance = np.einsum("ij,jk,ik->i", X, model.pcov, X)
        critical = t.ppf(0.975, model.dof)
        expected = np.array(
            [
                mean - critical * np.sqrt(variance),
                mean + critical * np.sqrt(variance),
                mean - critical * np.sqrt(variance + model.output_variance),
                mean + critical * np.sqrt(variance + model.output_variance),
            ]
        )
        self.assertArraysAlmostEqual(
            actual.flatten(), expected.flatten(), tol=1.0e-7, tol_zero=1.0e-10
        )
        self.assertTrue(np.all(actual[2] <= actual[0]))
        self.assertTrue(np.all(actual[1] <= actual[3]))

    def test_intervals_coordinate_and_output_units(self):
        query = np.array([[-1.0, 0.0], [2.0, 0.0]])
        other = query.copy()
        other[:, 0] *= 1.0e-9
        baseline = confidence_prediction_bands(
            LineModel(), query, 0.95, lambda x: x[-1]
        )
        scaled = confidence_prediction_bands(
            LineModel(1.0e-9), other, 0.95, lambda x: x[-1]
        )
        self.assertArraysAlmostEqual(
            baseline.flatten(), scaled.flatten(), tol=1.0e-7, tol_zero=1.0e-10
        )
        model = LineModel()
        standard = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        for factor in (2.0, 1000.0, -2.0):
            with self.subTest(factor=factor):
                converted = confidence_prediction_bands(
                    model, query, 0.95, lambda x: factor * x[-1] + 7.0
                )
                restored = (converted - 7.0) / factor
                if factor < 0.0:
                    restored = restored[[1, 0, 3, 2]]
                self.assertArraysAlmostEqual(
                    restored.flatten(), standard.flatten(), tol=1.0e-7
                )

    def test_unrelated_output_keeps_confidence_bands(self):
        model = LineModel()
        query = np.array([[-1.0, 0.0], [0.0, 0.0], [2.0, 0.0]])
        # The fitted slope has parameter uncertainty, but scatter in measured
        # y does not determine the noise of a new slope measurement.
        actual = confidence_prediction_bands(
            model, query, 0.95, lambda x: model.get_params()[0]
        )
        half_width = t.ppf(0.975, model.dof) * np.sqrt(model.pcov[0, 0])
        self.assertArraysAlmostEqual(actual[0], [model.params[0] - half_width] * 3)
        self.assertArraysAlmostEqual(actual[1], [model.params[0] + half_width] * 3)
        self.assertTrue(np.all(np.isnan(actual[2:])))

    def test_noise_uses_only_matching_response_flags(self):
        model = LineModel()
        query = np.array([[-1.0, 0.0], [2.0, 0.0]])
        expected = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        # A second fitted quantity can have arbitrary units. It must neither supply
        # residuals nor consume the first quantity's residual degrees of freedom.
        count = len(model.data)
        other = model.data.copy()
        other[:, -1] *= 1.0e9
        model.data = np.vstack((model.data, other))
        model.flags.extend(["other"] * count)
        model.dof += count
        original_function = model.function

        def function(x, flag=None):
            self.assertIsNone(flag)
            return original_function(x, flag)

        model.function = function
        actual = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        # The coverage factor still uses the fit's global degrees of freedom.
        centre = np.mean(expected[:2], axis=0)
        ratio = t.ppf(0.975, model.dof) / t.ppf(0.975, model.dof - count)
        expected = centre + (expected - centre) * ratio
        self.assertArraysAlmostEqual(actual.flatten(), expected.flatten())

    def test_no_matching_observations_keeps_confidence_bands(self):
        model = LineModel()
        model.flags = ["other"] * len(model.data)
        actual = confidence_prediction_bands(
            model, model.data[:2], 0.95, lambda x: x[-1]
        )
        self.assertTrue(np.all(np.isfinite(actual[:2])))
        self.assertTrue(np.all(np.isnan(actual[2:])))

    def test_no_additional_model_attributes_are_required(self):
        model = LineModel()
        query = model.data[:2]
        expected = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        del model.noise_variance
        del model.flags
        actual = confidence_prediction_bands(model, query, 0.95, lambda x: x[-1])
        self.assertArraysAlmostEqual(actual.flatten(), expected.flatten())

    def test_parameter_bounds_preserve_output_noise_conversion(self):
        for fraction in (0.0, 0.25):
            with self.subTest(remaining_step=fraction):
                model = LineModel()
                query = model.data[:2]
                expected = confidence_prediction_bands(
                    model, query, 0.95, lambda x: x[-1]
                )
                upper = model.get_params().copy() + fraction * model.delta_params
                original_setter = model.set_params
                model.set_params = lambda p: original_setter(np.minimum(p, upper))
                actual = confidence_prediction_bands(
                    model, query, 0.95, lambda x: 2.0 * x[-1]
                )
                self.assertArraysAlmostEqual(
                    (actual / 2.0).flatten(), expected.flatten(), tol=1.0e-7
                )

    def test_fixed_parameters_cannot_identify_response_conversion(self):
        model = LineModel()
        original = model.get_params().copy()
        model.pcov[:] = 0.0
        model.set_params = lambda p: None
        actual = confidence_prediction_bands(
            model, model.data[:2], 0.95, lambda x: x[-1]
        )
        self.assertTrue(np.all(np.isfinite(actual[:2])))
        self.assertTrue(np.all(np.isnan(actual[2:])))
        self.assertArraysAlmostEqual(model.get_params(), original, tol=0.0)

    def test_parameter_restoration_on_callback_error(self):
        model = LineModel()
        params = model.get_params().copy()

        def output(point):
            if not np.array_equal(model.get_params(), params):
                raise ValueError("evaluation failed")
            return point[-1]

        with self.assertRaises(ValueError):
            confidence_prediction_bands(model, model.data[:2], 0.95, output)
        self.assertArraysAlmostEqual(model.get_params(), params, tol=0.0, tol_zero=0.0)

    def test_validation_and_empty_queries(self):
        model = LineModel()
        for probability in (0.0, 1.0, np.nan):
            with self.assertRaises(ValueError):
                confidence_prediction_bands(
                    model, model.data[:2], probability, lambda x: x[-1]
                )
        self.assertEqual(
            confidence_prediction_bands(
                model, np.empty((0, 2)), 0.95, lambda x: x[-1]
            ).shape,
            (4, 0),
        )
