import unittest
from util import BurnManTest
import tempfile
import os
import numpy as np
from scipy.ndimage import gaussian_filter
from sympy import Matrix, Rational
import scipy.integrate as integrate

from burnman.utils.misc import extract_lines_between_markers
from burnman.utils.misc import run_cli_program_with_input
from burnman.utils.math import smooth_array
from burnman.utils.math import compare_l2_norm, l2_norm
from burnman.utils.math import compare_chisqr_factor
from burnman.utils.math import independent_row_indices
from burnman.utils.math import array_to_rational_matrix
from burnman.utils.math import complete_basis, generate_complete_basis


class test_utils(BurnManTest):
    def setUp(self):
        self.array = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        self.grid_spacing = np.array([1.0, 1.0])
        self.gaussian_rms_widths = np.array([1.0, 1.0])

    def test_run_cli(self):
        output = run_cli_program_with_input(["cat"], "Hello, CLI!\n", verbose=False)
        self.assertEqual(output, "Hello, CLI!\n")

    def test_extract_lines(self):
        test_content = """Header
    START_MARK
    Line 1
    Line 2
    END_MARK
    Footer"""

        # Create a temporary file with the test content
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
            tmp.write(test_content)
            tmp_path = tmp.name

        # Exclusive test (default)
        result = extract_lines_between_markers(tmp_path, "START_MARK", "END_MARK")
        self.assertEqual(result[0], "Line 1")
        self.assertEqual(result[1], "Line 2")

        # Inclusive test
        result = extract_lines_between_markers(
            tmp_path, "START_MARK", "END_MARK", inclusive=True
        )
        self.assertEqual(result[0], "START_MARK")
        self.assertEqual(result[1], "Line 1")
        self.assertEqual(result[2], "Line 2")
        self.assertEqual(result[3], "END_MARK")

        os.remove(tmp_path)  # Clean up

    def test_smooth_array_basic_gaussian_smoothing(self):
        # Test standard Gaussian smoothing using 'reflect'
        expected = gaussian_filter(self.array, sigma=(1.0, 1.0), mode="reflect")
        result = smooth_array(
            self.array, self.grid_spacing, self.gaussian_rms_widths, mode="reflect"
        )
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_smooth_array_invalid_mode_uses_scipy_fallback(self):
        # Ensure fallback to scipy behavior for standard modes
        result = smooth_array(
            self.array, self.grid_spacing, self.gaussian_rms_widths, mode="nearest"
        )
        expected = gaussian_filter(self.array, sigma=(1.0, 1.0), mode="nearest")
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_smooth_array_different_grid_spacing(self):
        # Test effect of non-uniform grid spacing
        grid_spacing = np.array([0.5, 2.0])
        sigma = tuple(self.gaussian_rms_widths / grid_spacing)
        expected = gaussian_filter(self.array, sigma=sigma, mode="reflect")
        result = smooth_array(
            self.array, grid_spacing, self.gaussian_rms_widths, mode="reflect"
        )
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_l2_norm_exact_match(self):
        x = np.linspace(0, 10, 100)
        f = np.sin(x)
        result = l2_norm(x, f, f)
        self.assertAlmostEqual(result, 0.0, places=10)

    def test_l2_norm_linear_difference(self):
        x = np.array([0, 1, 2])
        f1 = np.array([1, 2, 3])
        f2 = np.array([1, 1, 1])
        # diff = [0, 1, 2], squared = [0, 1, 4]
        # integrate: ∫0^2 of squared difference = trapezoid([0,1,4], x=[0,1,2]) = 3
        result = l2_norm(x, f1, f2)
        self.assertAlmostEqual(result, np.sqrt(3), places=4)

    def test_l2_norm_constant_vs_zero(self):
        x = np.linspace(0, 1, 5)
        f1 = np.ones_like(x)
        f2 = np.zeros_like(x)
        expected = np.sqrt(integrate.trapezoid(f1**2, x))
        result = l2_norm(x, f1, f2)
        self.assertAlmostEqual(result, expected, places=6)

    def test_l2_norm_nonuniform_grid(self):
        x = np.array([0.0, 0.1, 0.4, 1.0])
        f1 = np.array([1.0, 2.0, 3.0, 4.0])
        f2 = np.array([1.0, 2.0, 1.0, 0.0])
        result = l2_norm(x, f1, f2)
        self.assertTrue(result > 0)

    def test_l2_norm_input_shape_mismatch(self):
        x = np.array([0, 1, 2])
        f1 = np.array([1, 2, 3])
        f2 = np.array([1, 2])
        with self.assertRaises(ValueError):
            l2_norm(x, f1, f2)

    def test_compare_l2_norm_constant_vs_zero_profiles(self):
        depth = np.linspace(0, 1, 10)
        calc = [np.ones_like(depth), 2 * np.ones_like(depth)]
        obs = [np.zeros_like(depth), np.zeros_like(depth)]
        result = compare_l2_norm(depth, calc, obs)
        expected_1 = np.sqrt(integrate.trapezoid(np.ones_like(depth) ** 2, depth))
        expected_2 = np.sqrt(integrate.trapezoid((2 * np.ones_like(depth)) ** 2, depth))
        self.assertAlmostEqual(result[0], expected_1, places=6)
        self.assertAlmostEqual(result[1], expected_2, places=6)

    def test_compare_l2_norm_mismatched_profile_lengths(self):
        depth = np.array([0, 1, 2])
        calc = [np.array([1, 2, 3])]
        obs = [np.array([1, 2])]
        with self.assertRaises(ValueError):
            compare_l2_norm(depth, calc, obs)

    def test_chisqr_factor_perfect_match(self):
        profile1 = np.array([1.0, 2.0, 3.0])
        profile2 = np.array([4.0, 5.0, 6.0])
        calc = [profile1, profile2]
        obs = [profile1.copy(), profile2.copy()]
        result = compare_chisqr_factor(calc, obs)
        for val in result:
            self.assertAlmostEqual(val, 0.0, places=12)

    def test_chisqr_factor_known_difference(self):
        # Calculate chi-sqr manually for a simple case
        obs = np.array([10, 10, 10])
        calc = np.array([10, 11, 12])
        # chi-squared terms: ((calc - obs)/(0.01*mean(obs)))**2
        # mean(obs) = 10, 1% = 0.1
        # terms = [0, (1/0.1)^2=100, (2/0.1)^2=400] -> sum=500, avg=500/3=166.6667
        expected = 166.6667

        result = compare_chisqr_factor([calc], [obs])
        self.assertAlmostEqual(result[0], expected, places=3)

    def test_chisqr_factor_multiple_profiles(self):
        obs1 = np.array([1, 2, 3])
        obs2 = np.array([4, 5, 6])
        calc1 = np.array([1.1, 2.1, 3.1])
        calc2 = np.array([3.9, 4.9, 6.1])
        result = compare_chisqr_factor([calc1, calc2], [obs1, obs2])

        # Compute expected manually
        def manual_chi(calc, obs):
            diff = (calc - obs) / (0.01 * np.mean(obs))
            return np.mean(diff**2)

        expected = [manual_chi(calc1, obs1), manual_chi(calc2, obs2)]
        for r, e in zip(result, expected):
            self.assertAlmostEqual(r, e, places=6)

    def test_chisqr_factor_empty_input(self):
        # Empty input returns empty list
        result = compare_chisqr_factor([], [])
        self.assertEqual(result, [])

    def test_chisqr_factor_mismatched_lengths(self):
        # If calc and obs length mismatch, should raise IndexError or handle gracefully
        calc = [np.array([1, 2])]
        obs = []
        with self.assertRaises(IndexError):
            compare_chisqr_factor(calc, obs)

    def test_independent_row_indices(self):
        arr = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # All rows are independent
        result = independent_row_indices(arr)
        self.assertEqual(result, [0, 1, 2])

    def test_independent_row_indices_dependent_rows(self):
        arr = np.array([[1, 2], [2, 4], [3, 5]])  # This row is linearly dependent
        result = independent_row_indices(arr)
        # Expect two independent rows (rank 2)
        self.assertEqual(len(result), 2)
        self.assertTrue(0 in result)
        self.assertTrue(2 in result)

    def test_independent_row_indices_swap_logic(self):
        arr = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        # Row reduction should cause at least one row swap
        result = independent_row_indices(arr)
        self.assertEqual(len(result), 3)
        self.assertEqual(sorted(result), [0, 1, 2])

    def test_independent_row_indices_rank_deficient_matrix(self):
        arr = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3]])
        result = independent_row_indices(arr)
        self.assertEqual(result, [0])  # Only the first row is independent

    def test_independent_row_indices_floating_point_precision(self):
        arr = np.array([[1.0, 2.0], [2.0, 4.000001], [0.0, 1.0]])  # Nearly dependent
        result = independent_row_indices(arr)
        self.assertEqual(len(result), 2)

    def test_array_to_rational_matrix(self):
        arr = np.array([[0.5, 0.333], [0.25, 1.125]])
        result = array_to_rational_matrix(arr)

        self.assertIsInstance(result, Matrix)
        self.assertEqual(result.shape, (2, 2))

        # Check entries are Rational and values are close
        expected = [
            [Rational(1, 2), Rational(333, 1000)],
            [Rational(1, 4), Rational(9, 8)],
        ]

        for i in range(2):
            for j in range(2):
                self.assertEqual(result[i, j], expected[i][j])
                self.assertIsInstance(result[i, j], Rational)

    def test_complete_basis_flip(self):
        basis = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 0, 0, 0]])
        incomplete_basis = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
        basis2 = complete_basis(incomplete_basis, True)
        np.testing.assert_allclose(basis, basis2, rtol=1e-5)

    def test_complete_basis_no_flip(self):
        basis = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1]])
        incomplete_basis = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
        basis2 = complete_basis(incomplete_basis, False)
        np.testing.assert_allclose(basis, basis2, rtol=1e-5)

    def test_complete_basis_full(self):
        basis = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        basis2 = complete_basis(basis)
        np.testing.assert_allclose(basis, basis2, rtol=1e-5)

    def test_generate_complete_basis(self):
        incomplete_basis = np.array([[1, 0, 0], [1, 1, 0]])
        array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        c = np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]])
        c2 = generate_complete_basis(incomplete_basis, array)
        np.testing.assert_allclose(c, c2, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
