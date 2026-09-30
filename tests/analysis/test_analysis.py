import numpy as np
import pytest

from kineticanalysis.analysis.fit_functions import (
    correct_elongation_rate,
    correct_initiation_rate,
    fit_function_linear,
    function_approx,
    function_epitope,
    function_exact,
)


# ---------------------------------------------------------------------------
# function_approx: piecewise-linear autocorrelation model
#   ((t - x) / (c * t**2)) * heaviside(t - x, 0)
# ---------------------------------------------------------------------------
class TestFunctionApprox:
    def test_value_at_x_zero(self):
        # At x=0 the function reduces to 1/(c*t)
        t, c = 10.0, 0.5
        assert function_approx(0, t, c) == pytest.approx(1 / (c * t))

    def test_value_decreases_linearly_before_t(self):
        t, c = 10.0, 0.5
        v0 = function_approx(0, t, c)
        v5 = function_approx(5, t, c)
        # halfway to t, value should have halved (linear ramp to zero)
        assert v5 == pytest.approx(v0 / 2)

    def test_zero_at_and_after_t(self):
        t, c = 10.0, 0.5
        assert function_approx(10, t, c) == 0
        assert function_approx(15, t, c) == 0

    def test_vectorized_input(self):
        t, c = 10.0, 0.5
        x = np.array([0, 5, 10, 20])
        out = function_approx(x, t, c)
        np.testing.assert_allclose(out, [0.2, 0.1, 0.0, 0.0])


# ---------------------------------------------------------------------------
# correct_elongation_rate / correct_initiation_rate
# (same formula, duplicated between fit_functions.py and analysis_track.py)
# ---------------------------------------------------------------------------
class TestQueuingCorrection:
    @pytest.mark.parametrize("fn", [correct_elongation_rate, correct_initiation_rate])
    def test_no_correction_when_density_is_zero(self, fn):
        assert fn(2.0, 0) == pytest.approx(2.0)

    @pytest.mark.parametrize("fn", [correct_elongation_rate, correct_initiation_rate])
    def test_scales_up_with_density(self, fn):
        # rho_bar = 0.5 -> value should double
        assert fn(2.0, 0.5) == pytest.approx(4.0)

    @pytest.mark.parametrize("fn", [correct_elongation_rate, correct_initiation_rate])
    def test_returns_nan_when_density_is_none(self, fn):
        assert np.isnan(fn(2.0, None))

    @pytest.mark.parametrize("fn", [correct_elongation_rate, correct_initiation_rate])
    @pytest.mark.parametrize("rho_bar", [1, 1.2])
    def test_returns_nan_when_density_at_or_above_one(self, fn, rho_bar):
        # rho_bar >= 1 is unphysical (fully occupied or over-occupied mRNA)
        # and would otherwise divide by zero or flip sign.
        assert np.isnan(fn(2.0, rho_bar))


# ---------------------------------------------------------------------------
# function_exact / function_epitope: autocorrelation models used in the fit
# Kept to small N/M so the factorial-based series stays fast and exact.
# ---------------------------------------------------------------------------
class TestExactAndEpitopeFunctions:
    @pytest.mark.parametrize("func", [function_exact, function_epitope])
    def test_output_is_finite(self, func):
        k, c, N, M = 0.6, 0.1, 4, 8
        args = (0.0, k, c, N, M) if func is function_exact else (0.0, k, c, N)
        out = float(func(*args))
        assert np.isfinite(out)

    @pytest.mark.parametrize("func", [function_exact, function_epitope])
    def test_decays_monotonically_with_x(self, func):
        # An autocorrelation curve should decay as the lag x grows.
        k, c, N, M = 0.6, 0.1, 4, 8
        xs = [0, 1, 2, 5, 10]
        if func is function_exact:
            values = [float(func(x, k, c, N, M)) for x in xs]
        else:
            values = [float(func(x, k, c, N)) for x in xs]
        assert values == sorted(values, reverse=True)

    def test_exact_matches_epitope_regression_values(self):
        # Regression guard: locks in the current numeric output for a fixed
        # set of parameters so future refactors of the factorial series
        # (e.g. removing np.float128) don't silently change results.
        k, c, N, M = 0.6, 0.1, 4, 8
        assert float(function_exact(1.0, k, c, N, M)) == pytest.approx(
            0.5152923145492805, rel=1e-6
        )
        assert float(function_epitope(1.0, k, c, N)) == pytest.approx(
            1.4529458778607691, rel=1e-6
        )


# ---------------------------------------------------------------------------
# fit_function_linear
# ---------------------------------------------------------------------------
class TestFitFunctionLinear:
    def test_recovers_known_slope_and_intercept(self):
        rng = np.random.default_rng(0)
        x = np.arange(0, 20, dtype=float)
        true_slope, true_intercept = -0.5, 5.0
        y = true_slope * x + true_intercept + rng.normal(0, 0.01, size=x.size)

        slope, intercept, _ = fit_function_linear(x, y)

        assert slope == pytest.approx(true_slope, abs=0.05)
        assert intercept == pytest.approx(true_intercept, abs=0.1)

    def test_too_few_points_returns_sentinel(self):
        # If the sign change happens in the first point, fewer than 2
        # points are available to fit and the function should bail out
        # with its documented sentinel values rather than raising.
        x = np.array([0.0, 1.0])
        y = np.array([-1.0, -2.0])
        assert fit_function_linear(x, y) == (-1, -1, [-1, -1])

    def test_all_positive_curve_raises_indexerror(self):
        # Known limitation: the function assumes the autocorrelation curve
        # eventually crosses zero. If it never does (e.g. a noisy but
        # always-positive curve), `np.where(y_sign_value == -1)[0][0]`
        # indexes into an empty array and raises IndexError instead of
        # failing gracefully. This test documents the current behavior so
        # a future fix (returning a sentinel instead of crashing) is a
        # deliberate, visible change rather than an accidental one.
        x = np.arange(0, 20, dtype=float)
        y = np.abs(-0.5 * x + 5.0) + 1.0  # always positive
        with pytest.raises(IndexError):
            fit_function_linear(x, y)