import numpy as np
import pandas as pd
import pytest

from kineticanalysis.analysis import analysis_track


# ---------------------------------------------------------------------------
# check_continuous_time
# ---------------------------------------------------------------------------
class TestCheckContinuousTime:
    def test_evenly_spaced_points_are_continuous(self):
        x = np.arange(0, 10) * 0.5
        assert analysis_track.check_continuous_time(x, dt=0.5) is True

    def test_missing_point_is_not_continuous(self):
        x = np.array([0, 0.5, 1.0, 2.5, 3.0])
        assert analysis_track.check_continuous_time(x, dt=0.5) is False

    def test_small_jitter_within_rtol_is_continuous(self):
        x = np.array([0, 0.5001, 0.9999, 1.5002])
        assert analysis_track.check_continuous_time(x, dt=0.5, rtol=1e-2) is True


# ---------------------------------------------------------------------------
# check_track_validity
# ---------------------------------------------------------------------------
class TestCheckTrackValidity:
    def test_continuous_track_is_valid_and_unchanged_length(self):
        df = pd.DataFrame({
            "TRACK_ID": [1] * 6,
            "FRAME": [0, 1, 2, 3, 4, 5],
            "MEAN_INTENSITY_CH1": [10, 12, 11, 13, 14, 15],
        })
        valid, x_orig, y_orig, x_fixed, y_fixed = analysis_track.check_track_validity(
            df, id_track=1, delta_t=0.5
        )
        assert valid is True
        assert len(x_fixed) == len(y_fixed) == 6

    def test_small_gap_is_interpolated(self):
        # FRAME=3 is missing; gap is small enough to be filled in.
        df = pd.DataFrame({
            "TRACK_ID": [1] * 5,
            "FRAME": [0, 1, 2, 4, 5],
            "MEAN_INTENSITY_CH1": [10, 12, 11, 14, 15],
        })
        valid, x_orig, y_orig, x_fixed, y_fixed = analysis_track.check_track_validity(
            df, id_track=1, delta_t=0.5
        )
        assert valid is True
        assert len(x_orig) == 5          # original, unfixed data
        assert len(x_fixed) == 6         # one point interpolated in
        # interpolated intensity should sit between its two neighbours
        assert 11 < y_fixed[3] < 14

    def test_large_gap_is_rejected(self):
        df = pd.DataFrame({
            "TRACK_ID": [1] * 4,
            "FRAME": [0, 1, 2, 50],
            "MEAN_INTENSITY_CH1": [10, 12, 11, 14],
        })
        valid, *_ = analysis_track.check_track_validity(df, id_track=1, delta_t=0.5)
        assert valid is False

    def test_only_requested_track_id_is_used(self):
        df = pd.DataFrame({
            "TRACK_ID": [1, 1, 1, 2, 2, 2],
            "FRAME": [0, 1, 2, 0, 1, 2],
            "MEAN_INTENSITY_CH1": [10, 12, 11, 100, 100, 100],
        })
        _, _, y_orig, _, _ = analysis_track.check_track_validity(
            df, id_track=1, delta_t=0.5
        )
        assert list(y_orig) == [10, 12, 11]


# ---------------------------------------------------------------------------
# recommend_acquisition
# ---------------------------------------------------------------------------
class TestRecommendAcquisition:
    def test_known_values(self):
        result = analysis_track.recommend_acquisition(
            k_est=5.0, protein_length=1500, suntag_length=800,
            nb_full_prot=7, samples_per_ramp=24,
        )
        assert result["tau_c"] == pytest.approx((1500 + 800) / 5.0)
        assert result["tau_ramp"] == pytest.approx(800 / 5.0)
        assert result["T_recommended"] == pytest.approx(7 * result["tau_c"])
        assert result["dt_recommended"] == pytest.approx(result["tau_ramp"] / 24)
        assert result["dt_nyquist_limit"] == pytest.approx(result["tau_ramp"] / 2)

    def test_n_points_is_consistent_with_recommended_window(self):
        result = analysis_track.recommend_acquisition(
            k_est=5.0, protein_length=1500, suntag_length=800,
        )
        expected_n_points = int(result["T_recommended"] / result["dt_recommended"])
        assert result["n_points"] == expected_n_points

    def test_faster_elongation_rate_shortens_recommended_window(self):
        slow = analysis_track.recommend_acquisition(
            k_est=2.0, protein_length=1500, suntag_length=800
        )
        fast = analysis_track.recommend_acquisition(
            k_est=10.0, protein_length=1500, suntag_length=800
        )
        assert fast["T_recommended"] < slow["T_recommended"]


# ---------------------------------------------------------------------------
# correct_elongation_rate / correct_initiation_rate
#
# These are defined twice in the codebase (identically, so far) - once in
# fit_functions.py and once in analysis_track.py. These tests pin the
# behavior of *this* copy so the two implementations can't silently drift
# apart; see tests/analysis/test_fit_functions.py for the other copy.
# ---------------------------------------------------------------------------
class TestQueuingCorrectionInAnalysisTrack:
    @pytest.mark.parametrize(
        "fn", [analysis_track.correct_elongation_rate, analysis_track.correct_initiation_rate]
    )
    def test_no_correction_when_density_is_zero(self, fn):
        assert fn(2.0, 0) == pytest.approx(2.0)

    @pytest.mark.parametrize(
        "fn", [analysis_track.correct_elongation_rate, analysis_track.correct_initiation_rate]
    )
    def test_returns_nan_for_invalid_density(self, fn):
        assert np.isnan(fn(2.0, None))
        assert np.isnan(fn(2.0, 1))
        assert np.isnan(fn(2.0, 1.5))


# ---------------------------------------------------------------------------
# single_track_analysis
#
# autocorrelation() and the lmfit-based fit functions are monkeypatched so
# these tests exercise the dispatch/branching logic in
# single_track_analysis itself without depending on multipletau/lmfit
# actually converging on real data.
# ---------------------------------------------------------------------------
class TestSingleTrackAnalysis:
    @pytest.fixture(autouse=True)
    def _stub_autocorrelation(self, monkeypatch):
        x_auto = np.arange(10, dtype=float)
        y_auto = np.linspace(1.0, 0.1, 10)
        monkeypatch.setattr(
            analysis_track, "autocorrelation", lambda y, delta_t, normalize, mm: (x_auto, y_auto)
        )
        monkeypatch.setattr(
            analysis_track, "fit_autocorrelation_approx",
            lambda x, y, method="lm": (2.0, 0.5, [0.1, 0.05]),
        )

    def test_unknown_method_returns_all_nan(self):
        x = np.arange(10, dtype=float)
        y = np.linspace(1.0, 0.1, 10)
        result = analysis_track.single_track_analysis(x, y, method="not-a-real-method")
        _, _, k, c, elongation_r, translation_init_r, perr = result
        assert np.isnan(k)
        assert np.isnan(c)
        assert np.isnan(elongation_r)
        assert np.isnan(translation_init_r)
        assert perr == [np.nan, np.nan]

    def test_approx_method_uses_stubbed_fit(self):
        x = np.arange(10, dtype=float)
        y = np.linspace(1.0, 0.1, 10)
        _, _, k, c, elongation_r, translation_init_r, _ = analysis_track.single_track_analysis(
            x, y, method="approx", protein_size=1500, suntag_size=800, repetition_suntag=32,
        )
        assert k == 2.0
        assert c == 0.5
        assert translation_init_r == c
        # elongation_r = M / k * one_suntag_size, with the stubbed k
        one_suntag_size = 800 // 32
        M = 1500 // one_suntag_size
        assert elongation_r == pytest.approx(M / k * one_suntag_size)

    def test_correct_queuing_requires_mean_n_ribosome(self):
        x = np.arange(10, dtype=float)
        y = np.linspace(1.0, 0.1, 10)
        with pytest.raises(ValueError, match="mean_n_ribosome"):
            analysis_track.single_track_analysis(
                x, y, method="approx", correct_queuing=True, mean_n_ribosome=None,
            )

    def test_correct_queuing_rescales_rates(self):
        x = np.arange(10, dtype=float)
        y = np.linspace(1.0, 0.1, 10)

        uncorrected = analysis_track.single_track_analysis(x, y, method="approx")
        corrected = analysis_track.single_track_analysis(
            x, y, method="approx", correct_queuing=True, mean_n_ribosome=5,
            protein_size=1500, suntag_size=800,
        )

        _, _, _, _, elong_u, init_u, _ = uncorrected
        _, _, _, _, elong_c, init_c, _ = corrected

        # Correction divides by (1 - rho_bar) with rho_bar > 0, so the
        # corrected rates must be strictly larger than the uncorrected ones.
        assert elong_c > elong_u
        assert init_c > init_u