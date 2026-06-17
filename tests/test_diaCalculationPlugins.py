# This file is part of ap_association.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import warnings

from astropy.stats import median_absolute_deviation
import numpy as np
import pandas as pd
from scipy.stats import skew
import unittest

from lsst.meas.base import (
    MeanDiaPosition, MeanDiaPositionConfig,
    HTMIndexDiaPosition, HTMIndexDiaPositionConfig,
    NumDiaSourcesDiaPlugin, NumDiaSourcesDiaPluginConfig,
    SimpleSourceFlagDiaPlugin, SimpleSourceFlagDiaPluginConfig,
    WeightedMeanDiaPsfFlux, WeightedMeanDiaPsfFluxConfig,
    PercentileDiaPsfFlux, PercentileDiaPsfFluxConfig,
    SigmaDiaPsfFlux, SigmaDiaPsfFluxConfig,
    Chi2DiaPsfFlux, Chi2DiaPsfFluxConfig,
    MadDiaPsfFlux, MadDiaPsfFluxConfig,
    SkewDiaPsfFlux, SkewDiaPsfFluxConfig,
    MinMaxDiaPsfFlux, MinMaxDiaPsfFluxConfig,
    MaxSlopeDiaPsfFlux, MaxSlopeDiaPsfFluxConfig,
    ErrMeanDiaPsfFlux, ErrMeanDiaPsfFluxConfig,
    LinearFitDiaPsfFlux, LinearFitDiaPsfFluxConfig,
    StetsonJDiaPsfFlux, StetsonJDiaPsfFluxConfig,
    WeightedMeanDiaTotFlux, WeightedMeanDiaTotFluxConfig,
    SigmaDiaTotFlux, SigmaDiaTotFluxConfig,
    LombScarglePeriodogram, LombScarglePeriodogramConfig,
    LombScarglePeriodogramMulti, LombScarglePeriodogramMultiConfig,
    UnphysicalDiaSourceSeparation)
import lsst.utils.tests


def run_single_plugin(diaObjectCat,
                      diaObjectId,
                      diaSourceCat,
                      band,
                      plugin):
    """Wrapper for running single plugins.

    Reproduces some of the behavior of `lsst.ap.association.DiaCalcuation.run`

    Parameters
    ----------
    diaObjectCat : `pandas.DataFrame`
        Input object catalog to store data into and read from.
    diaSourcesCat : `pandas.DataFrame`
        DiaSource catalog to read data from and groupby on.
    fitlerName : `str`
        String name of the filter to process.
    plugin : `lsst.ap.association.DiaCalculationPlugin`
        Plugin to run.
    """
    diaObjectCat.set_index("diaObjectId", inplace=True, drop=False)
    diaSourceCat.set_index(
        ["diaObjectId", "band", "diaSourceId"],
        inplace=True,
        drop=False)

    objDiaSources = diaSourceCat.loc[diaObjectId]
    updatingFilterDiaSources = diaSourceCat.loc[
        (diaObjectId, band), :
    ]

    plugin.calculate(diaObjects=diaObjectCat,
                     diaObjectId=diaObjectId,
                     diaSources=objDiaSources,
                     filterDiaSources=updatingFilterDiaSources,
                     band=band)


def run_multi_plugin(diaObjectCat, diaSourceCat, band, plugin):
    """Wrapper for running multi plugins.

    Reproduces some of the behavior of `lsst.ap.association.DiaCalcuation.run`

    Parameters
    ----------
    diaObjectCat : `pandas.DataFrame`
        Input object catalog to store data into and read from.
    diaSourcesCat : `pandas.DataFrame`
        DiaSource catalog to read data from and groupby on.
    filterName : `str`
        String name of the filter to process.
    plugin : `lsst.ap.association.DiaCalculationPlugin`
        Plugin to run.
    """
    diaObjectCat.set_index("diaObjectId", inplace=True, drop=False)
    diaSourceCat.set_index(
        ["diaObjectId", "band", "diaSourceId"],
        inplace=True,
        drop=False)

    updatingFilterDiaSources = diaSourceCat.loc[
        (slice(None), band), :
    ]

    diaSourcesGB = diaSourceCat.groupby(level=0)
    filterDiaSourcesGB = updatingFilterDiaSources.groupby(level=0)

    plugin.calculate(diaObjects=diaObjectCat,
                     diaSources=diaSourcesGB,
                     filterDiaSources=filterDiaSourcesGB,
                     band=band)


def run_multiband_plugin(diaObjectCat, diaSourceCat, plugin):
    """Wrapper for running multi plugins.

    Reproduces some of the behavior of `lsst.ap.association.DiaCalcuation.run`

    Parameters
    ----------
    diaObjectCat : `pandas.DataFrame`
        Input object catalog to store data into and read from.
    diaSourcesCat : `pandas.DataFrame`
        DiaSource catalog to read data from and groupby on.
    plugin : `lsst.ap.association.DiaCalculationPlugin`
        Plugin to run.
    """
    diaObjectCat.set_index("diaObjectId", inplace=True, drop=False)
    diaSourceCat.set_index(
        ["diaObjectId", "band", "diaSourceId"],
        inplace=True,
        drop=False)

    diaSourcesGB = diaSourceCat.groupby(level=0)

    plugin.calculate(diaObjects=diaObjectCat,
                     diaSources=diaSourcesGB,
                     )


def make_diaObject_table(objId, plugin, default_value=np.nan, band=None):
    """Create a minimal diaObject table with columns required for the plugin

    Parameters
    ----------
    objId : `int`
        The diaObjectId
    plugin : `lsst.ap.association.DiaCalculationPlugin`
        The plugin that will be run.
    default_value : `float` or `int`, optional
        Value to set new columns to.
    band : `str`, optional
        Band designation to append to the plugin columns.

    Returns
    -------
    diaObjects : `pandas.DataFrame`
        Output catalog with the required columns for the plugin.
    """
    # Add an extra empty diaObject here. This ensures that
    # we properly test the source/object matching implicit
    # in the plugin calculations.
    diaObjects = {"diaObjectId": [objId, objId + 1]}
    for col in plugin.outputCols:
        if band is not None:
            diaObjects[f"{band}_{col}"] = default_value
        else:
            diaObjects[col] = default_value
    return pd.DataFrame(diaObjects)


class TestMeanPosition(unittest.TestCase):

    def testCalculate(self):
        """Test mean position calculation.

        DiaSources here are constructed without per-source coordinate
        errors, so each ``run_multi_plugin`` call legitimately triggers
        the no-errors warning from the plugin.  Suppress it here so the
        signal of an unexpected warning elsewhere is not drowned out.
        """
        n_sources = 10
        objId = 0

        # configure a 2 degree max separation
        plug = MeanDiaPosition(MeanDiaPositionConfig(MaxAllowedDiaSourceSeparation=7200.0),
                               "ap_meanPosition",
                               None)

        warnings.filterwarnings("ignore",
                                message="No DiaSources with finite coordinate errors",
                                category=UserWarning)
        self.addCleanup(warnings.resetwarnings)

        # Test expected means in RA.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1, 1, n_sources),
                                        "dec": np.zeros(n_sources),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "dec"], 0.0)

        # Test expected means in DEC.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "dec": np.linspace(-1, 1, n_sources),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "dec"], 0.0)

        # Test failure mode RA is nan.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.full(n_sources, np.nan),
                                        "dec": np.zeros(n_sources),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "dec"]))

        # Test failure mode DEC is nan.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "dec": np.full(n_sources, np.nan),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "dec"]))

        # configure the default 3 arcsecond separation
        plug = MeanDiaPosition(MeanDiaPositionConfig(MaxAllowedDiaSourceSeparation=3.0),
                               "ap_meanPosition",
                               None)

        # These 1 degree separations should raise
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1, 1, n_sources),
                                        "dec": np.zeros(n_sources),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        with self.assertRaises(UnphysicalDiaSourceSeparation):
            run_multi_plugin(diaObjects, diaSources, "g", plug)

        # 1 arcsecond separations should not raise
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1/3600., 1/3600., n_sources),
                                        "dec": np.zeros(n_sources),
                                        "midpointMjdTai": np.linspace(0, n_sources, n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "band": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

    def _makeDiaSourcesWithUncertainties(self, raErr, decErr, ra_dec_Cov, ra=None, dec=None, objId=0):
        """Build a tiny DiaSource DataFrame with per-source uncertainties.
        """
        n = len(raErr)
        if ra is None:
            ra = np.zeros(n)
        if dec is None:
            dec = np.zeros(n)
        return pd.DataFrame(data={
            "ra": ra,
            "dec": dec,
            "raErr": raErr,
            "decErr": decErr,
            "ra_dec_Cov": ra_dec_Cov,
            "midpointMjdTai": np.arange(n, dtype=float),
            "diaObjectId": n * [objId],
            "band": n * ["g"],
            "diaSourceId": np.arange(n, dtype=int),
        })

    def testUncertaintyDiagonalOnly(self):
        """Two coincident sources, no off-diagonal: chi^2 = 0, scale
        factor is 1, so the output is the diagonal inverse-variance
        weighted-mean covariance.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        raErr = np.array([1e-6, 2e-6])
        decErr = np.array([1e-6, 2e-6])
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(
            raErr=raErr,
            decErr=decErr,
            ra_dec_Cov=np.array([np.nan, np.nan]),
        )
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        # Var(weighted mean) = 1 / sum(1/sigma^2) per axis.
        expectedRaErr = 1.0/np.sqrt(np.sum(1.0/raErr**2))
        expectedDecErr = 1.0/np.sqrt(np.sum(1.0/decErr**2))
        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], expectedRaErr)
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], expectedDecErr)
        # No per-source ra_dec_Cov, so output ra_dec_Cov is NaN.
        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra_dec_Cov"]))

    def testUncertaintyFullCovariance(self):
        """Two coincident sources with off-diagonal covariance: chi^2 = 0,
        so the output is C_formal = (sum_i inv(C_i))^-1.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        raErr = np.array([1e-6, 2e-6])
        decErr = np.array([1.5e-6, 1.0e-6])
        rho = np.array([0.3, -0.2])  # correlation coefficient
        raDecCov = rho * raErr * decErr

        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(
            raErr=raErr, decErr=decErr, ra_dec_Cov=raDecCov)
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        n = len(raErr)
        cov = np.zeros((n, 2, 2))
        cov[:, 0, 0] = raErr**2
        cov[:, 1, 1] = decErr**2
        cov[:, 0, 1] = raDecCov
        cov[:, 1, 0] = raDecCov
        covObj = np.linalg.inv(np.linalg.inv(cov).sum(axis=0))

        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], np.sqrt(covObj[0, 0]))
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], np.sqrt(covObj[1, 1]))
        self.assertAlmostEqual(diaObjects.loc[objId, "ra_dec_Cov"], covObj[0, 1])

    def testUncertaintySingleSourceCopiesThrough(self):
        """A single-source group: outputs equal the source's uncertainties.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(
            raErr=np.array([1.5e-6]),
            decErr=np.array([0.8e-6]),
            ra_dec_Cov=np.array([2e-13]),
        )
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], 1.5e-6)
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], 0.8e-6)
        self.assertAlmostEqual(diaObjects.loc[objId, "ra_dec_Cov"], 2e-13)

    def testUncertaintyMissingColumnsFallsBackToScatter(self):
        """No raErr/decErr columns + N>=2 spread-out sources -> the
        weighted-mean path is unusable, so the DiaObject falls back to
        the unweighted mean position with the scatter-only SEM, and a
        warning is emitted.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        n = 3
        ra = np.linspace(-1e-4, 1e-4, n)
        dec = np.zeros(n)
        diaSources = pd.DataFrame(data={
            "ra": ra,
            "dec": dec,
            "midpointMjdTai": np.arange(n, dtype=float),
            "diaObjectId": n * [objId],
            "band": n * ["g"],
            "diaSourceId": np.arange(n, dtype=int),
        })
        with self.assertWarnsRegex(UserWarning, "No DiaSources with finite coordinate errors"):
            run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "dec"], 0.0)

        # Expected scatter term: at dec=0, the tangent-plane east offsets
        # equal the RA deltas (in degrees) to high precision, so the
        # standard error of the mean in RA is sample_std(ra, ddof=1) /
        # sqrt(N).  Dec is identically zero so its scatter is zero.
        expectedRaErr = np.std(ra, ddof=1)/np.sqrt(n)
        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], expectedRaErr)
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "ra_dec_Cov"], 0.0)

    def testUncertaintyMissingColumnsSingleSourceEmitsNaN(self):
        """With a single source and no raErr/decErr the uncertainty is
        undefined but the mean position is still computed; a warning is
        emitted that the per-source errors are unusable.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={
            "ra": [0.0],
            "dec": [0.0],
            "midpointMjdTai": [0.0],
            "diaObjectId": [objId],
            "band": ["g"],
            "diaSourceId": [0],
        })
        with self.assertWarnsRegex(UserWarning, "No DiaSources with finite coordinate errors"):
            run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "dec"], 0.0)
        self.assertTrue(np.isnan(diaObjects.loc[objId, "raErr"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "decErr"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra_dec_Cov"]))

    def testUncertaintyPartialCovarianceFallsBackToDiagonal(self):
        """Some sources have NaN ra_dec_Cov but valid raErr/decErr.

        When at least one included source lacks a finite ra_dec_Cov,
        diagonal weights are used for all included sources, the off-
        diagonal output is NaN, and the chi-squared falls back to its
        diagonal form.  Sources are coincident here, so chi^2 = 0 and
        there is no rescaling.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        raErr = np.array([1e-6, 2e-6, 3e-6])
        decErr = np.array([1e-6, 2e-6, 3e-6])
        # Only one source has a finite covariance
        raDecCov = np.array([1e-13, np.nan, np.nan])
        diaSources = self._makeDiaSourcesWithUncertainties(raErr=raErr, decErr=decErr, ra_dec_Cov=raDecCov)
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        # Diagonal inverse-variance weighted mean covariance.
        expectedRaErr = 1.0/np.sqrt(np.sum(1.0/raErr**2))
        expectedDecErr = 1.0/np.sqrt(np.sum(1.0/decErr**2))
        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], expectedRaErr)
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], expectedDecErr)
        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra_dec_Cov"]))

    def testUncertaintyScaleFactorInflatesWhenScatterExceedsErrors(self):
        """Two sources whose positional separation is much larger than
        their per-source errors: chi^2 >> dof, so the chi-squared scale
        factor inflates the formal covariance by chi^2 / dof.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        # Two sources symmetric around RA = 0, Dec = 0.  Spread of
        # 0.36 arcsec on each axis (within MaxAllowedDiaSourceSeparation
        # of 3 arcsec), still 100x the per-source sigma of 1e-6 deg = 3.6 mas,
        # so chi^2 / dof ~ 4e4 / 2 ~ 2e4.
        delta = 1e-4  # degrees
        ra = np.array([-delta, delta])
        dec = np.array([-delta, delta])
        raErr = np.array([1e-6, 1e-6])
        decErr = np.array([1e-6, 1e-6])
        raDecCov = np.array([np.nan, np.nan])

        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(
            raErr=raErr, decErr=decErr, ra_dec_Cov=raDecCov, ra=ra, dec=dec)
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        # With equal per-source errors the weighted mean equals the
        # unweighted mean = (0, 0), so the tangent-plane residuals are
        # essentially (ra, dec).
        # chi^2 = sum_i (r_east^2/sigma_ra^2 + r_north^2/sigma_dec^2)
        # dof = 2 * (N - 1) = 2.
        chi2 = np.sum(ra**2/raErr**2 + dec**2/decErr**2)
        dof = 2*(len(raErr) - 1)
        scale = max(1.0, chi2/dof)
        # Diagonal weighted-mean variance: 1/sum(1/sigma^2) per axis.
        varRaFormal = 1.0/np.sum(1.0/raErr**2)
        varDecFormal = 1.0/np.sum(1.0/decErr**2)
        expectedRaErr = np.sqrt(scale*varRaFormal)
        expectedDecErr = np.sqrt(scale*varDecFormal)

        # The inflation should be large (scale ~ 1e6):
        self.assertGreater(scale, 1e3)

        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"] / expectedRaErr, 1.0, places=5)
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"] / expectedDecErr, 1.0, places=5)
        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra_dec_Cov"]))

    def testUncertaintyScaleFactorNoInflationWhenConsistent(self):
        """When per-source positions are coincident, chi^2 = 0 and the
        chi-squared scale factor is exactly 1: the output equals
        C_formal = (sum_i inv(C_i))^-1.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        raErr = np.array([1e-6, 2e-6, 1.5e-6])
        decErr = np.array([1e-6, 1.5e-6, 2e-6])
        raDecCov = np.array([2e-13, -1e-13, 3e-14])

        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(raErr=raErr, decErr=decErr, ra_dec_Cov=raDecCov)
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        n = len(raErr)
        cov = np.zeros((n, 2, 2))
        cov[:, 0, 0] = raErr**2
        cov[:, 1, 1] = decErr**2
        cov[:, 0, 1] = raDecCov
        cov[:, 1, 0] = raDecCov
        covObj = np.linalg.inv(np.linalg.inv(cov).sum(axis=0))

        self.assertAlmostEqual(diaObjects.loc[objId, "raErr"], np.sqrt(covObj[0, 0]))
        self.assertAlmostEqual(diaObjects.loc[objId, "decErr"], np.sqrt(covObj[1, 1]))
        self.assertAlmostEqual(diaObjects.loc[objId, "ra_dec_Cov"], covObj[0, 1])

    def testWeightedMeanPositionDiffersFromUnweighted(self):
        """Two sources at asymmetric (ra) positions with very unequal
        per-source errors: the reported position is the inverse-variance
        weighted mean, which is much closer to the more-precise source
        than to the unweighted midpoint.
        """
        plug = MeanDiaPosition(MeanDiaPositionConfig(), "ap_meanPosition", None)
        objId = 0
        d = 1e-5  # degrees; well inside MaxAllowedDiaSourceSeparation.
        ra = np.array([-d, d])
        dec = np.array([0.0, 0.0])
        # Source 0 is 100x more precise than source 1.
        raErr = np.array([1e-6, 1e-4])
        decErr = np.array([1e-6, 1e-4])
        raDecCov = np.array([np.nan, np.nan])

        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = self._makeDiaSourcesWithUncertainties(
            raErr=raErr, decErr=decErr, ra_dec_Cov=raDecCov, ra=ra, dec=dec)
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        # Diagonal weighted mean of (ra_0, ra_1) with weights w_i = 1/raErr_i^2.
        w = 1.0/raErr**2
        expectedRa = np.sum(w*ra)/np.sum(w)
        # Output RA is wrapped to [0, 360); normalize back to a signed
        # offset near zero before comparing.
        outRa = ((diaObjects.loc[objId, "ra"] + 180.0) % 360.0) - 180.0
        # The unweighted mean would be 0; the weighted mean should be
        # very close to source 0 at ra = -d.
        self.assertAlmostEqual(outRa, expectedRa)
        self.assertLess(abs(outRa - (-d)), 0.01*d)
        self.assertAlmostEqual(diaObjects.loc[objId, "dec"], 0.0)


class TestHTMIndexPosition(unittest.TestCase):

    def testCalculate(self):
        """Test HTMPixel assignment calculation.
        """
        # Test expected pixelId at RA, DEC = 0
        objId = 0
        n_sources = 10
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaObjects.loc[objId, "ra"] = 0.
        diaObjects.loc[objId, "dec"] = 0.
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int)})
        plug = HTMIndexDiaPosition(HTMIndexDiaPositionConfig(),
                                   "ap_HTMIndex",
                                   None)

        run_single_plugin(diaObjectCat=diaObjects,
                          diaObjectId=objId,
                          diaSourceCat=diaSources,
                          band="g",
                          plugin=plug)
        self.assertEqual(diaObjects.at[objId, "pixelId"],
                         17042430230528)

        # Test expected pixelId at some value of RA and DEC.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaObjects.loc[objId, "ra"] = 45.37
        diaObjects.loc[objId, "dec"] = 13.67
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int)})
        run_single_plugin(diaObjectCat=diaObjects,
                          diaObjectId=objId,
                          diaSourceCat=diaSources,
                          band="g",
                          plugin=plug)
        self.assertEqual(diaObjects.at[objId, "pixelId"],
                         17450571968473)


class TestNDiaSourcesDiaPlugin(unittest.TestCase):

    def testCalculate(self):
        """Test that the number of DiaSources is correct.
        """

        for n_sources in [1, 8, 10]:
            # Test expected number of sources per object.
            objId = 0
            diaSources = pd.DataFrame(
                data={"diaObjectId": n_sources * [objId],
                      "band": n_sources * ["g"],
                      "diaSourceId": np.arange(n_sources, dtype=int)})
            plug = NumDiaSourcesDiaPlugin(NumDiaSourcesDiaPluginConfig(),
                                          "ap_nDiaSources",
                                          None)
            diaObjects = make_diaObject_table(objId, plug, default_value=int(-1))
            run_multi_plugin(diaObjects, diaSources, "g", plug)

            self.assertEqual(n_sources, diaObjects.at[objId, "nDiaSources"])
            self.assertEqual(diaObjects["nDiaSources"].dtype, np.int64)


class TestSimpleSourceFlagDiaPlugin(unittest.TestCase):

    def testCalculate(self):
        """Test that DiaObject flags are set.
        """
        objId = 0
        n_sources = 10

        # Test expected flags, no flags set.
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.zeros(n_sources, dtype=np.uint64)})
        plug = SimpleSourceFlagDiaPlugin(SimpleSourceFlagDiaPluginConfig(),
                                         "ap_diaObjectFlag",
                                         None)

        diaObjects = make_diaObject_table(objId, plug, default_value=np.uint64(0))
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 0)
        self.assertEqual(diaObjects["flags"].dtype, np.uint64)

        # Test expected flags, all flags set.
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.ones(n_sources, dtype=np.uint64)})
        diaObjects = make_diaObject_table(objId, plug, default_value=np.uint64(0))
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)
        self.assertEqual(diaObjects["flags"].dtype, np.uint64)

        # Test expected flags, random flags.
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.random.randint(0, 2 ** 16, size=n_sources)})

        diaObjects = make_diaObject_table(objId, plug, default_value=np.uint64(0))
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)
        self.assertEqual(diaObjects["flags"].dtype, np.uint64)

        # Test expected flags, one flag set.
        flag_array = np.zeros(n_sources, dtype=np.uint64)
        flag_array[4] = 256
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": flag_array})
        diaObjects = make_diaObject_table(objId, plug, default_value=np.uint64(0))
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)
        self.assertEqual(diaObjects["flags"].dtype, np.uint64)


class TestWeightedMeanDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test mean value calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected mean.
        # In the first test, we have only one object.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": np.linspace(-1, 1, n_sources),
                  "psfFluxErr": np.ones(n_sources)})

        plug = WeightedMeanDiaPsfFlux(WeightedMeanDiaPsfFluxConfig(),
                                      "ap_meanFlux",
                                      None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "u_psfFluxMean"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "u_psfFluxMeanErr"],
                               np.sqrt(1 / n_sources))
        self.assertEqual(diaObjects.loc[objId, "u_psfFluxNdata"], n_sources)
        # We expect this to be converted to float.
        # TODO DM-53254: This should be an integer (and should be checked
        # to be an integer).
        self.assertEqual(diaObjects["u_psfFluxNdata"].dtype, np.float64)

        # Test expected mean with a nan value.
        # In the second test, we have two objects (one empty).
        diaObjects = pd.DataFrame({"diaObjectId": [objId, objId + 1]})
        fluxes = np.linspace(-1, 1, n_sources)
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxMean"],
                               np.nanmean(fluxes))
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxMeanErr"],
                               np.sqrt(1 / (n_sources - 1)))
        self.assertEqual(diaObjects.loc[objId, "r_psfFluxNdata"], n_sources - 1)
        # We expect this to be converted to float.
        # TODO DM-53254: This should be an integer (and should be checked
        # to be an integer).
        self.assertEqual(diaObjects["r_psfFluxNdata"].dtype, np.float64)


class TestPercentileDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux percentile calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected percentile values.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = PercentileDiaPsfFlux(PercentileDiaPsfFluxConfig(),
                                    "ap_percentileFlux",
                                    None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObjects.at[objId, "u_psfFluxPercentile{:02d}".format(pTile)],
                testVal)

        # Test expected percentile values with a nan value.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObjects.at[objId, "r_psfFluxPercentile{:02d}".format(pTile)],
                testVal)


class TestSigmaDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux scatter calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected sigma scatter of fluxes.
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = SigmaDiaPsfFlux(SigmaDiaPsfFluxConfig(),
                               "ap_sigmaFlux",
                               None)
        diaObjects = make_diaObject_table(objId, plug, band='u')
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "u_psfFluxSigma"],
                               np.nanstd(fluxes, ddof=1))

        # test one input, returns nan.
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "band": 1 * ["g"],
                  "diaSourceId": [0],
                  "psfFlux": [fluxes[0]],
                  "psfFluxErr": [1.]})

        diaObjects = make_diaObject_table(objId, plug, band='g')
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "g_psfFluxSigma"]))

        # Test expected sigma scatter of fluxes with a nan value.
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        diaObjects = make_diaObject_table(objId, plug, band='r')
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxSigma"],
                               np.nanstd(fluxes, ddof=1))


class TestChi2DiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux chi2 calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected chi^2 value.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "u_psfFluxMean": [0.0]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = Chi2DiaPsfFlux(Chi2DiaPsfFluxConfig(),
                              "ap_chi2Flux",
                              None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "u_psfFluxChi2"],
            np.nansum(((diaSources["psfFlux"]
                        - np.nanmean(diaSources["psfFlux"]))
                       / diaSources["psfFluxErr"]) ** 2))

        # Test expected chi^2 value with a nan value set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "r_psfFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "r_psfFluxChi2"],
            np.nansum(((diaSources["psfFlux"]
                        - np.nanmean(diaSources["psfFlux"]))
                       / diaSources["psfFluxErr"]) ** 2))


class TestMadDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux median absolute deviation calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected MAD value.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = MadDiaPsfFlux(MadDiaPsfFluxConfig(),
                             "ap_madFlux",
                             None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "u_psfFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))

        # Test expected MAD value with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))


class TestSkewDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux skew calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected skew value.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = SkewDiaPsfFlux(SkewDiaPsfFluxConfig(),
                              "ap_skewFlux",
                              None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "u_psfFluxSkew"],
            skew_wrapper(fluxes))

        # Test expected skew value with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(
            diaObjects.at[objId, "r_psfFluxSkew"],
            skew_wrapper(fluxes))


class TestMinMaxDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux min/max calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected MinMax fluxes.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})

        plug = MinMaxDiaPsfFlux(MinMaxDiaPsfFluxConfig(),
                                "ap_minMaxFlux",
                                None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertEqual(diaObjects.loc[objId, "u_psfFluxMin"], -1)
        self.assertEqual(diaObjects.loc[objId, "u_psfFluxMax"], 1)

        # Test expected MinMax fluxes with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertEqual(diaObjects.loc[objId, "r_psfFluxMin"], -1)
        self.assertEqual(diaObjects.loc[objId, "r_psfFluxMax"], 1)


class TestMaxSlopeDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux maximum slope.
        """
        n_sources = 10
        objId = 0

        # Test max slope value.
        fluxes = np.linspace(-1, 1, n_sources)
        times = np.concatenate([np.linspace(0, 1, n_sources)[:-1], [1 - 1/90]])
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources),
                  "midpointMjdTai": times})

        plug = MaxSlopeDiaPsfFlux(MaxSlopeDiaPsfFluxConfig(),
                                  "ap_maxSlopeFlux",
                                  None)
        diaObjects = make_diaObject_table(objId, plug, band='u')
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "u_psfFluxMaxSlope"], 2 + 2/9)

        # Test max slope value returns nan on 1 input.
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "band": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "psfFlux": fluxes[0],
                  "psfFluxErr": np.ones(1),
                  "midpointMjdTai": times[0]})
        diaObjects = make_diaObject_table(objId, plug, band='g')
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "g_psfFluxMaxSlope"]))

        # Test max slope value inputing nan values.
        fluxes[4] = np.nan
        times[7] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources),
                  "midpointMjdTai": times})
        diaObjects = make_diaObject_table(objId, plug, band='r')
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxMaxSlope"], 2 + 2 / 9)


class TestErrMeanDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test error mean calculation.
        """
        n_sources = 10
        objId = 0

        # Test mean of the errors.
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.linspace(1, 2, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors})

        plug = ErrMeanDiaPsfFlux(ErrMeanDiaPsfFluxConfig(),
                                 "ap_errMeanFlux",
                                 None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "u_psfFluxErrMean"],
                               np.nanmean(errors).astype(np.float32))

        # Test mean of the errors with input nan value.
        errors[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxErrMean"],
                               np.nanmean(errors).astype(np.float32))


class TestLinearFitDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test a linear fit to flux vs time.
        """
        n_sources = 10
        objId = 0

        # Test best fit linear model.
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.linspace(1, 2, n_sources)
        times = np.linspace(0, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors,
                  "midpointMjdTai": times})

        plug = LinearFitDiaPsfFlux(LinearFitDiaPsfFluxConfig(),
                                   "ap_LinearFit",
                                   None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.loc[objId, "u_psfFluxLinearSlope"],
                               2.)
        self.assertAlmostEqual(diaObjects.loc[objId, "u_psfFluxLinearIntercept"],
                               -1.)

        # Test best fit linear model with input nans.
        fluxes[7] = np.nan
        errors[4] = np.nan
        times[2] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors,
                  "midpointMjdTai": times})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.loc[objId, "r_psfFluxLinearSlope"], 2.)
        self.assertAlmostEqual(diaObjects.loc[objId, "r_psfFluxLinearIntercept"],
                               -1.)


class TestStetsonJDiaPsfFlux(unittest.TestCase):

    def testCalculate(self):
        """Test the stetsonJ statistic.
        """
        n_sources = 10
        objId = 0

        # Test stetsonJ calculation.
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.ones(n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "u_psfFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors})

        plug = StetsonJDiaPsfFlux(StetsonJDiaPsfFluxConfig(),
                                  "ap_StetsonJ",
                                  None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        # Expected StetsonJ for the values created. Confirmed using Cesimum's
        # implementation. http://github.com/cesium-ml/cesium
        self.assertAlmostEqual(diaObjects.loc[objId, "u_psfFluxStetsonJ"],
                               -0.5958393936080928)

        # Test stetsonJ calculation returns nan on single input.
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "g_psfFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "band": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "psfFlux": fluxes[0],
                  "psfFluxErr": errors[0]})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "g_psfFluxStetsonJ"]))

        # Test stetsonJ calculation returns when nans are input.
        fluxes[7] = np.nan
        errors[4] = np.nan
        nonNanMask = np.logical_and(~np.isnan(fluxes),
                                    ~np.isnan(errors))
        diaObjects = pd.DataFrame(
            {"diaObjectId": [objId],
             "r_psfFluxMean": [np.average(fluxes[nonNanMask],
                                          weights=errors[nonNanMask])]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psfFlux": fluxes,
                  "psfFluxErr": errors})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_psfFluxStetsonJ"],
                               -0.5412797916187173)


class TestWeightedMeanDiaTotFlux(unittest.TestCase):

    def testCalculate(self):
        """Test mean value calculation.
        """
        n_sources = 10
        objId = 0

        # Test test mean on scienceFlux.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "scienceFlux": fluxes,
                  "scienceFluxErr": np.ones(n_sources)})

        plug = WeightedMeanDiaTotFlux(WeightedMeanDiaTotFluxConfig(),
                                      "ap_meanTotFlux",
                                      None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "u_scienceFluxMean"], 0.0)
        self.assertAlmostEqual(diaObjects.at[objId, "u_scienceFluxMeanErr"],
                               np.sqrt(1 / n_sources))

        # Test test mean on scienceFlux with input nans
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "scienceFlux": fluxes,
                  "scienceFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "r_scienceFluxMean"],
                               np.nanmean(fluxes))
        self.assertAlmostEqual(diaObjects.at[objId, "r_scienceFluxMeanErr"],
                               np.sqrt(1 / (n_sources - 1)))


def generatePeriodicData(n=10, period=10):
    """Generate noisy, sinusoidally-varying periodic data for testing Lomb-
    Scargle Periodogram.

    The returned fluxes will have, within the errors, the passed-in period and
    a power close to 1, because the fluxes are purely sinusoidal.

    Parameters
    ----------
    n : int
        Number of data points to generate.
    period : float
        Period of the periodic signal.

    Returns
    -------
    t : np.ndarray
        Time values.
    y_obs : np.ndarray
        Observed flux values.
    """
    np.random.seed(42)

    t = np.linspace(-2*np.pi, 2*np.pi, n) + 100*np.random.random(n)
    y = 10 + np.sin(2 * np.pi * t / period)
    y_obs = np.random.normal(y, 0.001)

    return t, y_obs


class TestMultiLombScarglePeriodogram(lsst.utils.tests.TestCase):

    def testCalculate(self):
        """Test Mulitband Lomb Scargle Periodogram."""
        n_sources = 10
        objId = 0

        # Create synthetic multi-band data
        times, fluxes = generatePeriodicData(n_sources, period=10)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources//2 * ["u"] + n_sources//2 * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "midpointMjdTai": times,
                  "psfFlux": fluxes,
                  "psfFluxErr": 1e-3+np.zeros(n_sources)})

        plugin = LombScarglePeriodogramMulti(LombScarglePeriodogramMultiConfig(),
                                             "ap_lombScarglePeriodogramMulti",
                                             None)

        run_multiband_plugin(diaObjects, diaSources, plugin)
        self.assertAlmostEqual(diaObjects.at[objId, "multiPeriod"], 10, delta=0.04)
        self.assertAlmostEqual(diaObjects.at[objId, "multiPower"], 1, delta=1e-2)
        # This implementation of LS returns a normalized power < 1.
        self.assertLess(diaObjects.at[objId, "multiPower"], 1)
        self.assertAlmostEqual(diaObjects.at[objId, "multiFap"], 0, delta=0.04)
        # Note: The below values are empirical, but seem reasonable, and
        # test that we get values for each band.
        self.assertAlmostEqual(diaObjects.at[objId, "u_multiAmp"], 0.029, delta=0.01)
        self.assertAlmostEqual(diaObjects.at[objId, "g_multiAmp"], 0.029, delta=0.01)
        self.assertAlmostEqual(diaObjects.at[objId, "u_multiPhase"], -2.0, delta=0.2)
        self.assertAlmostEqual(diaObjects.at[objId, "g_multiPhase"], 1.0, delta=0.1)

    def testCalculateTwoSources(self):
        """Test Mulitband Lomb Scargle Periodogram with 2 sources (minimum
        detections = 5), which will result in NaN output."""
        objId = 0
        n_sources = 2
        times, fluxes = generatePeriodicData(n_sources, period=10)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "midpointMjdTai": times,
                  "psfFlux": fluxes,
                  "psfFluxErr": 1e-3+np.zeros(n_sources)})

        plugin = LombScarglePeriodogramMulti(LombScarglePeriodogramMultiConfig(),
                                             "ap_lombScarglePeriodogramMulti",
                                             None)

        run_multi_plugin(diaObjects, diaSources, "u", plugin)
        self.assertTrue(np.isnan(diaObjects.at[objId, "multiPeriod"]))
        self.assertTrue(np.isnan(diaObjects.at[objId, "multiPower"]))
        self.assertTrue(np.isnan(diaObjects.at[objId, "multiFap"]))


class TestLombScarglePeriodogram(lsst.utils.tests.TestCase):

    def testCalculate(self):
        """Test Lomb Scargle Periodogram."""
        n_sources = 10
        objId = 0

        # Test period calculation.
        times, fluxes = generatePeriodicData(n_sources, period=10)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "midpointMjdTai": times,
                  "psfFlux": fluxes,
                  "psfFluxErr": 1e-3+np.zeros(n_sources)})

        plugin = LombScarglePeriodogram(LombScarglePeriodogramConfig(),
                                        "ap_lombScarglePeriodogram",
                                        None)

        run_multi_plugin(diaObjects, diaSources, "u", plugin)
        self.assertAlmostEqual(diaObjects.at[objId, "u_period"], 10, delta=0.04)
        # This implementation of LS returns a normalized power < 1.
        self.assertAlmostEqual(diaObjects.at[objId, "u_power"], 1, delta=1e-2)
        self.assertLess(diaObjects.at[objId, "u_power"], 1)

        # Test that we get the same result with a NaN flux.
        diaSources.loc[4, "psfFlux"] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "midpointMjdTai": times,
                  "psfFlux": fluxes,
                  "psfFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plugin)
        self.assertAlmostEqual(diaObjects.at[objId, "r_period"], 10, delta=0.04)
        self.assertAlmostEqual(diaObjects.at[objId, "r_power"], 1, delta=1e-2)
        # This implementation of LS returns a normalized power < 1.
        self.assertLess(diaObjects.at[objId, "r_power"], 1)


class TestSigmaDiaTotFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux scatter calculation.
        """
        n_sources = 10
        objId = 0

        # Test test scatter on scienceFlux.
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "scienceFlux": fluxes,
                  "scienceFluxErr": np.ones(n_sources)})

        plug = SigmaDiaTotFlux(SigmaDiaTotFluxConfig(),
                               "ap_sigmaTotFlux",
                               None)
        diaObjects = make_diaObject_table(objId, plug, band='u')
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "u_scienceFluxSigma"],
                               np.nanstd(fluxes, ddof=1))

        # Test test scatter on scienceFlux returns nan on 1 input.
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "band": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "scienceFlux": fluxes[0],
                  "scienceFluxErr": np.ones(1)})
        diaObjects = make_diaObject_table(objId, plug, band='g')
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "g_scienceFluxSigma"]))

        # Test test scatter on scienceFlux takes input nans.
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "band": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "scienceFlux": fluxes,
                  "scienceFluxErr": np.ones(n_sources)})
        diaObjects = make_diaObject_table(objId, plug, band='r')
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "r_scienceFluxSigma"],
                               np.nanstd(fluxes, ddof=1))


def skew_wrapper(values):
    """Compute scipy skew, omitting nans.

    This version works with both scipy<1.9 (where it erroneously returns a
    masked array) and scipy>=1.9 (where it correctly returns a float).

    Parameters
    ----------
    values : `np.ndarray`

    Returns
    -------
    skew_value : `float`
    """
    value = skew(values, bias=False, nan_policy="omit")
    if isinstance(value, np.ma.masked_array):
        return value.data
    else:
        return value


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
