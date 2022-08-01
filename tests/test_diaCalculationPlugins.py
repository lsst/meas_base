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
    WeightedMeanDiaPsFlux, WeightedMeanDiaPsFluxConfig,
    PercentileDiaPsFlux, PercentileDiaPsFluxConfig,
    SigmaDiaPsFlux, SigmaDiaPsFluxConfig,
    Chi2DiaPsFlux, Chi2DiaPsFluxConfig,
    MadDiaPsFlux, MadDiaPsFluxConfig,
    SkewDiaPsFlux, SkewDiaPsFluxConfig,
    MinMaxDiaPsFlux, MinMaxDiaPsFluxConfig,
    MaxSlopeDiaPsFlux, MaxSlopeDiaPsFluxConfig,
    ErrMeanDiaPsFlux, ErrMeanDiaPsFluxConfig,
    LinearFitDiaPsFlux, LinearFitDiaPsFluxConfig,
    StetsonJDiaPsFlux, StetsonJDiaPsFluxConfig,
    WeightedMeanDiaTotFlux, WeightedMeanDiaTotFluxConfig,
    SigmaDiaTotFlux, SigmaDiaTotFluxConfig)
import lsst.utils.tests


def run_single_plugin(diaObjectCat,
                      diaObjectId,
                      diaSourceCat,
                      filterName,
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
        ["diaObjectId", "filterName", "diaSourceId"],
        inplace=True,
        drop=False)

    objDiaSources = diaSourceCat.loc[diaObjectId]
    updatingFilterDiaSources = diaSourceCat.loc[
        (diaObjectId, filterName), :
    ]

    plugin.calculate(diaObjects=diaObjectCat,
                     diaObjectId=diaObjectId,
                     diaSources=objDiaSources,
                     filterDiaSources=updatingFilterDiaSources,
                     filterName=filterName)


def run_multi_plugin(diaObjectCat, diaSourceCat, filterName, plugin):
    """Wrapper for running multi plugins.

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
        ["diaObjectId", "filterName", "diaSourceId"],
        inplace=True,
        drop=False)

    updatingFilterDiaSources = diaSourceCat.loc[
        (slice(None), filterName), :
    ]

    diaSourcesGB = diaSourceCat.groupby(level=0)
    filterDiaSourcesGB = updatingFilterDiaSources.groupby(level=0)

    plugin.calculate(diaObjects=diaObjectCat,
                     diaSources=diaSourcesGB,
                     filterDiaSources=filterDiaSourcesGB,
                     filterName=filterName)


class TestMeanPosition(unittest.TestCase):

    def testCalculate(self):
        """Test mean position calculation.
        """
        n_sources = 10
        objId = 0

        plug = MeanDiaPosition(MeanDiaPositionConfig(),
                               "ap_meanPosition",
                               None)

        # Test expected means in RA.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1, 1, n_sources),
                                        "decl": np.zeros(n_sources),
                                        "midPointTai": np.linspace(0,
                                                                   n_sources,
                                                                   n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "filterName": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "decl"], 0.0)
        self.assertEqual(diaObjects.loc[objId, "radecTai"], 10)

        # Test expected means in DEC.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.linspace(-1, 1, n_sources),
                                        "midPointTai": np.linspace(0,
                                                                   n_sources,
                                                                   n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "filterName": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "ra"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "decl"], 0.0)
        self.assertEqual(diaObjects.loc[objId, "radecTai"], 10)

        # Test failure mode RA is nan.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.full(n_sources, np.nan),
                                        "decl": np.zeros(n_sources),
                                        "midPointTai": np.linspace(0,
                                                                   n_sources,
                                                                   n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "filterName": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "decl"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "radecTai"]))

        # Test failure mode DEC is nan.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.full(n_sources, np.nan),
                                        "midPointTai": np.linspace(0,
                                                                   n_sources,
                                                                   n_sources),
                                        "diaObjectId": n_sources * [objId],
                                        "filterName": n_sources * ["g"],
                                        "diaSourceId": np.arange(n_sources,
                                                                 dtype=int)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)

        self.assertTrue(np.isnan(diaObjects.loc[objId, "ra"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "decl"]))
        self.assertTrue(np.isnan(diaObjects.loc[objId, "radecTai"]))


class TestHTMIndexPosition(unittest.TestCase):

    def testCalculate(self):
        """Test HTMPixel assignment calculation.
        """
        # Test expected pixelId at RA, DEC = 0
        objId = 0
        n_sources = 10
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaObjects.loc[objId, "ra"] = 0.
        diaObjects.loc[objId, "decl"] = 0.
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int)})
        plug = HTMIndexDiaPosition(HTMIndexDiaPositionConfig(),
                                   "ap_HTMIndex",
                                   None)

        run_single_plugin(diaObjectCat=diaObjects,
                          diaObjectId=objId,
                          diaSourceCat=diaSources,
                          filterName="g",
                          plugin=plug)
        self.assertEqual(diaObjects.at[objId, "pixelId"],
                         17042430230528)

        # Test expected pixelId at some value of RA and DEC.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaObjects.loc[objId, "ra"] = 45.37
        diaObjects.loc[objId, "decl"] = 13.67
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int)})
        run_single_plugin(diaObjectCat=diaObjects,
                          diaObjectId=objId,
                          diaSourceCat=diaSources,
                          filterName="g",
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
            diaObjects = pd.DataFrame({"diaObjectId": [objId]})
            diaSources = pd.DataFrame(
                data={"diaObjectId": n_sources * [objId],
                      "filterName": n_sources * ["g"],
                      "diaSourceId": np.arange(n_sources, dtype=int)})
            plug = NumDiaSourcesDiaPlugin(NumDiaSourcesDiaPluginConfig(),
                                          "ap_nDiaSources",
                                          None)
            run_multi_plugin(diaObjects, diaSources, "g", plug)

            self.assertEqual(n_sources, diaObjects.at[objId, "nDiaSources"])


class TestSimpleSourceFlagDiaPlugin(unittest.TestCase):

    def testCalculate(self):
        """Test that DiaObject flags are set.
        """
        objId = 0
        n_sources = 10

        # Test expected flags, no flags set.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.zeros(n_sources, dtype=np.uint64)})
        plug = SimpleSourceFlagDiaPlugin(SimpleSourceFlagDiaPluginConfig(),
                                         "ap_diaObjectFlag",
                                         None)
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 0)

        # Test expected flags, all flags set.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.ones(n_sources, dtype=np.uint64)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)

        # Test expected flags, random flags.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": np.random.randint(0, 2 ** 16, size=n_sources)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)

        # Test expected flags, one flag set.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        flag_array = np.zeros(n_sources, dtype=np.uint64)
        flag_array[4] = 256
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["g"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "flags": flag_array})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertEqual(diaObjects.at[objId, "flags"], 1)


class TestWeightedMeanDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test mean value calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected mean.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": np.linspace(-1, 1, n_sources),
                  "psFluxErr": np.ones(n_sources)})

        plug = WeightedMeanDiaPsFlux(WeightedMeanDiaPsFluxConfig(),
                                     "ap_meanFlux",
                                     None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)

        self.assertAlmostEqual(diaObjects.loc[objId, "uPSFluxMean"], 0.0)
        self.assertAlmostEqual(diaObjects.loc[objId, "uPSFluxMeanErr"],
                               np.sqrt(1 / n_sources))
        self.assertEqual(diaObjects.loc[objId, "uPSFluxNdata"], n_sources)

        # Test expected mean with a nan value.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        fluxes = np.linspace(-1, 1, n_sources)
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxMean"],
                               np.nanmean(fluxes))
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxMeanErr"],
                               np.sqrt(1 / (n_sources - 1)))
        self.assertEqual(diaObjects.loc[objId, "rPSFluxNdata"], n_sources - 1)


class TestPercentileDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = PercentileDiaPsFlux(PercentileDiaPsFluxConfig(),
                                   "ap_percentileFlux",
                                   None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObjects.at[objId, "uPSFluxPercentile{:02d}".format(pTile)],
                testVal)

        # Test expected percentile values with a nan value.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObjects.at[objId, "rPSFluxPercentile{:02d}".format(pTile)],
                testVal)


class TestSigmaDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux scatter calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected sigma scatter of fluxes.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = SigmaDiaPsFlux(SigmaDiaPsFluxConfig(),
                              "ap_sigmaFlux",
                              None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "uPSFluxSigma"],
                               np.nanstd(fluxes, ddof=1))

        # test one input, returns nan.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "filterName": 1 * ["g"],
                  "diaSourceId": [0],
                  "psFlux": [fluxes[0]],
                  "psFluxErr": [1.]})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "gPSFluxSigma"]))

        # Test expected sigma scatter of fluxes with a nan value.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxSigma"],
                               np.nanstd(fluxes, ddof=1))


class TestChi2DiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux chi2 calculation.
        """
        n_sources = 10
        objId = 0

        # Test expected chi^2 value.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "uPSFluxMean": [0.0]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = Chi2DiaPsFlux(Chi2DiaPsFluxConfig(),
                             "ap_chi2Flux",
                             None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "uPSFluxChi2"],
            np.nansum(((diaSources["psFlux"]
                        - np.nanmean(diaSources["psFlux"]))
                       / diaSources["psFluxErr"]) ** 2))

        # Test expected chi^2 value with a nan value set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "rPSFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "rPSFluxChi2"],
            np.nansum(((diaSources["psFlux"]
                        - np.nanmean(diaSources["psFlux"]))
                       / diaSources["psFluxErr"]) ** 2))


class TestMadDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = MadDiaPsFlux(MadDiaPsFluxConfig(),
                            "ap_madFlux",
                            None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "uPSFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))

        # Test expected MAD value with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))


class TestSkewDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = SkewDiaPsFlux(SkewDiaPsFluxConfig(),
                             "ap_skewFlux",
                             None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(
            diaObjects.loc[objId, "uPSFluxSkew"],
            skew_wrapper(fluxes))

        # Test expected skew value with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(
            diaObjects.at[objId, "rPSFluxSkew"],
            skew_wrapper(fluxes))


class TestMinMaxDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})

        plug = MinMaxDiaPsFlux(MinMaxDiaPsFluxConfig(),
                               "ap_minMaxFlux",
                               None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertEqual(diaObjects.loc[objId, "uPSFluxMin"], -1)
        self.assertEqual(diaObjects.loc[objId, "uPSFluxMax"], 1)

        # Test expected MinMax fluxes with a nan set.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertEqual(diaObjects.loc[objId, "rPSFluxMin"], -1)
        self.assertEqual(diaObjects.loc[objId, "rPSFluxMax"], 1)


class TestMaxSlopeDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux maximum slope.
        """
        n_sources = 10
        objId = 0

        # Test max slope value.
        fluxes = np.linspace(-1, 1, n_sources)
        times = np.concatenate([np.linspace(0, 1, n_sources)[:-1], [1 - 1/90]])
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources),
                  "midPointTai": times})

        plug = MaxSlopeDiaPsFlux(MaxSlopeDiaPsFluxConfig(),
                                 "ap_maxSlopeFlux",
                                 None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "uPSFluxMaxSlope"], 2 + 2/9)

        # Test max slope value returns nan on 1 input.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "filterName": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "psFlux": fluxes[0],
                  "psFluxErr": np.ones(1),
                  "midPointTai": times[0]})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "gPSFluxMaxSlope"]))

        # Test max slope value inputing nan values.
        fluxes[4] = np.nan
        times[7] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": np.ones(n_sources),
                  "midPointTai": times})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxMaxSlope"], 2 + 2 / 9)


class TestErrMeanDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors})

        plug = ErrMeanDiaPsFlux(ErrMeanDiaPsFluxConfig(),
                                "ap_errMeanFlux",
                                None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "uPSFluxErrMean"],
                               np.nanmean(errors))

        # Test mean of the errors with input nan value.
        errors[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxErrMean"],
                               np.nanmean(errors))


class TestLinearFitDiaPsFlux(unittest.TestCase):

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
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors,
                  "midPointTai": times})

        plug = LinearFitDiaPsFlux(LinearFitDiaPsFluxConfig(),
                                  "ap_LinearFit",
                                  None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.loc[objId, "uPSFluxLinearSlope"],
                               2.)
        self.assertAlmostEqual(diaObjects.loc[objId, "uPSFluxLinearIntercept"],
                               -1.)

        # Test best fit linear model with input nans.
        fluxes[7] = np.nan
        errors[4] = np.nan
        times[2] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors,
                  "midPointTai": times})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.loc[objId, "rPSFluxLinearSlope"], 2.)
        self.assertAlmostEqual(diaObjects.loc[objId, "rPSFluxLinearIntercept"],
                               -1.)


class TestStetsonJDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test the stetsonJ statistic.
        """
        n_sources = 10
        objId = 0

        # Test stetsonJ calculation.
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.ones(n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "uPSFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors})

        plug = StetsonJDiaPsFlux(StetsonJDiaPsFluxConfig(),
                                 "ap_StetsonJ",
                                 None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        # Expected StetsonJ for the values created. Confirmed using Cesimum's
        # implementation. http://github.com/cesium-ml/cesium
        self.assertAlmostEqual(diaObjects.loc[objId, "uPSFluxStetsonJ"],
                               -0.5958393936080928)

        # Test stetsonJ calculation returns nan on single input.
        diaObjects = pd.DataFrame({"diaObjectId": [objId],
                                   "gPSFluxMean": [np.nanmean(fluxes)]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "filterName": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "psFlux": fluxes[0],
                  "psFluxErr": errors[0]})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "gPSFluxStetsonJ"]))

        # Test stetsonJ calculation returns when nans are input.
        fluxes[7] = np.nan
        errors[4] = np.nan
        nonNanMask = np.logical_and(~np.isnan(fluxes),
                                    ~np.isnan(errors))
        diaObjects = pd.DataFrame(
            {"diaObjectId": [objId],
             "rPSFluxMean": [np.average(fluxes[nonNanMask],
                                        weights=errors[nonNanMask])]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "psFlux": fluxes,
                  "psFluxErr": errors})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rPSFluxStetsonJ"],
                               -0.5412797916187173)


class TestWeightedMeanDiaTotFlux(unittest.TestCase):

    def testCalculate(self):
        """Test mean value calculation.
        """
        n_sources = 10
        objId = 0

        # Test test mean on totFlux.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "totFlux": fluxes,
                  "totFluxErr": np.ones(n_sources)})

        plug = WeightedMeanDiaTotFlux(WeightedMeanDiaTotFluxConfig(),
                                      "ap_meanTotFlux",
                                      None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "uTOTFluxMean"], 0.0)
        self.assertAlmostEqual(diaObjects.at[objId, "uTOTFluxMeanErr"],
                               np.sqrt(1 / n_sources))

        # Test test mean on totFlux with input nans
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "totFlux": fluxes,
                  "totFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)

        self.assertAlmostEqual(diaObjects.at[objId, "rTOTFluxMean"],
                               np.nanmean(fluxes))
        self.assertAlmostEqual(diaObjects.at[objId, "rTOTFluxMeanErr"],
                               np.sqrt(1 / (n_sources - 1)))


class TestSigmaDiaTotFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux scatter calculation.
        """
        n_sources = 10
        objId = 0

        # Test test scatter on totFlux.
        fluxes = np.linspace(-1, 1, n_sources)
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["u"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "totFlux": fluxes,
                  "totFluxErr": np.ones(n_sources)})

        plug = SigmaDiaTotFlux(SigmaDiaTotFluxConfig(),
                               "ap_sigmaTotFlux",
                               None)
        run_multi_plugin(diaObjects, diaSources, "u", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "uTOTFluxSigma"],
                               np.nanstd(fluxes, ddof=1))

        # Test test scatter on totFlux returns nan on 1 input.
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": 1 * [objId],
                  "filterName": 1 * ["g"],
                  "diaSourceId": np.arange(1, dtype=int),
                  "totFlux": fluxes[0],
                  "totFluxErr": np.ones(1)})
        run_multi_plugin(diaObjects, diaSources, "g", plug)
        self.assertTrue(np.isnan(diaObjects.at[objId, "gTOTFluxSigma"]))

        # Test test scatter on totFlux takes input nans.
        fluxes[4] = np.nan
        diaObjects = pd.DataFrame({"diaObjectId": [objId]})
        diaSources = pd.DataFrame(
            data={"diaObjectId": n_sources * [objId],
                  "filterName": n_sources * ["r"],
                  "diaSourceId": np.arange(n_sources, dtype=int),
                  "totFlux": fluxes,
                  "totFluxErr": np.ones(n_sources)})
        run_multi_plugin(diaObjects, diaSources, "r", plug)
        self.assertAlmostEqual(diaObjects.at[objId, "rTOTFluxSigma"],
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
