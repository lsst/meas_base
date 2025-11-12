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
    LombScarglePeriodogramMulti, LombScarglePeriodogramMultiConfig)
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
        """
        n_sources = 10
        objId = 0

        plug = MeanDiaPosition(MeanDiaPositionConfig(),
                               "ap_meanPosition",
                               None)

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
