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

"""Plugins for use in DiaSource summary statistics.

Output columns must be
as defined in the schema of the Apdb both in name and units.
"""

import functools
import warnings

from astropy.stats import median_absolute_deviation
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.sphgeom as sphgeom
from astropy.timeseries import LombScargle
from astropy.timeseries import LombScargleMultiband
import math
import statistics

from .diaCalculation import (
    DiaObjectCalculationPluginConfig,
    DiaObjectCalculationPlugin)
from .pluginRegistry import register


__all__ = ("MeanDiaPositionConfig", "MeanDiaPosition",
           "HTMIndexDiaPosition", "HTMIndexDiaPositionConfig",
           "NumDiaSourcesDiaPlugin", "NumDiaSourcesDiaPluginConfig",
           "SimpleSourceFlagDiaPlugin", "SimpleSourceFlagDiaPluginConfig",
           "WeightedMeanDiaPsfFluxConfig", "WeightedMeanDiaPsfFlux",
           "PercentileDiaPsfFlux", "PercentileDiaPsfFluxConfig",
           "SigmaDiaPsfFlux", "SigmaDiaPsfFluxConfig",
           "Chi2DiaPsfFlux", "Chi2DiaPsfFluxConfig",
           "MadDiaPsfFlux", "MadDiaPsfFluxConfig",
           "SkewDiaPsfFlux", "SkewDiaPsfFluxConfig",
           "MinMaxDiaPsfFlux", "MinMaxDiaPsfFluxConfig",
           "MaxSlopeDiaPsfFlux", "MaxSlopeDiaPsfFluxConfig",
           "ErrMeanDiaPsfFlux", "ErrMeanDiaPsfFluxConfig",
           "LinearFitDiaPsfFlux", "LinearFitDiaPsfFluxConfig",
           "StetsonJDiaPsfFlux", "StetsonJDiaPsfFluxConfig",
           "WeightedMeanDiaTotFlux", "WeightedMeanDiaTotFluxConfig",
           "SigmaDiaTotFlux", "SigmaDiaTotFluxConfig",
           "LombScarglePeriodogram", "LombScarglePeriodogramConfig",
           "LombScarglePeriodogramMulti", "LombScarglePeriodogramMultiConfig")


def catchWarnings(_func=None, *, warns=[]):
    """Decorator for generically catching numpy warnings.
    """
    def decoratorCatchWarnings(func):
        @functools.wraps(func)
        def wrapperCatchWarnings(*args, **kwargs):
            with warnings.catch_warnings():
                for val in warns:
                    warnings.filterwarnings("ignore", val)
                return func(*args, **kwargs)
        return wrapperCatchWarnings

    if _func is None:
        return decoratorCatchWarnings
    else:
        return decoratorCatchWarnings(_func)


def compute_optimized_periodogram_grid(x0, oversampling_factor=5, nyquist_factor=100):
    """
    Computes an optimized periodogram frequency grid for a given time series.

    Parameters
    ----------
    x0 : `array`
        The input time axis.
    oversampling_factor : `int`, optional
        The oversampling factor for frequency grid.
    nyquist_factor : `int`, optional
        The Nyquist factor for frequency grid.

    Returns
    -------
    frequencies : `array`
        The computed optimized periodogram frequency grid.
    """

    num_points = len(x0)
    baseline = np.max(x0) - np.min(x0)

    # Calculate the frequency resolution based on oversampling factor and baseline
    frequency_resolution = 1. / baseline / oversampling_factor

    num_frequencies = int(
        0.5 * oversampling_factor * nyquist_factor * num_points)
    frequencies = frequency_resolution + \
        frequency_resolution * np.arange(num_frequencies)

    return frequencies


class LombScarglePeriodogramConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_lombScarglePeriodogram")
class LombScarglePeriodogram(DiaObjectCalculationPlugin):
    """Compute the single-band period of a DiaObject given a set of DiaSources.
    """
    ConfigClass = LombScarglePeriodogramConfig

    plugType = "multi"
    outputCols = ["period", "power"]
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band):
        """Compute the periodogram.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Summary objects to store values in.
        diaSources : `pandas.DataFrame` or `pandas.DataFrameGroupBy`
            Catalog of DiaSources summarized by this DiaObject.
        """

        # Check and initialize output columns in diaObjects.
        if (periodCol := f"{band}_period") not in diaObjects.columns:
            diaObjects[periodCol] = np.nan
        if (powerCol := f"{band}_power") not in diaObjects.columns:
            diaObjects[powerCol] = np.nan

        def _calculate_period(df, min_detections=5, nterms=1, oversampling_factor=5, nyquist_factor=100):
            """Compute the Lomb-Scargle periodogram given a set of DiaSources.

            Parameters
            ----------
            df : `pandas.DataFrame`
                The input DataFrame.
            min_detections : `int`, optional
                The minimum number of detections.
            nterms : `int`, optional
                The number of terms in the Lomb-Scargle model.
            oversampling_factor : `int`, optional
                The oversampling factor for frequency grid.
            nyquist_factor : `int`, optional
                The Nyquist factor for frequency grid.

            Returns
            -------
            pd_tab : `pandas.Series`
                The output DataFrame with the Lomb-Scargle parameters.
            """
            tmpDf = df[~np.logical_or(np.isnan(df["psfFlux"]),
                                      np.isnan(df["midpointMjdTai"]))]

            if len(tmpDf) < min_detections:
                return pd.Series({periodCol: np.nan, powerCol: np.nan})

            time = tmpDf["midpointMjdTai"].to_numpy()
            flux = tmpDf["psfFlux"].to_numpy()
            flux_err = tmpDf["psfFluxErr"].to_numpy()

            lsp = LombScargle(time, flux, dy=flux_err, nterms=nterms)
            f_grid = compute_optimized_periodogram_grid(
                time, oversampling_factor=oversampling_factor, nyquist_factor=nyquist_factor)
            period = 1/f_grid
            power = lsp.power(f_grid)

            return pd.Series({periodCol: period[np.argmax(power)],
                              powerCol: np.max(power)})

        diaObjects.loc[:, [periodCol, powerCol]
                       ] = filterDiaSources.apply(_calculate_period)


class LombScarglePeriodogramMultiConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_lombScarglePeriodogramMulti")
class LombScarglePeriodogramMulti(DiaObjectCalculationPlugin):
    """Compute the multi-band LombScargle periodogram of a DiaObject given a set of DiaSources.
    """
    ConfigClass = LombScarglePeriodogramMultiConfig

    plugType = "multi"
    outputCols = ["multiPeriod", "multiPower",
                  "multiFap", "multiAmp", "multiPhase"]
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @staticmethod
    def calculate_baluev_fap(time, n, maxPeriod, zmax):
        """Calculate the False-Alarm probability using the Baluev approximation.

        Parameters
        ----------
        time : `array`
            The input time axis.
        n : `int`
            The number of detections.
        maxPeriod : `float`
            The maximum period in the grid.
        zmax : `float`
            The maximum power in the grid.

        Returns
        -------
        fap_estimate : `float`
            The False-Alarm probability Baluev approximation.

        Notes
        ----------
        .. [1] Baluev, R. V. 2008, MNRAS, 385, 1279
        .. [2] Süveges, M., Guy, L.P., Eyer, L., et al. 2015, MNRAS, 450, 2052
        """
        if n <= 2:
            return np.nan

        gam_ratio = math.factorial(int((n - 1)/2)) / math.factorial(int((n - 2)/2))
        fu = 1/maxPeriod
        return gam_ratio * np.sqrt(
            4*np.pi*statistics.variance(time)
        ) * fu * (1-zmax)**((n-4)/2) * np.sqrt(zmax)

    @staticmethod
    def generate_lsp_params(lsp_model, fbest, bands):
        """Generate the Lomb-Scargle parameters.
        Parameters
        ----------
        lsp_model : `astropy.timeseries.LombScargleMultiband`
            The Lomb-Scargle model.
        fbest : `float`
            The best period.
        bands : `array`
            The bands of the time series.

        Returns
        -------
        Amp : `array`
            The amplitude of the time series.
        Ph : `array`
            The phase of the time series.

        Notes
        ----------
        .. [1] VanderPlas, J. T., & Ivezic, Z. 2015, ApJ, 812, 18
        """
        best_params = lsp_model.model_parameters(fbest, units=True)

        name_params = [f"theta_base_{i}" for i in range(3)]
        name_params += [f"theta_band_{band}_{i}" for band in np.unique(bands) for i in range(3)]

        df_params = pd.DataFrame([best_params], columns=name_params)

        unique_bands = np.unique(bands)

        amplitude_band = [np.sqrt(df_params[f"theta_band_{band}_1"]**2
                                  + df_params[f"theta_band_{band}_2"]**2)
                          for band in unique_bands]
        phase_bands = [np.arctan2(df_params[f"theta_band_{band}_2"],
                                  df_params[f"theta_band_{band}_1"]) for band in unique_bands]

        amp = [a[0] for a in amplitude_band]
        ph = [p[0] for p in phase_bands]

        return amp, ph

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  **kwargs):
        """Compute the multi-band LombScargle periodogram of a DiaObject given
        a set of DiaSources.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Summary objects to store values in.
        diaSources : `pandas.DataFrame` or `pandas.DataFrameGroupBy`
            Catalog of DiaSources summarized by this DiaObject.
        **kwargs : `dict`
            Unused kwargs that are always passed to a plugin.
        """
        n_bands = len(diaSources["band"].unique())
        # Check and initialize output columns in diaObjects.
        if (periodCol := "multiPeriod") not in diaObjects.columns:
            diaObjects[periodCol] = np.nan
        if (powerCol := "multiPower") not in diaObjects.columns:
            diaObjects[powerCol] = np.nan
        if (fapCol := "multiFap") not in diaObjects.columns:
            diaObjects[fapCol] = np.nan
        if (ampCol := "multiAmp") not in diaObjects.columns:
            diaObjects[ampCol] = pd.Series([np.nan]*n_bands, dtype="object")
        if (phaseCol := "multiPhase") not in diaObjects.columns:
            diaObjects[phaseCol] = pd.Series([np.nan]*n_bands, dtype="object")

        def _calculate_period_multi(df, min_detections=9, oversampling_factor=5, nyquist_factor=100):
            """Calculate the multi-band Lomb-Scargle periodogram.

            Parameters
            ----------
            df : `pandas.DataFrame`
                The input DataFrame.
            min_detections : `int`, optional
                The minimum number of detections, including all bands.
            oversampling_factor : `int`, optional
                The oversampling factor for frequency grid.
            nyquist_factor : `int`, optional
                The Nyquist factor for frequency grid.

            Returns
            -------
            pd_tab : `pandas.Series`
                The output DataFrame with the Lomb-Scargle parameters.
            """
            tmpDf = df[~np.logical_or(np.isnan(df["psfFlux"]),
                                      np.isnan(df["midpointMjdTai"]))]

            if (len(tmpDf)) < min_detections:
                return pd.Series({periodCol: np.nan,
                                  powerCol: np.nan,
                                  fapCol: np.nan,
                                  ampCol: pd.Series([np.nan]*n_bands, dtype="object"),
                                  phaseCol: pd.Series([np.nan]*n_bands, dtype="object")})

            time = tmpDf["midpointMjdTai"].to_numpy()
            flux = tmpDf["psfFlux"].to_numpy()
            flux_err = tmpDf["psfFluxErr"].to_numpy()
            bands = tmpDf["band"].to_numpy()

            lsp = LombScargleMultiband(time, flux, bands, dy=flux_err,
                                       nterms_base=1, nterms_band=1)

            f_grid = compute_optimized_periodogram_grid(
                time, oversampling_factor=oversampling_factor, nyquist_factor=nyquist_factor)
            period = 1/f_grid
            power = lsp.power(f_grid)

            fap_estimate = self.calculate_baluev_fap(
                time, len(time), period[np.argmax(power)], np.max(power))

            params_table_new = self.generate_lsp_params(lsp, f_grid[np.argmax(power)], bands)

            pd_tab = pd.Series({periodCol: period[np.argmax(power)],
                                powerCol: np.max(power),
                                fapCol: fap_estimate,
                                ampCol: params_table_new[0],
                                phaseCol: params_table_new[1]
                                })

            return pd_tab

        diaObjects.loc[:, [periodCol, powerCol, fapCol, ampCol, phaseCol]
                       ] = diaSources.apply(_calculate_period_multi)


class MeanDiaPositionConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanPosition")
class MeanDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """

    ConfigClass = MeanDiaPositionConfig

    plugType = 'multi'

    outputCols = ["ra", "dec", "radecMjdTai"]
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Compute the mean ra/dec position of the diaObject given the
        diaSource locations.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Summary objects to store values in.
        diaSources : `pandas.DataFrame` or `pandas.DataFrameGroupBy`
            Catalog of DiaSources summarized by this DiaObject.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        for outCol in self.outputCols:
            if outCol not in diaObjects.columns:
                diaObjects[outCol] = np.nan

        def _computeMeanPos(df):
            aveCoord = geom.averageSpherePoint(
                list(geom.SpherePoint(src["ra"], src["dec"], geom.degrees)
                     for idx, src in df.iterrows()))
            ra = aveCoord.getRa().asDegrees()
            dec = aveCoord.getDec().asDegrees()
            if np.isnan(ra) or np.isnan(dec):
                radecMjdTai = np.nan
            else:
                radecMjdTai = df["midpointMjdTai"].max()

            return pd.Series({"ra": aveCoord.getRa().asDegrees(),
                              "dec": aveCoord.getDec().asDegrees(),
                              "radecMjdTai": radecMjdTai})

        ans = diaSources.apply(_computeMeanPos)
        diaObjects.loc[:, ["ra", "dec", "radecMjdTai"]] = ans


class HTMIndexDiaPositionConfig(DiaObjectCalculationPluginConfig):

    htmLevel = pexConfig.Field(
        dtype=int,
        doc="Level of the HTM pixelization.",
        default=20,
    )


@register("ap_HTMIndex")
class HTMIndexDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.

    Notes
    -----
    This plugin was implemented to satisfy requirements of old APDB interface
    which required ``pixelId`` column in DiaObject with HTM20 index. APDB
    interface had migrated to not need that information, but we keep this
    plugin in case it may be useful for something else.
    """
    ConfigClass = HTMIndexDiaPositionConfig

    plugType = 'single'

    inputCols = ["ra", "dec"]
    outputCols = ["pixelId"]
    needsFilter = False

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)
        self.pixelator = sphgeom.HtmPixelization(self.config.htmLevel)

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self, diaObjects, diaObjectId, **kwargs):
        """Compute the mean position of a DiaObject given a set of DiaSources

        Parameters
        ----------
        diaObjects : `pandas.dataFrame`
            Summary objects to store values in and read ra/dec from.
        diaObjectId : `int`
            Id of the diaObject to update.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        sphPoint = geom.SpherePoint(
            diaObjects.at[diaObjectId, "ra"] * geom.degrees,
            diaObjects.at[diaObjectId, "dec"] * geom.degrees)
        diaObjects.at[diaObjectId, "pixelId"] = self.pixelator.index(
            sphPoint.getVector())


class NumDiaSourcesDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_nDiaSources")
class NumDiaSourcesDiaPlugin(DiaObjectCalculationPlugin):
    """Compute the total number of DiaSources associated with this DiaObject.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["nDiaSources"]
    plugType = "multi"
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Compute the total number of DiaSources associated with this DiaObject.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/dec from.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        dtype = diaObjects["nDiaSources"].dtype
        diaObjects.loc[:, "nDiaSources"] = diaSources.diaObjectId.count().astype(dtype)


class SimpleSourceFlagDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_diaObjectFlag")
class SimpleSourceFlagDiaPlugin(DiaObjectCalculationPlugin):
    """Find if any DiaSource is flagged.

    Set the DiaObject flag if any DiaSource is flagged.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["flags"]
    plugType = "multi"
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Find if any DiaSource is flagged.

        Set the DiaObject flag if any DiaSource is flagged.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/dec from.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        dtype = diaObjects["flags"].dtype
        diaObjects.loc[:, "flags"] = diaSources.flags.any().astype(dtype)


class WeightedMeanDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanFlux")
class WeightedMeanDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    of the DiaSource measured on the difference image.

    Additionally store number of usable data points.
    """

    ConfigClass = WeightedMeanDiaPsfFluxConfig
    outputCols = ["psfFluxMean", "psfFluxMeanErr", "psfFluxNdata"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["invalid value encountered",
                          "divide by zero"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the weighted mean and mean error of the point source flux.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}_psfFluxMean".format(band)
        errName = "{}_psfFluxMeanErr".format(band)
        nDataName = "{}_psfFluxNdata".format(band)
        if meanName not in diaObjects.columns:
            diaObjects[meanName] = np.nan
        if errName not in diaObjects.columns:
            diaObjects[errName] = np.nan
        if nDataName not in diaObjects.columns:
            diaObjects[nDataName] = 0

        def _weightedMean(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psfFlux"]),
                                      np.isnan(df["psfFluxErr"]))]
            tot_weight = np.nansum(1 / tmpDf["psfFluxErr"] ** 2)
            fluxMean = np.nansum(tmpDf["psfFlux"]
                                 / tmpDf["psfFluxErr"] ** 2)
            fluxMean /= tot_weight
            if tot_weight > 0:
                fluxMeanErr = np.sqrt(1 / tot_weight)
            else:
                fluxMeanErr = np.nan
            nFluxData = len(tmpDf)

            return pd.Series({meanName: fluxMean,
                              errName: fluxMeanErr,
                              nDataName: nFluxData},
                             dtype="object")
        df = filterDiaSources.apply(_weightedMean).astype(diaObjects.dtypes[[meanName, errName, nDataName]])

        diaObjects.loc[:, [meanName, errName, nDataName]] = df


class PercentileDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    percentiles = pexConfig.ListField(
        dtype=int,
        default=[5, 25, 50, 75, 95],
        doc="Percentiles to calculate to compute values for. Should be "
            "integer values."
    )


@register("ap_percentileFlux")
class PercentileDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute percentiles of diaSource fluxes.
    """

    ConfigClass = PercentileDiaPsfFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = []
    plugType = "multi"
    needsFilter = True

    def __init__(self, config, name, metadata, **kwargs):
        DiaObjectCalculationPlugin.__init__(self,
                                            config,
                                            name,
                                            metadata,
                                            **kwargs)
        self.outputCols = ["psfFluxPercentile{:02d}".format(percent)
                           for percent in self.config.percentiles]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the percentile fluxes of the point source flux.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        pTileNames = []
        dtype = None
        for tilePercent in self.config.percentiles:
            pTileName = "{}_psfFluxPercentile{:02d}".format(band,
                                                            tilePercent)
            pTileNames.append(pTileName)
            if pTileName not in diaObjects.columns:
                diaObjects[pTileName] = np.nan
            elif dtype is None:
                dtype = diaObjects[pTileName].dtype

        def _fluxPercentiles(df):
            pTiles = np.nanpercentile(df["psfFlux"], self.config.percentiles)
            return pd.Series(
                dict((tileName, pTile)
                     for tileName, pTile in zip(pTileNames, pTiles)))
        if dtype is None:
            diaObjects.loc[:, pTileNames] = filterDiaSources.apply(_fluxPercentiles)
        else:
            diaObjects.loc[:, pTileNames] = filterDiaSources.apply(_fluxPercentiles).astype(dtype)


class SigmaDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaFlux")
class SigmaDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsfFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxSigma"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the sigma fluxes of the point source flux.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        column = "{}_psfFluxSigma".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.psfFlux.std().astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.psfFlux.std()


class Chi2DiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_chi2Flux")
class Chi2DiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute chi2 of diaSource fluxes.
    """

    ConfigClass = Chi2DiaPsfFluxConfig

    # Required input Cols
    inputCols = ["psfFluxMean"]
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxChi2"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the chi2 of the point source fluxes.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}_psfFluxMean".format(band)
        column = "{}_psfFluxChi2".format(band)

        def _chi2(df):
            delta = (df["psfFlux"]
                     - diaObjects.at[df.diaObjectId.iat[0], meanName])
            return np.nansum((delta / df["psfFluxErr"]) ** 2)

        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.apply(_chi2).astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.apply(_chi2)


class MadDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_madFlux")
class MadDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute median absolute deviation of diaSource fluxes.
    """

    ConfigClass = MadDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxMAD"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the median absolute deviation of the point source fluxes.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        column = "{}_psfFluxMAD".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = \
                filterDiaSources.psfFlux.apply(median_absolute_deviation,
                                               ignore_nan=True).astype(dtype)
        else:
            diaObjects.loc[:, column] = \
                filterDiaSources.psfFlux.apply(median_absolute_deviation,
                                               ignore_nan=True)


class SkewDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_skewFlux")
class SkewDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute the skew of diaSource fluxes.
    """

    ConfigClass = SkewDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxSkew"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the skew of the point source fluxes.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        column = "{}_psfFluxSkew".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.psfFlux.skew().astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.psfFlux.skew()


class MinMaxDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_minMaxFlux")
class MinMaxDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute min/max of diaSource fluxes.
    """

    ConfigClass = MinMaxDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxMin", "psfFluxMax"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute min/max of the point source fluxes.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        minName = "{}_psfFluxMin".format(band)
        if minName not in diaObjects.columns:
            diaObjects[minName] = np.nan
        maxName = "{}_psfFluxMax".format(band)
        if maxName not in diaObjects.columns:
            diaObjects[maxName] = np.nan

        dtype = diaObjects[minName].dtype
        diaObjects.loc[:, minName] = filterDiaSources.psfFlux.min().astype(dtype)
        diaObjects.loc[:, maxName] = filterDiaSources.psfFlux.max().astype(dtype)


class MaxSlopeDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_maxSlopeFlux")
class MaxSlopeDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute the maximum ratio time ordered deltaFlux / deltaTime.
    """

    ConfigClass = MinMaxDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxMaxSlope"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the maximum ratio time ordered deltaFlux / deltaTime.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """

        def _maxSlope(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psfFlux"]),
                                      np.isnan(df["midpointMjdTai"]))]
            if len(tmpDf) < 2:
                return np.nan
            times = tmpDf["midpointMjdTai"].to_numpy()
            timeArgs = times.argsort()
            times = times[timeArgs]
            fluxes = tmpDf["psfFlux"].to_numpy()[timeArgs]
            return (np.diff(fluxes) / np.diff(times)).max()

        column = "{}_psfFluxMaxSlope".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.apply(_maxSlope).astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.apply(_maxSlope)


class ErrMeanDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanErrFlux")
class ErrMeanDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute the mean of the dia source errors.
    """

    ConfigClass = ErrMeanDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxErrMean"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the mean of the dia source errors.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        column = "{}_psfFluxErrMean".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.psfFluxErr.mean().astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.psfFluxErr.mean()


class LinearFitDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_linearFit")
class LinearFitDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute fit a linear model to flux vs time.
    """

    ConfigClass = LinearFitDiaPsfFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxLinearSlope", "psfFluxLinearIntercept"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute fit a linear model to flux vs time.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """

        mName = "{}_psfFluxLinearSlope".format(band)
        if mName not in diaObjects.columns:
            diaObjects[mName] = np.nan
        bName = "{}_psfFluxLinearIntercept".format(band)
        if bName not in diaObjects.columns:
            diaObjects[bName] = np.nan
        dtype = diaObjects[mName].dtype

        def _linearFit(df):
            tmpDf = df[~np.logical_or(
                np.isnan(df["psfFlux"]),
                np.logical_or(np.isnan(df["psfFluxErr"]),
                              np.isnan(df["midpointMjdTai"])))]
            if len(tmpDf) < 2:
                return pd.Series({mName: np.nan, bName: np.nan})
            fluxes = tmpDf["psfFlux"].to_numpy()
            errors = tmpDf["psfFluxErr"].to_numpy()
            times = tmpDf["midpointMjdTai"].to_numpy()
            A = np.array([times / errors, 1 / errors]).transpose()
            m, b = lsq_linear(A, fluxes / errors).x
            return pd.Series({mName: m, bName: b}, dtype=dtype)

        diaObjects.loc[:, [mName, bName]] = filterDiaSources.apply(_linearFit)


class StetsonJDiaPsfFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_stetsonJ")
class StetsonJDiaPsfFlux(DiaObjectCalculationPlugin):
    """Compute the StetsonJ statistic on the DIA point source fluxes.
    """

    ConfigClass = LinearFitDiaPsfFluxConfig

    # Required input Cols
    inputCols = ["psfFluxMean"]
    # Output columns are created upon instantiation of the class.
    outputCols = ["psfFluxStetsonJ"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the StetsonJ statistic on the DIA point source fluxes.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}_psfFluxMean".format(band)

        def _stetsonJ(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psfFlux"]),
                                      np.isnan(df["psfFluxErr"]))]
            if len(tmpDf) < 2:
                return np.nan
            fluxes = tmpDf["psfFlux"].to_numpy()
            errors = tmpDf["psfFluxErr"].to_numpy()

            return self._stetson_J(
                fluxes,
                errors,
                diaObjects.at[tmpDf.diaObjectId.iat[0], meanName])

        column = "{}_psfFluxStetsonJ".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.apply(_stetsonJ).astype(dtype)
        else:
            diaObjects.loc[:, column] = filterDiaSources.apply(_stetsonJ)

    def _stetson_J(self, fluxes, errors, mean=None):
        """Compute the single band stetsonJ statistic.

        Parameters
        ----------
        fluxes : `numpy.ndarray` (N,)
            Calibrated lightcurve flux values.
        errors : `numpy.ndarray` (N,)
            Errors on the calibrated lightcurve fluxes.
        mean : `float`
            Starting mean from previous plugin.

        Returns
        -------
        stetsonJ : `float`
            stetsonJ statistic for the input fluxes and errors.

        References
        ----------
        .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
           Parameters for Cepheid Variables", PASP, 108, 851S, 1996
        """
        n_points = len(fluxes)
        flux_mean = self._stetson_mean(fluxes, errors, mean)
        delta_val = (
            np.sqrt(n_points / (n_points - 1)) * (fluxes - flux_mean) / errors)
        p_k = delta_val ** 2 - 1

        return np.mean(np.sign(p_k) * np.sqrt(np.fabs(p_k)))

    def _stetson_mean(self,
                      values,
                      errors,
                      mean=None,
                      alpha=2.,
                      beta=2.,
                      n_iter=20,
                      tol=1e-6):
        """Compute the stetson mean of the fluxes which down-weights outliers.

        Weighted biased on an error weighted difference scaled by a constant
        (1/``a``) and raised to the power beta. Higher betas more harshly
        penalize outliers and ``a`` sets the number of sigma where a weighted
        difference of 1 occurs.

        Parameters
        ----------
        values : `numpy.dnarray`, (N,)
            Input values to compute the mean of.
        errors : `numpy.ndarray`, (N,)
            Errors on the input values.
        mean : `float`
            Starting mean value or None.
        alpha : `float`
            Scalar down-weighting of the fractional difference. lower->more
            clipping. (Default value is 2.)
        beta : `float`
            Power law slope of the used to down-weight outliers. higher->more
            clipping. (Default value is 2.)
        n_iter : `int`
            Number of iterations of clipping.
        tol : `float`
            Fractional and absolute tolerance goal on the change in the mean
            before exiting early. (Default value is 1e-6)

        Returns
        -------
        mean : `float`
            Weighted stetson mean result.

        References
        ----------
        .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
           Parameters for Cepheid Variables", PASP, 108, 851S, 1996
        """
        n_points = len(values)
        n_factor = np.sqrt(n_points / (n_points - 1))
        inv_var = 1 / errors ** 2

        if mean is None:
            mean = np.average(values, weights=inv_var)
        for iter_idx in range(n_iter):
            chi = np.fabs(n_factor * (values - mean) / errors)
            tmp_mean = np.average(
                values,
                weights=inv_var / (1 + (chi / alpha) ** beta))
            diff = np.fabs(tmp_mean - mean)
            mean = tmp_mean
            if diff / mean < tol and diff < tol:
                break
        return mean


class WeightedMeanDiaTotFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanTotFlux")
class WeightedMeanDiaTotFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    forced photometered at the DiaSource location in the calibrated image.

    Additionally store number of usable data points.
    """

    ConfigClass = WeightedMeanDiaPsfFluxConfig
    outputCols = ["scienceFluxMean", "scienceFluxMeanErr"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["invalid value encountered",
                          "divide by zero"])
    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the weighted mean and mean error of the point source flux.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        totMeanName = "{}_scienceFluxMean".format(band)
        if totMeanName not in diaObjects.columns:
            diaObjects[totMeanName] = np.nan
        totErrName = "{}_scienceFluxMeanErr".format(band)
        if totErrName not in diaObjects.columns:
            diaObjects[totErrName] = np.nan

        def _meanFlux(df):
            tmpDf = df[~np.logical_or(np.isnan(df["scienceFlux"]),
                                      np.isnan(df["scienceFluxErr"]))]
            tot_weight = np.nansum(1 / tmpDf["scienceFluxErr"] ** 2)
            fluxMean = np.nansum(tmpDf["scienceFlux"]
                                 / tmpDf["scienceFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)

            return pd.Series({totMeanName: fluxMean,
                              totErrName: fluxMeanErr})

        df = filterDiaSources.apply(_meanFlux).astype(diaObjects.dtypes[[totMeanName, totErrName]])
        diaObjects.loc[:, [totMeanName, totErrName]] = df


class SigmaDiaTotFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaTotFlux")
class SigmaDiaTotFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsfFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["scienceFluxSigma"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """Compute the sigma fluxes of the point source flux measured on the
        calibrated image.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``band``.
        band : `str`
            Simple, string name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        column = "{}_scienceFluxSigma".format(band)
        if column in diaObjects:
            dtype = diaObjects[column].dtype
            diaObjects.loc[:, column] = filterDiaSources.scienceFlux.std().astype(dtype)
        else:

            diaObjects.loc[:, column] = filterDiaSources.scienceFlux.std()
