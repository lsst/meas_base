# This file is part of meas_base.
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

import astropy.table
import numpy as np

import lsst.afw.image
import lsst.pex.config
import lsst.pex.exceptions
import lsst.pipe.base
from lsst.utils.logging import PeriodicLogger
from .apCorrRegistry import getApCorrNameSet

# If True then scale instFlux error by apCorr; if False then use a more complex computation
# that over-estimates instFlux error (often grossly so) because it double-counts photon noise.
# This flag is intended to be temporary until we figure out a better way to compute
# the effects of aperture correction on instFlux uncertainty
UseNaiveFluxErr = True

__all__ = ("ApplyApCorrConfig", "ApplyApCorrTask")


class ApCorrInfo:
    """Catalog field names and keys needed to aperture correct a particular
    instrument flux.

    Parameters
    ----------
    schema : `lsst.afw.table`
        Source catalog schema. Three fields are used to generate keys:
        - ``{name}_instFlux``
        - ``{name}_instFluxErr``
        - ``{name}_flag``
        Three fields are added:
        - ``{name}_apCorr`` (only if not already added by proxy)
        - ``{name}_apCorrErr`` (only if not already added by proxy)
        - ``{name}_flag_apCorr``
    model : `str`
        Field name prefix for instFlux with aperture correction model, e.g.
        "base_PsfFlux"
    name : `str`
        Field name prefix for instFlux needing aperture correction; may be
        `None` if it is the same as ``model``

    Notes
    -----
    The aperture correction can be derived from the meaasurements in the
    column being aperture-corrected or from measurements in a different
    column (a "proxy"). In the first case, we will add columns to contain
    the aperture correction values; in the second case (using a proxy),
    we will add an alias to the proxy's aperture correction values. In
    all cases, we add a flag.

    Attributes
    ----------
    name : str
        Field name prefix for flux needing aperture correction.
    modelName : str
        Field name for aperture correction model for flux.
    modelSigmaName : str
        Field name for aperture correction model for fluxErr.
    doApCorrColumn : bool
        Should we write the aperture correction values?
        They should not be written if they're already being written by a proxy.
    instFluxName : str
        Name of ``instFlux`` field.
    instFluxErrName : str
        Name of ``instFlux`` sigma field.
    """

    name = None
    modelName = None
    modelSigmaName = None
    doApCorrColumn = None
    instFluxName = None
    instFluxErrName = None
    apCorrName = None

    def __init__(self, model, schema=None, name=None):
        if name is None:
            name = model
        self.name = name
        self.modelName = model + "_instFlux"
        self.modelSigmaName = model + "_instFluxErr"
        self.instFluxName = name + "_instFlux"
        self.instFluxErrName = name + "_instFluxErr"
        self.fluxFlagName = name + "_flag"
        self.apCorrName = name + "_apCorr"
        self.apCorrErrName = name + "_apCorrErr"
        self.apCorrFlagName = name + "_flag_apCorr"

        # No need to write the same aperture corrections multiple times
        if schema is not None:
            if name == model or model + "_apCorr" not in schema:
                schema.addField(
                    name + "_apCorr",
                    doc="aperture correction applied to %s" % (name,),
                    type=np.float64,
                )
                schema.addField(
                    name + "_apCorrErr",
                    doc="standard deviation of aperture correction applied to %s" % (name,),
                    type=np.float64,
                )
            else:
                aliases = schema.getAliasMap()
                aliases.set(name + "_apCorr", model + "_apCorr")
                aliases.set(name + "_apCorrErr", model + "_apCorrErr")

            if name + "_flag_apCorr" not in schema:
                schema.addField(
                    name + "_flag_apCorr",
                    doc="set if unable to aperture correct %s" % (name,),
                    type="Flag",
                )


class ApplyApCorrConfig(lsst.pex.config.Config):
    """Aperture correction configuration.
    """

    ignoreList = lsst.pex.config.ListField(
        doc="flux measurement algorithms in getApCorrNameSet() to ignore; "
            "if a name is listed that does not appear in getApCorrNameSet() then a warning is logged",
        dtype=str,
        optional=False,
        default=(),
    )
    doFlagApCorrFailures = lsst.pex.config.Field(
        doc="set the general failure flag for a flux when it cannot be aperture-corrected?",
        dtype=bool,
        default=True,
    )
    proxies = lsst.pex.config.DictField(
        doc="flux measurement algorithms to be aperture-corrected by reference to another algorithm; "
            "this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' "
            "is the algorithm supplying the corrections",
        keytype=str,
        itemtype=str,
        default={},
    )
    xColumn = lsst.pex.config.Field(
        doc="name of the x coordinate column in the catalog",
        dtype=str,
        default="slot_Centroid_x",
    )
    yColumn = lsst.pex.config.Field(
        doc="name of the y coordinate column in the catalog",
        dtype=str,
        default="slot_Centroid_y",
    )


class ApplyApCorrTask(lsst.pipe.base.Task):
    """Apply aperture corrections.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
    """
    ConfigClass = ApplyApCorrConfig
    _DefaultName = "applyApCorr"

    def __init__(self, schema=None, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)

        self.apCorrInfoDict = dict()
        apCorrNameSet = getApCorrNameSet()
        ignoreSet = set(self.config.ignoreList)
        missingNameSet = ignoreSet - set(apCorrNameSet)
        if missingNameSet:
            self.log.warning("Fields in ignoreList that are not in fluxCorrectList: %s",
                             sorted(missingNameSet))
        for name in sorted(apCorrNameSet - ignoreSet):
            if schema is not None and name + "_instFlux" not in schema:
                # if a field in the registry is missing from the schema, silently ignore it.
                continue
            self.apCorrInfoDict[name] = ApCorrInfo(schema=schema, model=name)

        for name, model in self.config.proxies.items():
            if name in apCorrNameSet:
                # Already done or ignored
                continue
            if schema is not None and name + "_instFlux" not in schema:
                # Silently ignore
                continue
            self.apCorrInfoDict[name] = ApCorrInfo(schema=schema, model=model, name=name)

    def run(self, catalog, apCorrMap):
        """Apply aperture corrections to a catalog of sources.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog` or `astropy.table.Table`
            Catalog of sources. Will be updated in place.
        apCorrMap : `lsst.afw.image.ApCorrMap`
            Aperture correction map

        Notes
        -----
        If you show debug-level log messages then you will see statistics for
        the effects of aperture correction.
        """
        self.log.info("Applying aperture corrections to %d instFlux fields", len(self.apCorrInfoDict))
        if UseNaiveFluxErr:
            self.log.debug("Use naive instFlux sigma computation")
        else:
            self.log.debug("Use complex instFlux sigma computation that double-counts photon noise "
                           "and thus over-estimates instFlux uncertainty")

        # Wrap the task logger to a periodic logger.
        periodicLog = PeriodicLogger(self.log)

        if isinstance(catalog, astropy.table.Table):
            # Prepare the astropy table for use with this task.
            self._prepAstropyTable(catalog)
        xColumn = self.config.xColumn
        yColumn = self.config.yColumn

        for apCorrIndex, apCorrInfo in enumerate(self.apCorrInfoDict.values()):
            apCorrModel = apCorrMap.get(apCorrInfo.modelName)
            apCorrErrModel = apCorrMap.get(apCorrInfo.modelSigmaName)
            if None in (apCorrModel, apCorrErrModel):
                missingNames = [(apCorrInfo.modelName, apCorrInfo.modelSigmaName)[i]
                                for i, model in enumerate((apCorrModel, apCorrErrModel)) if model is None]
                self.log.warning("Cannot aperture correct %s because could not find %s in apCorrMap",
                                 apCorrInfo.name, " or ".join(missingNames))
                catalog[apCorrInfo.apCorrFlagName] = np.ones(len(catalog), dtype=np.bool_)
                continue

            # Say we've failed before we start; we'll unset these flags on success.
            catalog[apCorrInfo.apCorrFlagName] = True
            oldFluxFlagState = np.zeros(len(catalog), dtype=np.bool_)
            if self.config.doFlagApCorrFailures:
                oldFluxFlagState = catalog[apCorrInfo.fluxFlagName]
                catalog[apCorrInfo.fluxFlagName] = True

            apCorr = apCorrModel.evaluate(catalog[xColumn], catalog[yColumn])
            if not UseNaiveFluxErr:
                apCorrErr = apCorrErrModel.evaluate(
                    catalog[xColumn],
                    catalog[yColumn],
                )
            else:
                apCorrErr = np.zeros(len(catalog))

            if apCorrInfo.doApCorrColumn:
                catalog[apCorrInfo.apCorrName] = apCorr
                catalog[apCorrInfo.apCorrErrName] = apCorrErr

            # Check bad values that are not finite (possible for coadds).
            good = np.isfinite(apCorr) & (apCorr > 0.0) & (apCorrErr >= 0.0)

            if good.sum() == 0:
                continue

            instFlux = catalog[apCorrInfo.instFluxName]
            instFluxErr = catalog[apCorrInfo.instFluxErrName]
            catalog[apCorrInfo.instFluxName][good] = instFlux[good]*apCorr[good]
            if UseNaiveFluxErr:
                catalog[apCorrInfo.instFluxErrName][good] = instFluxErr[good]*apCorr[good]
            else:
                a = instFluxErr[good]/instFlux[good]
                b = apCorrErr[good]/apCorr[good]
                err = np.abs(instFlux[good]*apCorr[good])*np.sqrt(a*a + b*b)
                catalog[apCorrInfo.instFluxErrName][good] = err

            flags = catalog[apCorrInfo.apCorrFlagName].copy()
            flags[good] = False
            catalog[apCorrInfo.apCorrFlagName] = flags
            if self.config.doFlagApCorrFailures:
                fluxFlags = catalog[apCorrInfo.fluxFlagName].copy()
                fluxFlags[good] = oldFluxFlagState[good]
                catalog[apCorrInfo.fluxFlagName] = fluxFlags

            # Log a message if it has been a while since the last log.
            periodicLog.log("Aperture corrections applied to %d fields out of %d",
                            apCorrIndex + 1, len(self.apCorrInfoDict))

            if self.log.isEnabledFor(self.log.DEBUG):
                # log statistics on the effects of aperture correction
                apCorrArr = np.asarray([s[apCorrInfo.instFluxName] for s in catalog])
                apCorrErrArr = np.asarray([s[apCorrInfo.instFluxErrName] for s in catalog])
                self.log.debug("For instFlux field %r: mean apCorr=%s, stdDev apCorr=%s, "
                               "mean apCorrErr=%s, stdDev apCorrErr=%s for %s sources",
                               apCorrInfo.name, apCorrArr.mean(), apCorrArr.std(),
                               apCorrErrArr.mean(), apCorrErrArr.std(), len(catalog))

    def _prepAstropyTable(self, table):
        """Prepare an astropy table for use with this task.

        Parameters
        ----------
        table : `astropy.table.Table`
            Table to prepare.

        Returns
        -------
        `astropy.table.Table`
            Prepared table.
        """
        for apCorrInfo in self.apCorrInfoDict.values():
            if apCorrInfo.apCorrName not in table.colnames:
                table[apCorrInfo.apCorrName] = np.zeros(len(table), dtype=np.float64)
            if apCorrInfo.apCorrErrName not in table.colnames:
                table[apCorrInfo.apCorrErrName] = np.zeros(len(table), dtype=np.float64)
            if apCorrInfo.apCorrFlagName not in table.colnames:
                table[apCorrInfo.apCorrFlagName] = np.zeros(len(table), dtype=bool)
            if apCorrInfo.fluxFlagName not in table.colnames:
                table[apCorrInfo.fluxFlagName] = np.zeros(len(table), dtype=bool)
