from __future__ import absolute_import, division
#
# LSST Data Management System
# Copyright 2015 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math

import numpy

import lsst.pex.config
import lsst.pex.exceptions
import lsst.afw.image
import lsst.pipe.base
from .apCorrRegistry import getApCorrNameSet

# If True then scale flux sigma by apCorr; if False then use a more complex computation
# that over-estimates flux error (often grossly so) because it double-counts photon noise.
# This flag is intended to be temporary until we figure out a better way to compute
# the effects of aperture correction on flux uncertainty
UseNaiveFluxSigma = True

__all__ = ("ApplyApCorrConfig", "ApplyApCorrTask")


class ApCorrInfo(object):
    """!Catalog field names and keys needed to aperture correct a particular flux
    """

    def __init__(self, schema, name):
        """!Construct an ApCorrInfo and add fields to the schema

        @param[in,out] schema  source catalog schema;
            three fields are used to generate keys:
            - {name}_flux
            - {name}_fluxSigma
            - {name}_flag
            three fields are added:
            - {name}_apCorr
            - {name}_apCorrSigma
            - {name}_flag_apCorr
        @param[in] name  field name prefix for flux needing aperture correction, e.g. "base_PsfFlux"

        ApCorrInfo has the following attributes:
        - name: field name prefix for flux needing aperture correction
        - fluxName: name of flux field
        - fluxSigmaName: name of flux sigma field
        - fluxKey: key to flux field
        - fluxSigmaKey: key to flux sigma field
        - fluxFlagKey: key to flux flag field
        - apCorrKey: key to new aperture correction field
        - apCorrSigmaKey: key to new aperture correction sigma field
        - apCorrFlagKey: key to new aperture correction flag field
        """
        self.name = name
        self.fluxName = name + "_flux"
        self.fluxSigmaName = name + "_fluxSigma"
        self.fluxKey = schema.find(self.fluxName).key
        self.fluxSigmaKey = schema.find(self.fluxSigmaName).key
        self.fluxFlagKey = schema.find(name + "_flag").key
        self.apCorrKey = schema.addField(
            name + "_apCorr",
            doc="aperture correction applied to %s" % (name,),
            type=float,
        )
        self.apCorrSigmaKey = schema.addField(
            name + "_apCorrSigma",
            doc="aperture correction applied to %s" % (name,),
            type=float,
        )
        self.apCorrFlagKey = schema.addField(
            name + "_flag_apCorr",
            doc="set if unable to aperture correct %s" % (name,),
            type="Flag",
        )


class ApplyApCorrConfig(lsst.pex.config.Config):
    ignoreList = lsst.pex.config.ListField(
        doc="flux measurement algorithms in getApCorrNameSet() to ignore;" +
        " if a name is listed that does not appear in getApCorrNameSet() then a warning is logged",
        dtype=str,
        optional=False,
        default=(),
    )
    doFlagApCorrFailures = lsst.pex.config.Field(
        doc="set the general failure flag for a flux when it cannot be aperture-corrected?",
        dtype=bool,
        default=True,
    )


class ApplyApCorrTask(lsst.pipe.base.Task):
    """!Apply aperture corrections
    """
    ConfigClass = ApplyApCorrConfig
    _DefaultName = "applyApCorr"

    def __init__(self, schema, **kwds):
        """Construct an instance of this task
        """
        lsst.pipe.base.Task.__init__(self, **kwds)

        self.apCorrInfoDict = dict()
        apCorrNameSet = getApCorrNameSet()
        ignoreSet = set(self.config.ignoreList)
        missingNameSet = ignoreSet - set(apCorrNameSet)
        if missingNameSet:
            self.log.warn("Fields in ignoreList that are not in fluxCorrectList: %s" %
                          (sorted(list(missingNameSet)),))
        for name in apCorrNameSet - ignoreSet:
            if name + "_flux" in schema:
                self.apCorrInfoDict[name] = ApCorrInfo(schema=schema, name=name)
            else:
                # if a field in the registry is missing from the schema, silently ignore it.
                pass

    def run(self, catalog, apCorrMap):
        """Apply aperture corrections to a catalog of sources

        @param[in,out] catalog  catalog of sources
        @param[in] apCorrMap  aperture correction map (an lsst.afw.image.ApCorrMap)

        If you show debug-level log messages then you will see statistics for the effects of
        aperture correction.
        """
        self.log.info("Applying aperture corrections to %d flux fields" % (len(self.apCorrInfoDict),))
        if UseNaiveFluxSigma:
            self.log.info("Use naive flux sigma computation")
        else:
            self.log.info("Use complex flux sigma computation that double-counts photon noise "
                          " and thus over-estimates flux uncertainty")
        for apCorrInfo in self.apCorrInfoDict.itervalues():
            apCorrModel = apCorrMap.get(apCorrInfo.fluxName)
            apCorrSigmaModel = apCorrMap.get(apCorrInfo.fluxSigmaName)
            if None in (apCorrModel, apCorrSigmaModel):
                missingNames = [(apCorrInfo.fluxName, apCorrInfo.fluxSigmaName)[i]
                                for i, model in enumerate((apCorrModel, apCorrSigmaModel)) if model is None]
                self.log.warn("Could not find %s in apCorrMap" % (" or ".join(missingNames),))
                for source in catalog:
                    source.set(apCorrInfo.apCorrFlagKey, True)
                continue

            for source in catalog:
                center = source.getCentroid()
                # say we've failed when we start; we'll unset these flags when we succeed
                source.set(apCorrInfo.apCorrFlagKey, True)
                oldFluxFlagState = False
                if self.config.doFlagApCorrFailures:
                    oldFluxFlagState = source.get(apCorrInfo.fluxFlagKey)
                    source.set(apCorrInfo.fluxFlagKey, True)

                apCorr = 1.0
                apCorrSigma = 0.0
                try:
                    apCorr = apCorrModel.evaluate(center)
                    if not UseNaiveFluxSigma:
                        apCorrSigma = apCorrSigmaModel.evaluate(center)
                except lsst.pex.exceptions.DomainError:
                    continue

                source.set(apCorrInfo.apCorrKey, apCorr)
                source.set(apCorrInfo.apCorrSigmaKey, apCorrSigma)
                if apCorr <= 0.0 or apCorrSigma < 0.0:
                    continue

                flux = source.get(apCorrInfo.fluxKey)
                fluxSigma = source.get(apCorrInfo.fluxSigmaKey)
                source.set(apCorrInfo.fluxKey, flux*apCorr)
                if UseNaiveFluxSigma:
                    source.set(apCorrInfo.fluxSigmaKey, fluxSigma*apCorr)
                else:
                    a = fluxSigma/flux
                    b = apCorrSigma/apCorr
                    source.set(apCorrInfo.fluxSigmaKey, abs(flux*apCorr)*math.sqrt(a*a + b*b))
                source.set(apCorrInfo.apCorrFlagKey, False)
                if self.config.doFlagApCorrFailures:
                    source.set(apCorrInfo.fluxFlagKey, oldFluxFlagState)

            if self.log.getThreshold() <= self.log.DEBUG:
                # log statistics on the effects of aperture correction
                apCorrArr = numpy.array([s.get(apCorrInfo.apCorrKey) for s in catalog])
                apCorrSigmaArr = numpy.array([s.get(apCorrInfo.apCorrSigmaKey) for s in catalog])
                self.log.logdebug("For flux field %r: mean apCorr=%s, stdDev apCorr=%s,"
                                  " mean apCorrSigma=%s, stdDev apCorrSigma=%s for %s sources" %
                                  (apCorrInfo.name, apCorrArr.mean(), apCorrArr.std(),
                                   apCorrSigmaArr.mean(), apCorrSigmaArr.std(), len(catalog)))
