#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2014 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""Class to measure a Ccd calexp exposure, using a reference catalog from a Coadd
   See the comment in the parent class in forced.py for more information.
"""
import argparse

from lsst.pex.config import Config, ConfigurableField, DictField, Field
from lsst.pipe.base import ArgumentParser,ButlerInitializedTaskRunner,DataIdContainer
import lsst.afw.image

from .forcedImage import *
from .base import *

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("ProcessForcedCcdConfig", "ProcessForcedCcdTask")

class ProcessForcedCcdDataIdContainer(DataIdContainer):
    """A version of DataIdContainer specialized for forced photometry on CCDs.

    Required because we need to add "tract" to the raw data ID keys, and that's tricky.
    This IdContainer assumes that a calexp is being measured using the detection information
    from the set of coadds which intersect with the calexp.  a set of reference catalog
    from a coadd which overlaps it.  It needs the calexp id (visit, raft, sensor), but it
    also uses the tract to decide what set of coadds to use.  The references from the tract
    whose patches intersect with the calexp are used.
    """
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        for dataId in self.idList:
            if "tract" not in dataId:
                raise argparse.ArgumentError(None, "--id must include tract")
            tract = dataId.pop("tract")
            # making a DataRef for src fills out any missing keys and allows us to iterate
            for srcDataRef in namespace.butler.subset("src", dataId=dataId):
                forcedDataId = srcDataRef.dataId.copy()
                forcedDataId['tract'] = tract
                dataRef = namespace.butler.dataRef(
                    datasetType = "forced_src",
                    dataId = forcedDataId,
                    )
                self.refList.append(dataRef)

class ProcessForcedCcdConfig(ProcessImageForcedConfig):
    doApplyUberCal = Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = False
    )


class ProcessForcedCcdTask(ProcessImageForcedTask):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclasses for running on Ccds).
    """

    ConfigClass = ProcessForcedCcdConfig
    RunnerClass = ButlerInitializedTaskRunner
    _DefaultName = "forcedCcdTask"
    dataPrefix = ""

    def makeIdFactory(self, dataRef):
        """For the output table, make a source id for each output from the # exposure and ccd ids

        @param dataRef       Data reference from butler
        """
        expBits = dataRef.get("ccdExposureId_bits")
        expId = long(dataRef.get("ccdExposureId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def fetchReferences(self, dataRef, exposure):
        """Return an iterable of reference sources which overlap the exposure

        @param dataRef       Data reference from butler
        @param exposure      afw Exposure
        """
        bbox = exposure.getBBox(lsst.afw.image.PARENT)
        return self.references.fetchInBox(dataRef, bbox, exposure.getWcs())

    def getExposure(self, dataRef):
        """Read input exposure to measure

        @param dataRef       Data reference from butler
        """
        exposure = ProcessImageForcedTask.getExposure(self, dataRef)
        if not self.config.doApplyUberCal:
            return exposure
        if applyMosaicResults is None:
            raise RuntimeError(
                "Cannot use improved calibrations for %s because meas_mosaic could not be imported."
                % dataRef.dataId
                )
        else:
            applyMosaicResults(dataRef, calexp=exposure)
        return exposure

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=ProcessForcedCcdDataIdContainer)
        return parser


