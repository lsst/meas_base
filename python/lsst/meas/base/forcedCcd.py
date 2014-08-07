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

"""Class to measure a CCD calexp exposure, using a reference catalog from a Coadd.
"""

import argparse

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.table

from .forcedImage import *
from .base import *

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("ProcessForcedCcdConfig", "ProcessForcedCcdTask")

class ProcessForcedCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A version of DataIdContainer specialized for forced photometry on CCDs.

    Required because we need to add "tract" to the raw data ID keys, and that's tricky.
    This IdContainer assumes that a calexp is being measured using the detection information
    from the set of coadds which intersect with the calexp.  a set of reference catalog
    from a coadd which overlaps it.  It needs the calexp id (e.g. visit, raft, sensor), but it
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
    doApplyUberCal = lsst.pex.config.Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = True
    )

## \addtogroup LSST_task_documentation
## \{
## \page processForcedCcdTask
## \ref ProcessForcedCcdTask_ "ProcessForcedCcdTask"
## \copybrief ProcessForcedCcdTask
## \}

class ProcessForcedCcdTask(ProcessImageForcedTask):
    """!
    \anchor ProcessForcedCcdTask_
    \brief The ProcessForcedCcdTask is a command-line driver for measuring the properties of sources on a
           single ccd.

    \section meas_base_processForcedCcdTask_Contents Contents

     - \ref meas_base_processForcedCcdTask_Purpose
     - \ref meas_base_processForcedCcdTask_IO
     - \ref meas_base_processForcedCcdTask_Config

    \section meas_base_processForcedCcdTask_Purpose	Description

    \copybrief ProcessForcedCcdTask

    This task is a subclass of ProcessForcedImageTask which is specifically for doing forced
    measurement on a single exposure, using as a reference catalog the detections which
    were made on overlapping coadds.

    \section meas_base_processForcedCcdTask_IO		Inputs/Outputs to the run method

    See \link lsst.meas.base.forcedImage.ProcessImageForcedTask#run ProcessImageForcedTask.run\endlink
    for more information.

    The run method (inherited from ProcessForcedImageTask) takes a lsst.daf.persistence.ButlerDataRef
    argument that corresponds to a single CCD.  This should contain the data ID keys that correspond to
    the "forced_src" dataset (the output dataset for ProcessForcedCcdTask), which are typically all those
    used to specify the "calexp" dataset (e.g. visit, raft, sensor for LSST data) as well as a coadd
    tract.  The tract is used to look up the appropriate coadd measurement catalogs to use as references
    (e.g. deepCoadd_src; see CoaddSrcReferencesTask for more information). While the tract must be given
    as part of the dataRef, the patches are determined automatically from the bounding box and WCS of the
    calexp to be measured.

    In addition to the run method, ProcessForcedCcdTask overrides several methods of ProcessForcedImageTask
    to specialize it for single-CCD processing, including makeIdFactory(), fetchReferences(), and
    getExposure().  None of these should be called directly by the user, though it may be useful
    to override them further in subclasses.

    \section meas_base_processForcedCcdTask_Config       Configuration parameters

    See \ref ProcessForcedCcdConfig
    """

    ConfigClass = ProcessForcedCcdConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
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
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=ProcessForcedCcdDataIdContainer)
        return parser


