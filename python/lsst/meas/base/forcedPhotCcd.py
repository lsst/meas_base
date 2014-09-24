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

import argparse

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.table

from .forcedPhotImage import ProcessImageForcedTask, ProcessImageForcedConfig

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("ForcedPhotCcdConfig", "ForcedPhotCcdTask")

class ForcedPhotCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
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

class ForcedPhotCcdConfig(ProcessImageForcedConfig):
    doApplyUberCal = lsst.pex.config.Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = False
    )

## @addtogroup LSST_task_documentation
## @{
## @page processForcedCcdTask
## ForcedPhotCcdTask
## @copybrief ForcedPhotCcdTask
## @}

class ForcedPhotCcdTask(ProcessImageForcedTask):
    """!
    A command-line driver for performing forced measurement on CCD images

    This task is a subclass of ForcedPhotImageTask which is specifically for doing forced
    measurement on a single CCD exposure, using as a reference catalog the detections which
    were made on overlapping coadds.

    The run method (inherited from ForcedPhotImageTask) takes a lsst.daf.persistence.ButlerDataRef
    argument that corresponds to a single CCD.  This should contain the data ID keys that correspond to
    the "forced_src" dataset (the output dataset for ForcedPhotCcdTask), which are typically all those
    used to specify the "calexp" dataset (e.g. visit, raft, sensor for LSST data) as well as a coadd
    tract.  The tract is used to look up the appropriate coadd measurement catalogs to use as references
    (e.g. deepCoadd_src; see CoaddSrcReferencesTask for more information). While the tract must be given
    as part of the dataRef, the patches are determined automatically from the bounding box and WCS of the
    calexp to be measured, and the filter used to fetch references is set via config
    (BaseReferencesConfig.filter).

    In addition to the run method, ForcedPhotCcdTask overrides several methods of ForcedPhotImageTask
    to specialize it for single-CCD processing, including makeIdFactory(), fetchReferences(), and
    getExposure().  None of these should be called directly by the user, though it may be useful
    to override them further in subclasses.
    """

    ConfigClass = ForcedPhotCcdConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcdTask"
    dataPrefix = ""

    def makeIdFactory(self, dataRef):
        """Create an object that generates globally unique source IDs from per-CCD IDs and the CCD ID.

        @param dataRef       Data reference from butler.  The "ccdExposureId_bits" and "ccdExposureId"
                             datasets are accessed.  The data ID must have the keys that correspond
                             to ccdExposureId, which is generally the same that correspond to "calexp"
                             (e.g. visit, raft, sensor for LSST data).
        """
        expBits = dataRef.get("ccdExposureId_bits")
        expId = long(dataRef.get("ccdExposureId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def fetchReferences(self, dataRef, exposure):
        """Return an iterable of reference sources which overlap the exposure

        @param dataRef       Data reference from butler corresponding to the image to be measured;
                             should have tract, patch, and filter keys.
        @param exposure      lsst.afw.image.Exposure to be measured (used only to obtain a Wcs and bounding
                             box).

        All work is delegated to the references subtask; see CoaddSrcReferencesTask for information
        about the default behavior.
        """
        bbox = exposure.getBBox()
        return self.references.fetchInBox(dataRef, bbox, exposure.getWcs())

    def getExposure(self, dataRef):
        """Read input exposure to measure

        @param dataRef       Data reference from butler.  Only the 'calexp' dataset is used,
                             unless config.doApplyUberCal is true, in which case the corresponding
                             meas_mosaic outputs are used as well.
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
                               ContainerClass=ForcedPhotCcdDataIdContainer)
        return parser


