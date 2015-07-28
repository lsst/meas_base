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

"""Base command-line driver task for forced measurement.  Must be inherited to specialize for
a specific dataset to be used (see ForcedPhotCcdTask, ForcedPhotCoaddTask).
"""

import lsst.pex.config
import lsst.daf.base
import lsst.pipe.base
import lsst.pex.config

from .references import MultiBandReferencesTask
from .forcedMeasurement import ForcedMeasurementTask

__all__ = ("ProcessImageForcedConfig", "ProcessImageForcedTask")

class ProcessImageForcedConfig(lsst.pex.config.Config):
    """!Config class for forced measurement driver task."""

    references = lsst.pex.config.ConfigurableField(
        target=MultiBandReferencesTask,
        doc="subtask to retrieve reference source catalog"
        )
    measurement = lsst.pex.config.ConfigurableField(
        target=ForcedMeasurementTask,
        doc="subtask to do forced measurement"
        )
    coaddName = lsst.pex.config.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
        )
    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent": "parentObjectId", "deblend_nChild": "deblend_nChild"}
        )

    def setDefaults(self):
        # coadds do not yet include include aperture correction data,
        # but that is planned, so might as well warn when it occurs
        self.measurement.doApplyApCorr = "noButWarn"

## @addtogroup LSST_task_documentation
## @{
## @page ProcessImageForcedTask
## ProcessImageForcedTask
## @copybrief ProcessImageForcedTask
## @}

class ProcessImageForcedTask(lsst.pipe.base.CmdLineTask):
    """!A base class for command-line forced measurement drivers.

    This is a an abstract class, which is the common ancestor for ForcedPhotCcdTask
    and ForcedPhotCoaddTask.  It provides the run() method that does most of the
    work, while delegating a few customization tasks to other methods that are
    overridden by subclasses.

    This task is not directly usable as a CmdLineTask; subclasses must:
     - Set the _DefaultName class attribute
     - Implement makeIdFactory
     - Implement fetchReferences
     - (optional) Implement attachFootprints
    """
    ConfigClass = ProcessImageForcedConfig
    _DefaultName = "processImageForcedTask"

    def __init__(self, butler=None, refSchema=None, **kwds):
        """Initialize the task.

        ForcedPhotImageTask takes two keyword arguments beyond the usual CmdLineTask arguments:
         - refSchema: the Schema of the reference catalog, passed to the constructor of the references
           subtask
         - butler: a butler that will be passed to the references subtask to allow it to load its Schema
           from disk
        At least one of these arguments must be present; if both are, schema takes precedence.
        """
        super(lsst.pipe.base.CmdLineTask, self).__init__(**kwds)
        self.makeSubtask("references", butler=butler, schema=refSchema)
        if refSchema is None:
            refSchema = self.references.schema
        self.makeSubtask("measurement", refSchema=refSchema)

    def run(self, dataRef):
        """!Measure a single exposure using forced detection for a reference catalog.

        @param[in]  dataRef   An lsst.daf.persistence.ButlerDataRef. It is passed to the
                              references subtask to obtain the reference WCS, the getExposure()
                              method (implemented by derived classes) to read the measurement
                              image, and the fetchReferences() method (implemented by derived
                              classes) to get the exposure and load the reference catalog (see
                              the CoaddSrcReferencesTask for more information).  Sources are
                              generated with generateSources() in the measurement subtask.  These
                              are passed to measurement's run method which fills the source
                              catalog with the forced measurement results.  The sources are then
                              passed to the writeOutputs() method (implemented by derived classes)
                              which writes the outputs.  See derived class documentation for which
                              datasets and data ID keys are used.
        """
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)
        refCat = self.fetchReferences(dataRef, exposure)
        sources = self.measurement.generateSources(exposure, refCat, refWcs,
                                                   idFactory=self.makeIdFactory(dataRef))
        self.log.info("Performing forced measurement on %s" % dataRef.dataId)
        self.attachFootprints(dataRef, sources, refCat, exposure, refWcs)
        self.measurement.run(exposure, sources, refCat, refWcs)
        self.writeOutput(dataRef, sources)

    def makeIdFactory(self, dataRef):
        """!Hook for derived classes to define how to make an IdFactory for forced sources.

        Note that this is for forced source IDs, not object IDs, which are usually handled by
        the copyColumns config option.
        """
        raise NotImplementedError()

    def fetchReferences(self, dataRef, exposure):
        """!Hook for derived classes to define how to get references objects.

        Derived classes should call one of the fetch* methods on the references subtask,
        but which one they call depends on whether the region to get references for is a
        easy to describe in patches (as it would be when doing forced measurements on a
        coadd), or is just an arbitrary box (as it would be for CCD forced measurements).
        """
        raise NotImplementedError()

    def attachFootprints(self, dataRef, sources, refCat, exposure, refWcs):
        """!Hook for derived classes to define how to attach Footprints to blank sources prior to measurement

        Footprints for forced photometry must be in the pixel coordinate system of the image being
        measured, while the actual detections may start out in a different coordinate system.

        Subclasses for ForcedPhotImageTask must implement this method to define how those Footprints
        should be generated.

        The default implementation transforms the Footprints from the reference catalog from the refWcs
        to the exposure's Wcs, which downgrades HeavyFootprints into regular Footprints, destroying
        deblend information.
        """
        exposureWcs = exposure.getWcs()
        region = exposure.getBBox(lsst.afw.image.PARENT)
        for srcRecord, refRecord in zip(sources, refCat):
            srcRecord.setFootprint(refRecord.getFootprint().transform(refWcs, exposureWcs, region))

    def getExposure(self, dataRef):
        """!Read input exposure on which to perform the measurements

        @param dataRef       Data reference from butler.
        """
        return dataRef.get(self.dataPrefix + "calexp", immediate=True)

    def writeOutput(self, dataRef, sources):
        """!Write forced source table

        @param dataRef  Data reference from butler; the forced_src dataset (with self.dataPrefix included)
                        is all that will be modified.
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "forced_src")

    def getSchemaCatalogs(self):
        """!Get a dict of Schema catalogs that will be used by this Task.

        In the case of forced taks, there is only one schema for each type of forced measurement.
        The dataset type for this measurement is defined in the mapper.
        """
        catalog = lsst.afw.table.SourceCatalog(self.measurement.mapper.getOutputSchema())
        catalog.getTable().setMetadata(self.measurement.algMetadata)
        datasetType = self.dataPrefix + "forced_src"
        return {datasetType:catalog}

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return self.dataPrefix + "forced_config"

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return self.dataPrefix + "forced_metadata"

