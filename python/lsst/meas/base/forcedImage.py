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
a specific dataset to be used (see ForcedCcdTask, ForcedCoaddTask).
"""

import lsst.pex.config
import lsst.daf.base
import lsst.pipe.base
import lsst.pex.config

from .base import *
from .references import CoaddSrcReferencesTask
from .forcedMeasurement import *

__all__ = ("ProcessImageForcedConfig", "ProcessImageForcedTask")

class ProcessImageForcedConfig(lsst.pex.config.Config):
    """Config class for forced measurement driver task."""

    references = lsst.pex.config.ConfigurableField(
        target=CoaddSrcReferencesTask,
        doc="Retrieve reference source catalog"
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

class ProcessImageForcedTask(lsst.pipe.base.CmdLineTask):
    ConfigClass = ProcessImageForcedConfig
    _DefaultName = "processImageForcedTask"

    def __init__(self, butler=None, refSchema=None, **kwds):
        super(lsst.pipe.base.CmdLineTask, self).__init__(**kwds)
        self.makeSubtask("references")
        if not refSchema:
            refSchema = self.references.getSchema(butler)
        flags = MeasurementDataFlags()  # just a placeholder for now
        self.makeSubtask("measurement", refSchema=refSchema, flags=flags)


    def run(self, dataRef):
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)
        refCat = list(self.fetchReferences(dataRef, exposure))
        retStruct = self.measurement.run(exposure, refCat, refWcs,
                                         idFactory=self.makeIdFactory(dataRef),
                                         algMetadata=self.measurement.algMetadata)
        self.writeOutput(dataRef, retStruct.sources)

    def makeIdFactory(self, dataRef):
        """Hook for derived classes to define how to make an IdFactory for forced sources.

        Note that this is for forced source IDs, not object IDs, which are usually handled by
        the copyColumns config option.
        """
        raise NotImplementedError()

    def fetchReferences(self, dataRef, exposure):
        """Hook for derived classes to define how to get references objects.

        Derived classes should call one of the fetch* methods on the references subtask,
        but which one they call depends on whether the region to get references for is a
        easy to describe in patches (as it would be when doing forced measurements on a
        coadd), or is just an arbitrary box (as it would be for CCD forced measurements).
        """
        raise NotImplementedError()

    def getExposure(self, dataRef):
        """Read input exposure on which to perform the measurements

        @param dataRef       Data reference from butler
        """
        return dataRef.get(self.dataPrefix + "calexp", immediate=True)

    def writeOutput(self, dataRef, sources):
        """Write forced source table

        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "forced_src")

    def getSchemaCatalogs(self):
        """Get a dict of Schema catalogs that will be used by this Task.
        In the case of forced taks, there is only one schema for each type of forced measurement.
        The dataset type for this measurement is defined in the mapper.
        """
        catalog = lsst.afw.table.SourceCatalog(self.forcedMeasurement.mapper.getOutputSchema())
        catalog.getTable().setMetadata(self.measurement.algMetadata)
        datasetType = self.dataPrefix + "forced"
        return {datasetType:catalog}

    def _getConfigName(self):
        """Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return self.dataPrefix + "forced_config"

    def _getMetadataName(self):
        """Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return self.dataPrefix + "forced_metadata"

