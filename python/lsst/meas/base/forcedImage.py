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
"""Base classes for forced measurement plugins and the driver task for these.

In forced measurement, a reference catalog is used to define restricted measurements (usually just fluxes)
on an image.  As the reference catalog may be deeper than the detection limit of the measurement image, we
do not assume that we can use detection and deblend information from the measurement image.  Instead, we
assume this information is present in the reference catalog and can be "transformed" in some sense to
the measurement frame.  At the very least, this means that Footprints from the reference catalog should
be transformed and installed as Footprints in the output measurement catalog.  If we have a procedure that
can transform HeavyFootprints, we can then proceed with measurement as usual, but using the reference
catalog's id and parent fields to define deblend families.  If this transformation does not preserve
HeavyFootprints (this is currently the case), then we will only be able to replace objects with noise
one deblend family at a time, and hence measurements run in single-object mode may be contaminated by
neighbors when run on objects with parent != 0.

Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).  Note that
the reference catalog may be in a different coordinate system; it is the responsibility of plugins
to transform the data they need themselves, using the reference WCS provided.  However, for plugins
that only require a position or shape, they may simply use output SourceCatalog's centroid or shape slots,
which will generally be set to the transformed position of the reference object before any other plugins are
run, and hence avoid using the reference catalog at all.
"""

import math

import lsst.pex.config
import lsst.daf.base
from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
from lsst.pex.config import DictField, ConfigurableField

from .base import *
from .references import CoaddSrcReferencesTask, BaseReferencesTask

__all__ = ("ForcedPluginConfig", "ForcedPlugin", "ForcedMeasurementConfig", "ForcedMeasurementTask")

class ForcedPluginConfig(BasePluginConfig):
    """Base class for configs of forced measurement plugins."""
    pass

class ForcedPlugin(BasePlugin):

    # All subclasses of ForcedPlugin should be registered here
    registry = PluginRegistry(ForcedPluginConfig)

    ConfigClass = ForcedPluginConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        @param[in,out]  schemaMapper  A SchemaMapper that maps reference catalog fields to output
                                      catalog fields.  Output fields should be added to the
                                      output schema.  While most plugins will not need to map
                                      fields from the reference schema, if they do so, those fields
                                      will be transferred before any plugins are run.
        @param[in]  flags        A set of bitflags describing the data that the plugin
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       A PluginMap of previously-initialized plugins
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        self.config = config
        self.name = name

    def measure(self, measRecord, exposure, refRecord, refWcs):
        """Measure the properties of a source on a single image, given data from a
        reference record.

        @param[in] exposure       lsst.afw.image.ExposureF, containing the pixel data to
                                  be measured and the associated Psf, Wcs, etc.  All
                                  other sources in the image will have been replaced by
                                  noise according to deblender outputs.
        @param[in,out] measRecord lsst.afw.table.SourceRecord to be filled with outputs,
                                  and from which previously-measured quantities can be
                                  retreived.
        @param[in] refRecord      lsst.afw.table.SimpleRecord that contains additional
                                  parameters to define the fit, as measured elsewhere.
        @param[in] refWcs         The coordinate system for the reference catalog values.
                                  An afw.geom.Angle may be passed, indicating that a
                                  local tangent Wcs should be created for each object
                                  using afw.image.makeLocalWcs and the given angle as
                                  a pixel scale.

        In the normal mode of operation, the source centroid will be set to the
        WCS-transformed position of the reference object, so plugins that only
        require a reference position should not have to access the reference object
        at all.
        """
        raise NotImplementedError()

    def measureN(self, measCat, exposure, refCat, refWcs):
        """Measure the properties of a group of blended sources on a single image,
        given data from a reference record.

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  Sources
                                 not in the blended hierarchy to be measured will have
                                 been replaced with noise using deblender outputs.
        @param[in,out] measCat   lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be
                                 retrieved, containing only the sources that should be
                                 measured together in this call.
        @param[in] refCat        lsst.afw.table.SimpleCatalog that contains additional
                                 parameters to define the fit, as measured elsewhere.
                                 Ordered such that zip(sources, references) may be used.
        @param[in] refWcs        The coordinate system for the reference catalog values.
                                 An afw.geom.Angle may be passed, indicating that a
                                 local tangent Wcs should be created for each object
                                 using afw.image.makeLocalWcs and the given Angle as
                                 a pixel scale.

        In the normal mode of operation, the source centroids will be set to the
        WCS-transformed position of the reference object, so plugins that only
        require a reference position should not have to access the reference object
        at all.
        """
        raise NotImplementedError()

class ForcedMeasurementConfig(BaseMeasurementConfig):
    """Config class for forced measurement driver task."""

    plugins = ForcedPlugin.registry.makeField(
        multi=True,
        default=["centroid.peak",
                 ],
        doc="Plugins to be run and their configuration"
        )
    references = ConfigurableField(
        target=CoaddSrcReferencesTask,
        doc="Retrieve reference source catalog"
        )
    copyColumns = DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent":"parentObjectId"}
        )
    coaddName = lsst.pex.config.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
    )

    def setDefaults(self):
        self.slots.shape = None
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = None
        self.slots.instFlux = None

class ForcedMeasurementTask(CmdLineTask):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclassed for running on CCDs and Coadds).
    """

    ConfigClass = ForcedMeasurementConfig
    RunnerClass = ButlerInitializedTaskRunner
    # this must be defined by the a child class: dataPrefix =
    _DefaultName = "forcedMeasurementTask"

    #  The primary input to this init is the butler, which is a keyword argument,
    #
    def __init__(self, butler=None, refSchema=None, **kwds):
        CmdLineTask.__init__(self, **kwds)
        self.makeSubtask("references")
        self.algMetadata = lsst.daf.base.PropertyList()
        self.plugins = PluginMap()
        if not refSchema:
            refSchema = self.references.getSchema(butler)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        minimalSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.mapper.addMinimalSchema(minimalSchema)
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)

        flags = MeasurementDataFlags()
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, self.mapper, flags=flags,
                                                   others=self.plugins, metadata=self.algMetadata)

    def run(self, dataRef):
        """Main driver method called once for each measurement image by the command-line interface.
        """
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)
        refCat = list(self.fetchReferences(dataRef, exposure))
        retStruct = self.forcedMeasure(exposure, refCat, refWcs, dataRef=dataRef)
        self.writeOutput(dataRef, retStruct.sources)

    def forcedMeasure(self, exposure, refCat, refWcs, dataRef=None):
        """This is a service routine which allows a forced measurement to be made without
        constructing the reference list.  The reference list and refWcs can then be passed in.
        """
        # First create a refList from the original which excludes children when a member
        # of the parent chain is not within the list.  This can occur at boundaries when
        # the parent is outside and one of the children is within.  Currently, the parent
        # chain is always only one deep, but just in case, this code removes any reference
        # when the parent chain to its topmost parent is broken

        # Construct a footprints dict which does not contain the footprint, just the parentID
        # We will add the footprint from the transformed sources from generateSources later.
        footprints = {ref.getId(): ref.getParent() for ref in refCat}
        refList = list()
        for ref in refCat:
            refId = ref.getId()
            topId = refId
            while(topId > 0):
                if not topId in footprints.keys():
                    footprints.pop(refId)
                    break
                topId = footprints[topId]
            if topId == 0:
                refList.append(ref)

        # now generate transformed source corresponding to the cleanup up refLst
        sources = self.generateSources(dataRef, exposure, refList, refWcs)

        # Steal the transformed source footprint and use it to complete the footprints dict,
        # which then looks like {ref.getId(): (ref.getParent(), source.getFootprint())}
        for (reference, source) in zip(refList, sources):
            footprints[reference.getId()] = (reference.getParent(), source.getFootprint())

        self.log.info("Performing forced measurement on %d sources" % len(refList))

        # Build a catalog of just the references we intend to measure
        referenceCat = lsst.afw.table.SourceCatalog(self.mapper.getInputSchema())
        referenceCat.extend(refList)

        self.config.slots.setupTable(sources.table)

        # convert the footprints to the coordinate system of the exposure
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)

        # Create parent cat which slices both the referenceCat and measCat (sources)
        # first, get the reference and source records which have no parent
        refParentCat, measParentCat = referenceCat.getChildren(0, sources)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat,measParentCat)):

            # first process the records which have the current parent as children
            refChildCat, measChildCat = referenceCat.getChildren(refParentRecord.getId(), sources)
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                noiseReplacer.insertSource(refChildRecord.getId())
                for plugin in self.plugins.iter():
                    plugin.measure(measChildRecord, exposure, refChildRecord, refWcs)
                noiseReplacer.removeSource(refChildRecord.getId())

            # then process the parent record
            noiseReplacer.insertSource(refParentRecord.getId())
            for plugin in self.plugins.iter():
                plugin.measure(measParentRecord, exposure, refParentRecord, refWcs)
            for plugin in self.plugins.iterN():
                plugin.measureN(measChildCat, exposure, refChildCat)
                plugin.measureN(measParentCat[parentIdx:parentIdx+1], exposure,
                                refParentCat[parentIdx:parentIdx+1])
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()
        return Struct(sources=sources)

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
        raise NotImplementedError()

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
        catalog = lsst.afw.table.SourceCatalog(self.mapper.getOutputSchema())
        datasetType = self.dataPrefix + "forced"
        return {datasetType:catalog}

    def generateSources(self, dataRef, exposure, refCat, refWcs):
        """Generate sources to be measured, copying any fields in self.config.copyColumns
        Also, transform footprints to the measurement coordinate system, noting that
        we do not currently have the ability to transform heavy footprints because they
        need deblender information that we don't have.  So when the reference and measure
        WCS are different, we won't be able to use the heavyFootprint child info at all.

        @param dataRef     Data reference from butler
        @param exposure    Exposure to be measured
        @param refCat      Sequence (not necessarily a SourceCatalog) of reference sources
        @param idFactory   Factory to generate unique ids for forced sources
        @return Source catalog ready for measurement
        """
        idFactory = self.makeIdFactory(dataRef)
        table = lsst.afw.table.SourceTable.make(self.mapper.getOutputSchema(), idFactory)
        sources = lsst.afw.table.SourceCatalog(table)
        table = sources.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(refCat))
        expregion = exposure.getBBox(lsst.afw.image.PARENT)
        targetWcs = exposure.getWcs()
        for ref in refCat:
            newsource = sources.addNew()
            newsource.assign(ref, self.mapper)
            footprint = newsource.getFootprint()
            # if heavy, just transform the "light footprint" and leave the rest behind
            if footprint.isHeavy():
                footprint = lsst.afw.detection.Footprint(footprint)
            if not refWcs == targetWcs:
                footprint = footprint.transform(refWcs, targetWcs, expregion, True)
            newsource.setFootprint(footprint)
        return sources

    def _getConfigName(self):
        """Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return self.dataPrefix + "forced_config"

    def _getMetadataName(self):
        """Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return self.dataPrefix + "forced_metadata"

