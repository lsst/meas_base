#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST Corporation.
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

Command-line driver tasks for forced measurement can be found in forcedPhotImage.py, including
ForcedPhotImageTask, ForcedPhotCcdTask, and ForcedPhotCoaddTask.
"""

import lsst.pex.config
import lsst.pipe.base

from .base import *

__all__ = ("ForcedPluginConfig", "ForcedPlugin",
           "ForcedMeasurementConfig", "ForcedMeasurementTask")

class ForcedPluginConfig(BasePluginConfig):
    """Base class for configs of forced measurement plugins."""
    pass

class ForcedPlugin(BasePlugin):

    # All subclasses of ForcedPlugin should be registered here
    registry = PluginRegistry(ForcedPluginConfig)

    ConfigClass = ForcedPluginConfig

    def __init__(self, config, name, schemaMapper, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        @param[in,out]  schemaMapper  A SchemaMapper that maps reference catalog fields to output
                                      catalog fields.  Output fields should be added to the
                                      output schema.  While most plugins will not need to map
                                      fields from the reference schema, if they do so, those fields
                                      will be transferred before any plugins are run.
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        BasePlugin.__init__(self, config, name)

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
        default=["base_TransformedCentroid",
                 "base_TransformedShape",
                 "base_GaussianFlux",
                 "base_CircularApertureFlux",
                 "base_PsfFlux",
                 ],
        doc="Plugins to be run and their configuration"
        )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")

    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent":"parentObjectId"}
        )

    def setDefaults(self):
        self.slots.centroid = "base_TransformedCentroid"
        self.slots.shape = "base_TransformedShape"
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = None
        self.slots.instFlux = None

## @addtogroup LSST_task_documentation
## @{
## @page pageForcedMeasurementTask ForcedMeasurementTask
## ForcedMeasurementTask
## @copybrief ForcedMeasurementTask
## @}

class ForcedMeasurementTask(BaseMeasurementTask):
    """!
	A subtask for measuring the properties of sources on a single exposure, using an existing
    "reference" catalog to constrain some aspects of the measurement.

    The task is configured with a list of "plugins": each plugin defines the values it
    measures (i.e. the columns in a table it will fill) and conducts that measurement
    on each detected source (see ForcedPlugin).  The job of the
    measurement task is to initialize the set of plugins (which includes setting up the
    catalog schema) from their configuration, and then invoke each plugin on each
    source.

    Most of the time, ForcedMeasurementTask will be used via one of the subclasses of
    ForcedPhotImageTask, ForcedPhotCcdTask and ForcedPhotCoaddTask.  These combine
    this measurement subtask with a "references" subtask (see BaseReferencesTask and
    CoaddSrcReferencesTask) to perform forced measurement using measurements performed on
    another image as the references.  There is generally little reason to use
    ForcedMeasurementTask outside of one of these drivers, unless it is necessary to avoid
    using the Butler for I/O.

    ForcedMeasurementTask has only three methods: __init__(), run(), and generateSources().
    For configuration options, see SingleFrameMeasurementConfig.
    """

    ConfigClass = ForcedMeasurementConfig

    def __init__(self, refSchema, algMetadata=None, **kwds):
        """!
        Initialize the task.  Set up the execution order of the plugins and initialize
        the plugins, giving each plugin an opportunity to add its measurement fields to
        the output schema and to record information in the task metadata.

        Note that while SingleFrameMeasurementTask is passed an initial Schema that is
        appended to in order to create the output Schema, ForcedMeasurementTask is
        initialized with the Schema of the reference catalog, from which a new Schema
        for the output catalog is created.  Fields to be copied directly from the
        reference Schema are added before Plugin fields are added.

        @param[in]  refSchema      Schema of the reference catalog.  Must match the catalog
                                   later passed to generateSources() and/or run().
        @param[in,out] algMetadata lsst.daf.base.PropertyList used to record information about
                                   each algorithm.  An empty PropertyList will be created if None.
        @param[in]     **kwds      Keyword arguments passed from lsst.pipe.base.Task.__init__
        """
        super(ForcedMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema, lsst.afw.table.Schema(1))
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)
        self.initializePlugins(schemaMapper=self.mapper)
        self.schema = self.mapper.getOutputSchema()
        self.config.slots.setupSchema(self.schema)

    def run(self, exposure, refCat, refWcs, idFactory=None):
        """!
        Perform forced measurement.

        @param[in]  exposure     lsst.afw.image.ExposureF to be measured; must have at least a Wcs attached.
        @param[in]  refCat       A sequence of SourceRecord objects that provide reference information
                                 for the measurement.  These will be passed to each Plugin in addition
                                 to the output SourceRecord.
        @param[in]  refWcs       Wcs that defines the X,Y coordinate system of refCat
        @param[in]  idFactory    factory for creating IDs for sources

        @return SourceCatalog containing measurement results.

        Delegates creating the initial empty SourceCatalog to generateSources(), then fills it.
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
        sources = self.generateSources(exposure, refList, refWcs, idFactory)

        # Steal the transformed source footprint and use it to complete the footprints dict,
        # which then looks like {ref.getId(): (ref.getParent(), source.getFootprint())}
        for (reference, source) in zip(refList, sources):
            footprints[reference.getId()] = (reference.getParent(), source.getFootprint())

        self.log.info("Performing forced measurement on %d sources" % len(refList))

        # Build a catalog of just the references we intend to measure
        referenceCat = lsst.afw.table.SourceCatalog(self.mapper.getInputSchema())
        referenceCat.extend(refList)

        # convert the footprints to the coordinate system of the exposure
        if self.config.doReplaceWithNoise:
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure, footprints, log=self.log)
        else:
            noiseReplacer = DummyNoiseReplacer()

        # Create parent cat which slices both the referenceCat and measCat (sources)
        # first, get the reference and source records which have no parent
        refParentCat, measParentCat = referenceCat.getChildren(0, sources)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat,measParentCat)):

            # first process the records which have the current parent as children
            refChildCat, measChildCat = referenceCat.getChildren(refParentRecord.getId(), sources)
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                noiseReplacer.insertSource(refChildRecord.getId())
                self.callMeasure(measChildRecord, exposure, refChildRecord, refWcs)
                noiseReplacer.removeSource(refChildRecord.getId())

            # then process the parent record
            noiseReplacer.insertSource(refParentRecord.getId())
            self.callMeasure(measParentRecord, exposure, refParentRecord, refWcs)
            self.callMeasureN(measParentCat[parentIdx:parentIdx+1], exposure,
                              refParentCat[parentIdx:parentIdx+1])
            # measure all the children simultaneously
            self.callMeasureN(measChildCat, exposure, refChildCat)
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()
        return lsst.pipe.base.Struct(sources=sources)

    def generateSources(self, exposure, refCat, refWcs, idFactory=None):
        """!
        Initialize an output SourceCatalog using information from the reference catalog.

        This generate a new blank SourceRecord for each record in refCat, copying any
        fields in ForcedMeasurementConfig.copyColumns.  This also transforms the Footprints
        in refCat to the measurement coordinate system if it differs from refWcs, and attaches
        these to the new SourceRecords.  Note that we do not currently have the ability to
        transform heavy footprints, so when the reference and measure WCSs are different,
        HeavyFootprints will be converted to regular Footprints, which makes it impossible
        to properly measure blended objects.

        @param[in] exposure    Exposure to be measured
        @param[in] refCat      Sequence (not necessarily a SourceCatalog) of reference SourceRecords.
        @param[in] refWcs      Wcs that defines the X,Y coordinate system of refCat
        @param[in] idFactory   factory for creating IDs for sources

        @return Source catalog ready for measurement
        """
        if idFactory == None:
            idFactory = lsst.afw.table.IdFactory.makeSimple()
        table = lsst.afw.table.SourceTable.make(self.schema, idFactory)
        sources = lsst.afw.table.SourceCatalog(table)
        table = sources.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(refCat))
        expRegion = exposure.getBBox()
        targetWcs = exposure.getWcs()
        for ref in refCat:
            newSource = sources.addNew()
            newSource.assign(ref, self.mapper)
            footprint = newSource.getFootprint()
            if footprint is not None:
                # if heavy, just transform the "light footprint" and leave the rest behind
                if footprint.isHeavy():
                    footprint = lsst.afw.detection.Footprint(footprint)
                if not refWcs == targetWcs:
                    footprint = footprint.transform(refWcs, targetWcs, expRegion, True)
                newSource.setFootprint(footprint)
        return sources
