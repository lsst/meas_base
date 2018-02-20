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
HeavyFootprints (this is currently the case, at least for CCD forced photometry), then we will only
be able to replace objects with noise one deblend family at a time, and hence measurements run in
single-object mode may be contaminated by neighbors when run on objects with parent != 0.

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
from builtins import zip

import lsst.pex.config
import lsst.pipe.base

from .pluginRegistry import PluginRegistry
from .baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
                              BaseMeasurementConfig, BaseMeasurementTask)
from .noiseReplacer import NoiseReplacer, DummyNoiseReplacer

__all__ = ("ForcedPluginConfig", "ForcedPlugin",
           "ForcedMeasurementConfig", "ForcedMeasurementTask")


class ForcedPluginConfig(BaseMeasurementPluginConfig):
    """Base class for configs of forced measurement plugins."""
    pass


class ForcedPlugin(BaseMeasurementPlugin):

    # All subclasses of ForcedPlugin should be registered here
    registry = PluginRegistry(ForcedPluginConfig)

    ConfigClass = ForcedPluginConfig

    def __init__(self, config, name, schemaMapper, metadata, logName=None):
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
        BaseMeasurementPlugin.__init__(self, config, name, logName=logName)

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
        default=["base_PixelFlags",
                 "base_TransformedCentroid",
                 "base_SdssCentroid",
                 "base_TransformedShape",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_CircularApertureFlux",
                 "base_PsfFlux",
                 "base_LocalBackground",
                 ],
        doc="Plugins to be run and their configuration"
    )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")
    undeblended = ForcedPlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )

    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent": "parentObjectId", "deblend_nChild": "deblend_nChild",
                 "coord_ra": "coord_ra", "coord_dec": "coord_dec"}
    )

    checkUnitsParseStrict = lsst.pex.config.Field(
        doc="Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'",
        dtype=str,
        default="raise",
    )

    def setDefaults(self):
        self.slots.centroid = "base_TransformedCentroid"
        self.slots.shape = "base_TransformedShape"
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = None
        self.slots.instFlux = None
        self.slots.calibFlux = None

## \addtogroup LSST_task_documentation
## \{
## \page ForcedMeasurementTask
## \ref ForcedMeasurementTask_ "ForcedMeasurementTask"
## \copybrief ForcedMeasurementTask
## \}


class ForcedMeasurementTask(BaseMeasurementTask):
    """!
    \anchor ForcedMeasurementTask_

    \brief A subtask for measuring the properties of sources on a single
    exposure, using an existing "reference" catalog to constrain some aspects
    of the measurement.

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

    ForcedMeasurementTask has only three methods: __init__(), run(), and generateMeasCat().
    For configuration options, see SingleFrameMeasurementConfig.

    Notes
    -----

    *Forced* measurement means that the plugins are provided with a reference
    source containing centroid and/or shape measurements that they may use
    however they see fit. Some plugins can use these to set the location and
    size of apertures, but others may choose to ignore this information,
    essentially performing an unforced measurement starting at the position
    of the reference source (which may nevertheless be useful for certain
    investigations). Knowing how the plugin uses the reference information is
    essential to interpreting its resulting measurements. Typically, centroid
    and shape measurement plugins (e.g., ``SdssCentroid`` and ``SdssShape``)
    are performing unforced measurements.
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
                                   later passed to generateMeasCat() and/or run().
        @param[in,out] algMetadata lsst.daf.base.PropertyList used to record information about
                                   each algorithm.  An empty PropertyList will be created if None.
        @param[in]     **kwds      Keyword arguments passed from lsst.pipe.base.Task.__init__
        """
        super(ForcedMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(), False)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        self.initializePlugins(schemaMapper=self.mapper)
        self.schema = self.mapper.getOutputSchema()
        self.schema.checkUnits(parse_strict=self.config.checkUnitsParseStrict)

    def run(self, measCat, exposure, refCat, refWcs, exposureId=None, beginOrder=None, endOrder=None):
        """!
        Perform forced measurement.

        @param[in]  exposure     lsst.afw.image.ExposureF to be measured; must have at least a Wcs attached.
        @param[in]  measCat      Source catalog for measurement results; must be initialized with empty
                                 records already corresponding to those in refCat (via e.g. generateMeasCat).
        @param[in]  refCat       A sequence of SourceRecord objects that provide reference information
                                 for the measurement.  These will be passed to each Plugin in addition
                                 to the output SourceRecord.
        @param[in]  refWcs       Wcs that defines the X,Y coordinate system of refCat
        @param[in]  exposureId   optional unique exposureId used to calculate random number
                                 generator seed in the NoiseReplacer.
        @param[in]  beginOrder   beginning execution order (inclusive): measurements with
                                 executionOrder < beginOrder are not executed. None for no limit.
        @param[in]  endOrder     ending execution order (exclusive): measurements with
                                 executionOrder >= endOrder are not executed. None for no limit.

        Fills the initial empty SourceCatalog with forced measurement results.  Two steps must occur
        before run() can be called:
         - generateMeasCat() must be called to create the output measCat argument.
         - Footprints appropriate for the forced sources must be attached to the measCat records.  The
           attachTransformedFootprints() method can be used to do this, but this degrades HeavyFootprints
           to regular Footprints, leading to non-deblended measurement, so most callers should provide
           Footprints some other way.  Typically, calling code will have access to information that will
           allow them to provide HeavyFootprints - for instance, ForcedPhotCoaddTask uses the HeavyFootprints
           from deblending run in the same band just before non-forced is run measurement in that band.
        """
        # First check that the reference catalog does not contain any children for which
        # any member of their parent chain is not within the list.  This can occur at
        # boundaries when the parent is outside and one of the children is within.
        # Currently, the parent chain is always only one deep, but just in case, this
        # code checks for any case where when the parent chain to a child's topmost
        # parent is broken and raises an exception if it occurs.
        #
        # I.e. this code checks that this precondition is satisfied by whatever reference
        # catalog provider is being paired with it.
        refCatIdDict = {ref.getId(): ref.getParent() for ref in refCat}
        for ref in refCat:
            refId = ref.getId()
            topId = refId
            while(topId > 0):
                if topId not in refCatIdDict:
                    raise RuntimeError("Reference catalog contains a child for which at least "
                                       "one parent in its parent chain is not in the catalog.")
                topId = refCatIdDict[topId]

        # Construct a footprints dict which looks like
        # {ref.getId(): (ref.getParent(), source.getFootprint())}
        # (i.e. getting the footprint from the transformed source footprint)
        footprints = {ref.getId(): (ref.getParent(), measRecord.getFootprint())
                      for (ref, measRecord) in zip(refCat, measCat)}

        self.log.info("Performing forced measurement on %d source%s", len(refCat),
                      "" if len(refCat) == 1 else "s")

        if self.config.doReplaceWithNoise:
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure,
                                          footprints, log=self.log, exposureId=exposureId)
            algMetadata = measCat.getTable().getMetadata()
            if algMetadata is not None:
                algMetadata.addInt("NOISE_SEED_MULTIPLIER", self.config.noiseReplacer.noiseSeedMultiplier)
                algMetadata.addString("NOISE_SOURCE", self.config.noiseReplacer.noiseSource)
                algMetadata.addDouble("NOISE_OFFSET", self.config.noiseReplacer.noiseOffset)
                if exposureId is not None:
                    algMetadata.addLong("NOISE_EXPOSURE_ID", exposureId)
        else:
            noiseReplacer = DummyNoiseReplacer()

        # Create parent cat which slices both the refCat and measCat (sources)
        # first, get the reference and source records which have no parent
        refParentCat, measParentCat = refCat.getChildren(0, measCat)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat, measParentCat)):

            # first process the records which have the current parent as children
            refChildCat, measChildCat = refCat.getChildren(refParentRecord.getId(), measCat)
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                noiseReplacer.insertSource(refChildRecord.getId())
                self.callMeasure(measChildRecord, exposure, refChildRecord, refWcs,
                                 beginOrder=beginOrder, endOrder=endOrder)
                noiseReplacer.removeSource(refChildRecord.getId())

            # then process the parent record
            noiseReplacer.insertSource(refParentRecord.getId())
            self.callMeasure(measParentRecord, exposure, refParentRecord, refWcs,
                             beginOrder=beginOrder, endOrder=endOrder)
            self.callMeasureN(measParentCat[parentIdx:parentIdx+1], exposure,
                              refParentCat[parentIdx:parentIdx+1],
                              beginOrder=beginOrder, endOrder=endOrder)
            # measure all the children simultaneously
            self.callMeasureN(measChildCat, exposure, refChildCat,
                              beginOrder=beginOrder, endOrder=endOrder)
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()

        # Undeblended plugins only fire if we're running everything
        if endOrder is None:
            for measRecord, refRecord in zip(measCat, refCat):
                for plugin in self.undeblendedPlugins.iter():
                    self.doMeasurement(plugin, measRecord, exposure, refRecord, refWcs)

    def generateMeasCat(self, exposure, refCat, refWcs, idFactory=None):
        """!Initialize an output SourceCatalog using information from the reference catalog.

        This generates a new blank SourceRecord for each record in refCat.  Note that this
        method does not attach any Footprints.  Doing so is up to the caller (who may
        call attachedTransformedFootprints or define their own method - see run() for more
        information).

        @param[in] exposure    Exposure to be measured
        @param[in] refCat      Sequence (not necessarily a SourceCatalog) of reference SourceRecords.
        @param[in] refWcs      Wcs that defines the X,Y coordinate system of refCat
        @param[in] idFactory   factory for creating IDs for sources

        @return    Source catalog ready for measurement
        """
        if idFactory is None:
            idFactory = lsst.afw.table.IdFactory.makeSimple()
        table = lsst.afw.table.SourceTable.make(self.schema, idFactory)
        measCat = lsst.afw.table.SourceCatalog(table)
        table = measCat.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(refCat))
        for ref in refCat:
            newSource = measCat.addNew()
            newSource.assign(ref, self.mapper)
        return measCat

    def attachTransformedFootprints(self, sources, refCat, exposure, refWcs):
        """!Default implementation for attaching Footprints to blank sources prior to measurement

        Footprints for forced photometry must be in the pixel coordinate system of the image being
        measured, while the actual detections may start out in a different coordinate system.
        This default implementation transforms the Footprints from the reference catalog from the
        refWcs to the exposure's Wcs, which downgrades HeavyFootprints into regular Footprints,
        destroying deblend information.

        Note that ForcedPhotImageTask delegates to this method in its own attachFootprints method.
        attachFootprints can then be overridden by its subclasses to define how their Footprints
        should be generated.

        See the documentation for run() for information about the relationships between run(),
        generateMeasCat(), and attachTransformedFootprints().
        """
        exposureWcs = exposure.getWcs()
        region = exposure.getBBox(lsst.afw.image.PARENT)
        for srcRecord, refRecord in zip(sources, refCat):
            srcRecord.setFootprint(refRecord.getFootprint().transform(refWcs, exposureWcs, region))
