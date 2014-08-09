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

Command-line driver tasks for forced measurement can be found in forcedImage.py, including ForcedImageTask,
ForcedCcdTask, and ForcedCoaddTask.
"""

import lsst.pex.config
import lsst.pipe.base

from .base import *

__all__ = ("ForcedPluginConfig", "ForcedPlugin", "WrappedForcedPlugin",
           "ForcedMeasurementConfig", "ForcedMeasurementTask")

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

class WrappedForcedPlugin(ForcedPlugin):
    """A base class for ForcedPlugins that delegates the algorithmic work to a C++
    Algorithm class.

    Derived classes of WrappedForcedPlugin must set the AlgClass class attribute
    to the C++ class being wrapped, which must meet the requirements defined in the
    "Implementing New Plugins and Algorithms" section of the meas_base documentation.  This is usually done
    by calling the generate() class method.
    """

    AlgClass = None

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.resultMapper = self.AlgClass.makeResultMapper(schema, name, config.makeControl())
        # TODO: check flags

    def measure(self, measRecord, exposure, refRecord, refWcs):
        inputs = self.AlgClass.Input(measRecord)
        result = self.AlgClass.apply(exposure, inputs, self.config.makeControl())
        self.resultMapper.apply(measRecord, result)

    def measureN(self, measCat, exposure, refCat, refWcs):
        assert hasattr(AlgClass, "applyN")  # would be better if we could delete this method somehow
        inputs = self.AlgClass.Input.Vector(measCat)
        results = self.AlgClass.applyN(exposure, inputs, self.config.makeControl())
        for result, measRecord in zip(results, measCat):
            self.resultMapper.apply(measRecord, result)

    def fail(self, measRecord, error=None):
        # The ResultMapper will set detailed flag bits describing the error if error is not None,
        # and set a general failure bit otherwise.
        self.resultMapper.fail(measRecord, error)

    @classmethod
    def generate(Base, AlgClass, name=None, doRegister=True, ConfigClass=None, executionOrder=None):
        """Create a new derived class of WrappedForcedPlugin from a C++ Algorithm class.

        @param[in]   AlgClass   The name of the (Swigged) C++ Algorithm class this Plugin will delegate to.
        @param[in]   name       The name to use when registering the Plugin (ignored if doRegister=False).
                                Defaults to the result of generateAlgorithmName(AlgClass).
        @param[in]   doRegister   If True (default), register the new Plugin so it can be configured to be
                                  run by ProcessImageForcedTask.
        @param[in]   ConfigClass  The ConfigClass associated with the new Plugin.  This should have a
                                  makeControl() method that returns the Control object used by the C++
                                  Algorithm class.
        @param[in]   executionOrder   If not None, a float that overrides the default executionOrder
                                      for this algorithm (see BasePluginConfig.executionOrder).

        For more information, please see the "Adding New Algorithms" section of the main meas_base
        documentation.
        """
        if ConfigClass is None:
            ConfigClass = lsst.pex.config.makeConfigClass(AlgClass.Control, base=Base.ConfigClass,
                                                          module=AlgClass.__module__)
            if hasattr(AlgClass, "applyN"):
                ConfigClass.doMeasureN = lsst.pex.config.Field(
                    dtype=bool, default=True,
                    doc="whether to run this plugin multi-object mode"
                    )
        if executionOrder is not None:
            ConfigClass.executionOrder.default = float(executionOrder)
        PluginClass = type(AlgClass.__name__ + "ForcedPlugin", (Base,),
                           dict(AlgClass=AlgClass, ConfigClass=ConfigClass))
        if doRegister:
            if name is None:
                name = generateAlgorithmName(AlgClass)
            Base.registry.register(name, PluginClass)
        return PluginClass

class ForcedMeasurementConfig(BaseMeasurementConfig):
    """Config class for forced measurement driver task."""

    plugins = ForcedPlugin.registry.makeField(
        multi=True,
        default=["centroid.peak",
                 ],
        doc="Plugins to be run and their configuration"
        )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")

    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent":"parentObjectId"}
        )

    def setDefaults(self):
        self.slots.centroid = "centroid.peak"
        self.slots.shape = None
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = None
        self.slots.instFlux = None

class ForcedMeasurementTask(BaseMeasurementTask):
    """Forced measurement driver task

    This task is intended to be a subtask of another command line task, such as ProcessImageForcedTask
    It should be subclassed for running on CCDs and Coadds.  See ProcessCcd(Coadd)ForcedTask.

    Note that the __init__ method takes the refSchema, while the run method takes the refCat.
    The two schemas must match.
    """

    ConfigClass = ForcedMeasurementConfig

    def __init__(self, refSchema, algMetadata=None, flags=None, **kwds):
        BaseMeasurementTask.__init__(self, algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema())
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, self.mapper, flags=flags,
                                             others=self.plugins, metadata=self.algMetadata)

    def run(self, exposure, refCat, refWcs, idFactory=None):
        """The reference list and refWcs are passed in to the run method, along with the exposure.
        It creates its own source catalog, which is returned in the result structure.
        The IdFactory is assumed to be known by the parentTask (see ProcessImageForcedTask.run)
        and will default to IdFactory.makeSimple() if not specified.
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
        self.config.slots.setupTable(sources.table)

        # convert the footprints to the coordinate system of the exposure
        noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure, footprints, log=self.log)

        # Create parent cat which slices both the referenceCat and measCat (sources)
        # first, get the reference and source records which have no parent
        refParentCat, measParentCat = referenceCat.getChildren(0, sources)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat,measParentCat)):

            # first process the records which have the current parent as children
            refChildCat, measChildCat = referenceCat.getChildren(refParentRecord.getId(), sources)
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                noiseReplacer.insertSource(refChildRecord.getId())
                callMeasure(self, measChildRecord, exposure, refChildRecord, refWcs)
                noiseReplacer.removeSource(refChildRecord.getId())

            # then process the parent record
            noiseReplacer.insertSource(refParentRecord.getId())
            callMeasure(self, measParentRecord, exposure, refParentRecord, refWcs)
            callMeasureN(self, measChildCat, exposure, refChildCat)
            callMeasureN(self, measParentCat[parentIdx:parentIdx+1], exposure,
                         refParentCat[parentIdx:parentIdx+1])
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()
        return lsst.pipe.base.Struct(sources=sources)

    def generateSources(self, exposure, refCat, refWcs, idFactory=None):
        """Generate sources to be measured, copying any fields in self.config.copyColumns
        Also, transform footprints to the measurement coordinate system, noting that
        we do not currently have the ability to transform heavy footprints because they
        need deblender information that we don't have.  So when the reference and measure
        WCS are different, we won't be able to use the heavyFootprint child info at all.

        @param exposure    Exposure to be measured
        @param refCat      Sequence (not necessarily a SourceCatalog) of reference sources
        @param refWcs      Wcs that goes with the X,Y cooridinates of the refCat
        @param idFactory   factory for making new ids from reference catalog ids
        @param algMetadata algorithm metadata to be attached to the catalog (PropertySet)
        @return Source catalog ready for measurement
        """
        if idFactory == None:
            idFactory = lsst.afw.table.IdFactory.makeSimple()
        table = lsst.afw.table.SourceTable.make(self.mapper.getOutputSchema(), idFactory)
        table.setVersion(self.TableVersion)
        sources = lsst.afw.table.SourceCatalog(table)
        table = sources.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(refCat))
        expRegion = exposure.getBBox(lsst.afw.image.PARENT)
        targetWcs = exposure.getWcs()
        for ref in refCat:
            newSource = sources.addNew()
            newSource.assign(ref, self.mapper)
            footprint = newSource.getFootprint()
            # if heavy, just transform the "light footprint" and leave the rest behind
            if footprint.isHeavy():
                footprint = lsst.afw.detection.Footprint(footprint)
            if not refWcs == targetWcs:
                footprint = footprint.transform(refWcs, targetWcs, expRegion, True)
            newSource.setFootprint(footprint)
        return sources
