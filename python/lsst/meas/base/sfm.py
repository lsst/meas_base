"""Base classes for single-frame measurement plugins and the driver task for these.

In single-frame measurement, we assumes that detection and probably deblending have already been run on
the same frame, so a SourceCatalog has already been created with Footprints (which may be HeavyFootprints).
Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).
"""

import lsst.pex.config
import lsst.pipe.base
import lsst.daf.base
import lsst.meas.algorithms
from .base import *
__all__ = ("SingleFramePluginConfig", "SingleFramePlugin",
"SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")

class SingleFramePluginConfig(BasePluginConfig):
    """Base class for configs of single-frame plugins."""
    pass

class SingleFramePlugin(BasePlugin):
    """Base class for single-frame plugins."""


    # All subclasses of SingleFramePlugin should be registered here
    registry = PluginRegistry(SingleFramePluginConfig)
    ConfigClass = SingleFramePluginConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        @param[in,out]  schema   The Source schema.  New fields should be added here to
                                 hold measurements produced by this plugin.
        @param[in]  flags        A set of bitflags describing the data that the plugin
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       An PluginMap of previously-initialized plugins
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        self.config = config

    def measureSingle(self, exposure, source):
        """Measure the properties of a source on a single image
        (single-epoch image or coadd).

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  All
                                 other sources in the image will have been replaced by
                                 noise according to deblender outputs.

        @param[in,out] measRecord  lsst.afw.table.SourceRecord to be filled with outputs,
                                   and from which previously-measured quantities can be
                                   retreived.

        """
        raise NotImplementedError()

    def measureMulti(self, exposure, sources):
        """Measure the properties of a group of blended sources on a single image
        (single-epoch image or coadd).

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  Sources
                                 not in the blended hierarchy to be measured will have
                                 been replaced with noise using deblender outputs.

        @param[in,out] measCat   lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be
                                 retrieved, containing only the sources that should be
                                 measured together in this call.

        """
        raise NotImplementedError()

class SingleFrameMeasurementConfig(BaseMeasurementConfig):

    plugins = SingleFramePlugin.registry.makeField(
        multi=True,
        default=["centroid.peak",
                 "flags.pixel",
                 "centroid.gaussian",
                 "centroid.naive",
                 "shape.sdss",
                 "flux.gaussian",
                 "flux.naive",
                 "flux.psf",
                 "flux.sinc",
                 "classification.extendedness",
                 "skycoord",
                 ],
        doc="Plugin plugins to be run and their configuration"
        )

class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """Single-frame measurement driver task"""

    ConfigClass = SingleFrameMeasurementConfig
    _DefaultName = "measurement"

    """ Initialize the task, including setting up the execution order of the plugins
        and providing the task with the metadata and schema objects

        @param[in] schema      lsst.afw.table.Schema, which should have been initialized
                               to include the measurement fields from the plugins already
    """
    #   The algMetadata parameter is currently required by the pipe_base running mechanism
    #   This is a termporary state untile the plugins are converted.
    def __init__(self, schema, algMetadata=None, flags=None, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.schema = schema
        self.algMetadata = lsst.daf.base.PropertyList()
        self.plugins = PluginMap()
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, schema=schema, flags=flags,
                others=self.plugins, metadata=self.algMetadata)

    """ Run single frame measurement over an exposure and source catalog

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.

        @param[in, out] measCat  lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be retreived.

    """
    def run(self, exposure, measCat):
        assert measCat.getSchema().contains(self.schema)
        self.config.slots.setupTable(measCat.table, prefix=self.config.prefix)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
            for measRecord in measCat}

        # noiseReplacer is used to fill the footprints with noise and save off heavy footprints
        # of what was in the exposure beforehand.
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)

        # loop through all the parent sources, processing the children, then the parent
        measParentCat = measCat.getChildren(0)
        self.log.info("There are %d parent sources"%len(measParentCat))
        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                for plugin in self.plugins.iterSingle():
                    plugin.measureSingle(exposure, measChildRecord)
                noiseReplacer.removeSource(measChildRecord.getId())
            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            for plugin in self.plugins.iterSingle():
                plugin.measureSingle(exposure, measParentRecord)

            for plugin in self.plugins.iterMulti():
                plugin.measureMulti(exposure, measParentCat[parentIndex:parentIndex+1])
                plugin.measureMulti(exposure, measChildCat)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()

