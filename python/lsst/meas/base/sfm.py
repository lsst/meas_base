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
"""Base classes for single-frame measurement plugins and the driver task for these.

In single-frame measurement, we assume that detection and probably deblending have already been run on
the same frame, so a SourceCatalog has already been created with Footprints (which may be HeavyFootprints).
Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).
"""

import lsst.pipe.base as pipeBase

from .pluginRegistry import PluginRegistry
from .baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
                              BaseMeasurementConfig, BaseMeasurementTask)
from .noiseReplacer import NoiseReplacer, DummyNoiseReplacer

__all__ = ("SingleFramePluginConfig", "SingleFramePlugin",
           "SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")


class SingleFramePluginConfig(BaseMeasurementPluginConfig):
    """
    Base class for configs of single-frame plugin algorithms.
    """
    pass


class SingleFramePlugin(BaseMeasurementPlugin):
    """Base class for single-frame plugin algorithms.

    Parameters
    ----------
    config :
        An instance of this class's ConfigClass.
    name :
        The string the plugin was registered with.
    schema :
        The Source schema.  New fields should be added here to
        hold measurements produced by this plugin.
    metadata :
        Plugin metadata that will be attached to the output catalog
    logName :
        May include logName to use for this plugin

    Notes
    -----
    New Plugins can be created in Python by inheriting directly from this class
    and implementing measure(), fail() (from BasePlugin), and optionally __init__
    and measureN().  Plugins can also be defined in C++ via the WrappedSingleFramePlugin
    class.
    """

    # All subclasses of SingleFramePlugin should be registered here
    registry = PluginRegistry(SingleFramePluginConfig)
    ConfigClass = SingleFramePluginConfig

    def __init__(self, config, name, schema, metadata, logName=None, **kwds):
        BaseMeasurementPlugin.__init__(self, config, name, logName=logName)

    def measure(self, measRecord, exposure):
        """Measure the properties of a source on a single image (single-epoch image or coadd).

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            to be filled with outputs,
            and from which previously-measured quantities can be
            retreived.

        exposure : `lsst.afw.image.ExposureF`
            containing the pixel data to be measured and the associated Psf, Wcs, etc. All
            other sources in the image will have been replaced by noise according to deblender outputs.
        """
        raise NotImplementedError()

    def measureN(self, measCat, exposure):
        """Measure the properties of a group of blended sources on a single image
        (single-epoch image or coadd).

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            to be filled with outputs, and from which previously-measured quantities can be
            retrieved, containing only the sources that should be measured together in this call.

        exposure : `lsst.afw.image.ExposureF`
            containing the pixel data to be measured and the associated Psf, Wcs, etc.  Sources
            not in the blended hierarchy to be measured will have been
            replaced with noise using deblender outputs.

        Notes
        -----
        Derived classes that do not implement measureN() should just inherit this
        disabled version.  Derived classes that do implement measureN() should additionally
        add a bool doMeasureN config field to their config class to signal that measureN-mode
        is available.
        """
        raise NotImplementedError()


class SingleFrameMeasurementConfig(BaseMeasurementConfig):
    """
    Config class for single frame measurement driver task.
    """

    plugins = SingleFramePlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_SdssCentroid",
                 "base_NaiveCentroid",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_PsfFlux",
                 "base_CircularApertureFlux",
                 "base_SkyCoord",
                 "base_Variance",
                 "base_Blendedness",
                 "base_LocalBackground",
                 ],
        doc="Plugins to be run and their configuration"
    )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")
    undeblended = SingleFramePlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )

## @addtogroup LSST_task_documentation
## @{
## @page SingleFrameMeasurementTask
## @ref SingleFrameMeasurementTask_ "SingleFrameMeasurementTask"
## @copybrief SingleFrameMeasurementTask
## @}


class SingleFrameMeasurementTask(BaseMeasurementTask):
    """A subtask for measuring the properties of sources on a single exposure.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        to be initialized to include the
        measurement fields from the plugins already
    algMetadata : `lsst.daf.base.PropertyList`
        used to record information about
        each algorithm.  An empty PropertyList will be created if None.
    kwds:
        Keyword arguments forwarded to lsst.pipe.base.Task.__init__

    Notes
    -----
    The task is configured with a list of "plugins": each plugin defines the values it
    measures (i.e. the columns in a table it will fill) and conducts that measurement
    on each detected source (see SingleFramePlugin).  The job of the
    measurement task is to initialize the set of plugins (which includes setting up the
    catalog schema) from their configuration, and then invoke each plugin on each
    source.

    When run after the deblender (see lsst.meas.deblender.SourceDeblendTask),
    SingleFrameMeasurementTask also replaces each source's neighbors with noise before
    measuring each source, utilizing the HeavyFootprints created by the deblender (see
    NoiseReplacer).

    SingleFrameMeasurementTask has only two methods: __init__() and run().  For configuration
    options, see SingleFrameMeasurementConfig.

    A complete example of using SingleFrameMeasurementTask

    The code below is in examples/runSingleFrameTask.py

    See meas_algorithms_detection_Example for more information on SourceDetectionTask.

    First, import the required tasks (there are some other standard imports;
    read the file if you're confused):

    .. code-block:py
        from lsst.meas.algorithms.detection import SourceDetectionTask
        from lsst.meas.base import SingleFrameMeasurementTask

    We need to create our tasks before processing any data as the task constructors
    can add extra columns to the schema.  The most important argument we pass these to these
    is an lsst.afw.table.Schema object, which contains information about the fields (i.e. columns) of the
    measurement catalog we'll create, including names, types, and additional documentation.
    Tasks that operate on a catalog are typically passed a Schema upon construction, to which
    they add the fields they'll fill later when run.  We construct a mostly empty Schema that
    contains just the fields required for a SourceCatalog like this:

    .. code-block:py
        schema = afwTable.SourceTable.makeMinimalSchema()

    Now we can configure and create the SourceDetectionTask:

    .. code-block:py
        # Create the detection task
        #
        config = SourceDetectionTask.ConfigClass()
        config.thresholdPolarity = "both"
        config.background.isNanSafe = True
        config.thresholdValue = 3
        detectionTask = SourceDetectionTask(config=config, schema=schema)

    We then move on to configuring the measurement task:

    .. code-block:py
        config = SingleFrameMeasurementTask.ConfigClass()

    While a reasonable set of plugins is configured by default, we'll customize the list.
    We also need to unset one of the slots at the same time, because we're
    not running the algorithm that it's set to by default, and that would cause problems later:

    .. code-block:py
        config.plugins.names.clear()
        for plugin in ["base_SdssCentroid", "base_SdssShape",
        "base_CircularApertureFlux", "base_GaussianFlux"]:
        config.plugins.names.add(plugin)
        config.slots.psfFlux = None

    Now, finally, we can construct the measurement task:

    .. code-block:py
        measureTask = SingleFrameMeasurementTask(schema, config=config)

    After constructing all the tasks, we can inspect the Schema we've created:

    .. code-block:py

        print(schema)

    All of the fields in the
    schema can be accessed via the get() method on a record object.  See afwTable for more
    information.

    We're now ready to process the data (we could loop over multiple exposures/catalogs using the same
    task objects).  First create the output table and process the image to find sources:

    .. code-block:py

        tab = afwTable.SourceTable.make(schema)
        sources = result.sources

    Then measure them:

    .. code-block:py

        measureTask.run(sources, exposure)

    We then might plot the results (e.g. if you set `--ds9` on the command line)

    .. code-block:py

        if display:
        # display on ds9 (see also --debug argparse option)
            frame = 1
            ds9.mtv(exposure, frame=frame)

            with ds9.Buffering():
                for s in sources:
                    xy = s.getCentroid()
                    ds9.dot('+', *xy, ctype=ds9.CYAN if s.get("flags_negative") else ds9.GREEN, frame=frame)
                    ds9.dot(s.getShape(), *xy, ctype=ds9.RED, frame=frame)

    and end up with something like
    html runSingleFrameTask-ds9.png
    """

    ConfigClass = SingleFrameMeasurementConfig

    NOISE_SEED_MULTIPLIER = "NOISE_SEED_MULTIPLIER"
    NOISE_SOURCE = "NOISE_SOURCE"
    NOISE_OFFSET = "NOISE_OFFSET"
    NOISE_EXPOSURE_ID = "NOISE_EXPOSURE_ID"

    def __init__(self, schema, algMetadata=None, **kwds):
        super(SingleFrameMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.schema = schema
        self.config.slots.setupSchema(self.schema)
        self.initializePlugins(schema=self.schema)

        # Check to see if blendedness is one of the plugins
        if 'base_Blendedness' in self.plugins:
            self.doBlendedness = True
            self.blendPlugin = self.plugins['base_Blendedness']
        else:
            self.doBlendedness = False

    @pipeBase.timeMethod
    def run(self, measCat, exposure, noiseImage=None, exposureId=None, beginOrder=None, endOrder=None):
        """Run single frame measurement over an exposure and source catalog

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            to be filled with outputs. Must contain all the SourceRecords to be measured (with Footprints
            attached), and have a schema that is a superset of self.schema.
        exposure : `lsst.afw.image.ExposureF`
            containing the pixel data to
            be measured and the associated Psf, Wcs, etc.
        noiseImage : `lsst.afw.image.ImageF`
            for test which need to control
            noiseReplacement
        exposureId :
            optional unique exposureId used to calculate random number
            generator seed in the NoiseReplacer.
        beginOrder :
            beginning execution order (inclusive): measurements with
            executionOrder < beginOrder are not executed. None for no limit.
        endOrder :
            ending execution order (exclusive): measurements with
            executionOrder >= endOrder are not executed. None for no limit.
        """
        assert measCat.getSchema().contains(self.schema)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in measCat}

        # noiseReplacer is used to fill the footprints with noise and save heavy footprints
        # of the source pixels so that they can be restored one at a time for measurement.
        # After the NoiseReplacer is constructed, all pixels in the exposure.getMaskedImage()
        # which belong to objects in measCat will be replaced with noise
        if self.config.doReplaceWithNoise:
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure, footprints,
                                          noiseImage=noiseImage, log=self.log, exposureId=exposureId)
            algMetadata = measCat.getMetadata()
            if algMetadata is not None:
                algMetadata.addInt(self.NOISE_SEED_MULTIPLIER, self.config.noiseReplacer.noiseSeedMultiplier)
                algMetadata.addString(self.NOISE_SOURCE, self.config.noiseReplacer.noiseSource)
                algMetadata.addDouble(self.NOISE_OFFSET, self.config.noiseReplacer.noiseOffset)
                if exposureId is not None:
                    algMetadata.addLong(self.NOISE_EXPOSURE_ID, exposureId)
        else:
            noiseReplacer = DummyNoiseReplacer()

        self.runPlugins(noiseReplacer, measCat, exposure, beginOrder, endOrder)

    def runPlugins(self, noiseReplacer, measCat, exposure, beginOrder=None, endOrder=None):
        """Function which calls the defined measument plugins on an exposure

        Parameters
        ----------
        noiseReplacer : lsst.meas.base.NoiseReplacer
            noiseReplacer to fill sources not being measured with noise.

        measCat : lsst.afw.table.SourceCatalog
            SourceCatalog to be filled with outputs. Must contain all the SourceRecords to be measured (with
            Footprints attached), and have a schema that is a superset of self.schema.

        exposure : lsst.afw.image.ExposureF
            Exposure contaning the pixel data to be measured and the associated PSF, WCS, etc.

        beginOrder : float
            beginning execution order (inclusive): measurements with executionOrder < beginOrder are not
            executed. None for no limit.

        endOrder : float
            ending execution order (exclusive): measurements with executionOrder >= endOrder are not
            executed. None for no limit.
        """
        # First, create a catalog of all parentless sources
        # Loop through all the parent sources, first processing the children, then the parent
        measParentCat = measCat.getChildren(0)

        nMeasCat = len(measCat)
        nMeasParentCat = len(measParentCat)
        self.log.info("Measuring %d source%s (%d parent%s, %d child%s) ",
                      nMeasCat, ("" if nMeasCat == 1 else "s"),
                      nMeasParentCat, ("" if nMeasParentCat == 1 else "s"),
                      nMeasCat - nMeasParentCat, ("" if nMeasCat - nMeasParentCat == 1 else "ren"))

        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                self.callMeasure(measChildRecord, exposure, beginOrder=beginOrder, endOrder=endOrder)

                if self.doBlendedness:
                    self.blendPlugin.cpp.measureChildPixels(exposure.getMaskedImage(), measChildRecord)

                noiseReplacer.removeSource(measChildRecord.getId())

            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            self.callMeasure(measParentRecord, exposure, beginOrder=beginOrder, endOrder=endOrder)

            if self.doBlendedness:
                self.blendPlugin.cpp.measureChildPixels(exposure.getMaskedImage(), measParentRecord)

                # Finally, process both the parent and the child set through measureN
            self.callMeasureN(measParentCat[parentIdx:parentIdx+1], exposure,
                              beginOrder=beginOrder, endOrder=endOrder)
            self.callMeasureN(measChildCat, exposure, beginOrder=beginOrder, endOrder=endOrder)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()

        # Undeblended plugins only fire if we're running everything
        if endOrder is None:
            for source in measCat:
                for plugin in self.undeblendedPlugins.iter():
                    self.doMeasurement(plugin, source, exposure)

        # Now we loop over all of the sources one more time to compute the blendedness metrics
        if self.doBlendedness:
            for source in measCat:
                self.blendPlugin.cpp.measureParentPixels(exposure.getMaskedImage(), source)

    def measure(self, measCat, exposure):
        """
        Backwards-compatibility alias for run()
        """
        self.run(measCat, exposure)
