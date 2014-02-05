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

__all__ = ("SingleFramePluginConfig", "SingleFramePlugin", "SingleFrameMeasurementConfig",
           "SingleFrameMeasurementTask")

class SingleFramePluginConfig(BasePluginConfig):
    """Base class for configs of single-frame plugin algorithms."""
    pass

class SingleFramePlugin(BasePlugin):
    """Base class for single-frame plugin algorithms."""

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
        @param[in]  others       A PluginMap of previously-initialized plugins
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        self.config = config
        self.name = name

    def measure(self, exposure, source):
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

    def measureN(self, exposure, sources):
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
                 ],
        doc="Plugin plugins to be run and their configuration"
        )

class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """Single-frame measurement driver task"""

    ConfigClass = SingleFrameMeasurementConfig
    _DefaultName = "measurement"

    #   The algMetadata parameter is currently required by the pipe_tasks running mechanism
    #   This is a temporary state until pipe_tasks is converted to the new plugin framework.
    def __init__(self, schema, algMetadata=None, flags=None, **kwds):
        """Initialize the task, including setting up the execution order of the plugins
        and providing the task with the metadata and schema objects

        @param[in] schema      lsst.afw.table.Schema, which should have been initialized
                               to include the measurement fields from the plugins already
        """
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.schema = schema
        self.algMetadata = lsst.daf.base.PropertyList()
        self.plugins = PluginMap()
        # Init the plugins, sorted by execution order.  At the same time add to the schema
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, schema=schema, flags=flags,
                others=self.plugins, metadata=self.algMetadata)

    def run(self, exposure, measCat):
        """ Run single frame measurement over an exposure and source catalog

            @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                     be measured and the associated Psf, Wcs, etc.

            @param[in, out] measCat  lsst.afw.table.SourceCatalog to be filled with outputs,
                                     and from which previously-measured quantities can be retreived.

        """
        assert measCat.getSchema().contains(self.schema)
        self.config.slots.setupTable(measCat.table, prefix=self.config.prefix)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
            for measRecord in measCat}

        # noiseReplacer is used to fill the footprints with noise and save off heavy footprints
        # of the source pixels so that they can re restored one at a time for measurement.
        # After the NoiseReplacer is constructed, all pixels in the exposure.getMaskedImage()
        # which belong to objects in measCat will be replaced with noise
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)

        # First, create a catalog of all parentless sources
        # Loop through all the parent sources, first processing the children, then the parent
        measParentCat = measCat.getChildren(0)
        self.log.info("There are %d parent sources"%len(measParentCat))
        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                for plugin in self.plugins.iter():
                    plugin.measure(exposure, measChildRecord)
                noiseReplacer.removeSource(measChildRecord.getId())
            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            for plugin in self.plugins.iter():
                plugin.measure(exposure, measParentRecord)
            # Finally, process both the parent and the child set through measureN
            for plugin in self.plugins.iterN():
                plugin.measureN(exposure, measParentCat[parentIndex:parentIndex+1])
                plugin.measureN(exposure, measChildCat)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()

