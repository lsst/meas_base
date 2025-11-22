# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

r"""Base classes for single-frame measurement plugins and the associated task.

In single-frame measurement, we assume that detection and probably deblending
have already been run on the same frame, so a `~lsst.afw.table.SourceCatalog`
has already been created with `lsst.afw.detection.Footprint`\ s (which may be
"heavy" — that is, include pixel data). Measurements are generally recorded in
the coordinate system of the image being measured (and all slot-eligible
fields must be), but non-slot fields may be recorded in other coordinate
systems if necessary to avoid information loss (this should, of course, be
indicated in the field documentation).
"""

from lsst.utils.logging import PeriodicLogger
from lsst.utils.timer import timeMethod

from .pluginRegistry import PluginRegistry
from .baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
                              BaseMeasurementConfig, BaseMeasurementTask)

__all__ = ("SingleFramePluginConfig", "SingleFramePlugin",
           "SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")


class SingleFramePluginConfig(BaseMeasurementPluginConfig):
    """Base class for single-frame plugin configuration classes.
    """
    pass


class SingleFramePlugin(BaseMeasurementPlugin):
    """Base class for single-frame measurement plugin.

    Parameters
    ----------
    config : `SingleFramePlugin.ConfigClass`
        Configuration for this plugin.
    name : `str`
        The string with which the plugin was registered.
    schema : `lsst.afw.table.Schema`
        The schema for the source table . New fields are added here to
        hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    logName : `str`, optional
        Name to use when logging errors.

    Notes
    -----
    New plugins can be created in Python by inheriting directly from this
    class and implementing the `measure`, `fail` (from `BasePlugin`), and
    optionally `__init__` and `measureN` methods.  Plugins can also be defined
    in C++ via the `WrappedSingleFramePlugin` class.
    """

    registry = PluginRegistry(SingleFramePluginConfig)
    """Registry of subclasses of `SingleFramePlugin` (`PluginRegistry`).
    """

    ConfigClass = SingleFramePluginConfig

    def __init__(self, config, name, schema, metadata, logName=None, **kwds):
        BaseMeasurementPlugin.__init__(self, config, name, logName=logName)

    def measure(self, measRecord, exposure):
        """Measure the properties of a source on a single image.

        The image may be from a single epoch, or it may be a coadd.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record describing the object being measured. Previously-measured
            quantities may be retrieved from here, and it will be updated
            in-place tih the outputs of this plugin.
        exposure : `lsst.afw.image.ExposureF`
            The pixel data to be measured, together with the associated PSF,
            WCS, etc. All other sources in the image should have been replaced
            by noise according to deblender outputs.
        """
        raise NotImplementedError()


class SingleFrameMeasurementConfig(BaseMeasurementConfig):
    """Config class for single frame measurement driver task.
    """

    plugins = SingleFramePlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_SdssCentroid",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_PsfFlux",
                 "base_CircularApertureFlux",
                 "base_SkyCoord",
                 "base_Variance",
                 "base_Blendedness",
                 "base_LocalBackground",
                 "base_CompensatedTophatFlux",
                 "base_ClassificationSizeExtendedness",
                 ],
        doc="Plugins to be run and their configuration"
    )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")
    undeblended = SingleFramePlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )


class SingleFrameMeasurementTask(BaseMeasurementTask):
    """A subtask for measuring the properties of sources on a single exposure.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Schema of the output resultant catalog. Will be updated to provide
        fields to accept the outputs of plugins which will be executed by this
        task.
    algMetadata : `lsst.daf.base.PropertyList`, optional
        Used to record metadaa about algorithm execution. An empty
        `lsst.daf.base.PropertyList` will be created if `None`.
    **kwds
        Keyword arguments forwarded to `BaseMeasurementTask`.
    """

    ConfigClass = SingleFrameMeasurementConfig

    def __init__(self, schema, algMetadata=None, **kwds):
        super(SingleFrameMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.schema = schema
        self.config.slots.setupSchema(self.schema)
        self.initializePlugins(schema=self.schema)
        self.addInvalidPsfFlag(self.schema)

        # Check to see if blendedness is one of the plugins
        if 'base_Blendedness' in self.plugins:
            self.doBlendedness = True
            self.blendPlugin = self.plugins['base_Blendedness']
        else:
            self.doBlendedness = False

    @timeMethod
    def run(
        self,
        measCat,
        exposure,
        noiseImage=None,
        exposureId=None,
        beginOrder=None,
        endOrder=None,
        footprints=None,
    ):
        r"""Run single frame measurement over an exposure and source catalog.

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog to be filled with the results of measurement. Must contain
            all the `lsst.afw.table.SourceRecord`\ s to be measured (with
            `lsst.afw.detection.Footprint`\ s attached), and have a schema
            that is a superset of ``self.schema``.
        exposure : `lsst.afw.image.ExposureF`
            Image containing the pixel data to be measured together with
            associated PSF, WCS, etc.
        noiseImage : `lsst.afw.image.ImageF`, optional
            Can be used to specify the a predictable noise replacement field
            for testing purposes.
        exposureId : `int`, optional
            Unique exposure identifier used to calculate the random number
            generator seed during noise replacement.
        beginOrder : `float`, optional
            Start execution order (inclusive): measurements with
            ``executionOrder < beginOrder`` are not executed. `None` for no
            limit.
        endOrder : `float`, optional
            Final execution order (exclusive): measurements with
            ``executionOrder >= endOrder`` are not executed. `None` for no
            limit.
        footprints : `dict` {`int`: `lsst.afw.detection.Footprint`}, optional
            List of footprints to use for noise replacement. If this is not
            supplied then the footprints from the measCat are used.
        """
        assert measCat.getSchema().contains(self.schema)
        if footprints is None:
            footprints = self.getFootprintsFromCatalog(measCat)

        # noiseReplacer is used to fill the footprints with noise and save
        # heavy footprints of the source pixels so that they can be restored
        # one at a time for measurement.  After the NoiseReplacer is
        # constructed, all pixels in the exposure.getMaskedImage() which
        # belong to objects in measCat will be replaced with noise

        # Initialize the noise replacer
        noiseReplacer = self.initNoiseReplacer(exposure, measCat, footprints, exposureId, noiseImage)

        self.runPlugins(noiseReplacer, measCat, exposure, beginOrder, endOrder)

    def runPlugins(
        self,
        noiseReplacer,
        measCat,
        exposure,
        beginOrder=None,
        endOrder=None,
    ):
        r"""Call the configured measument plugins on an image.

        Parameters
        ----------
        noiseReplacer : `NoiseReplacer`
            Used to fill sources not being measured with noise.
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog to be filled with the results of measurement. Must contain
            all the `lsst.afw.table.SourceRecord`\ s to be measured (with
            `lsst.afw.detection.Footprint`\ s attached), and have a schema
            that is a superset of ``self.schema``.
        exposure : `lsst.afw.image.ExposureF`
            Image containing the pixel data to be measured together with
            associated PSF, WCS, etc.
        beginOrder : `float`, optional
            Start execution order (inclusive): measurements with
            ``executionOrder < beginOrder`` are not executed. `None` for no
            limit.
        endOrder : `float`, optional
            Final execution order (exclusive): measurements with
            ``executionOrder >= endOrder`` are not executed. `None` for no
            limit.
        """
        nMeasCat = len(measCat)
        self.log.info("Measuring %d source%s", nMeasCat, ("" if nMeasCat == 1 else "s"))

        # Wrap the task logger into a period logger
        periodicLog = PeriodicLogger(self.log)

        for recordIndex, measRecord in enumerate(measCat):
            noiseReplacer.insertSource(measRecord.getId())
            self.callMeasure(measRecord, exposure, beginOrder=beginOrder, endOrder=endOrder)

            if self.doBlendedness:
                self.blendPlugin.cpp.measureChildPixels(exposure.getMaskedImage(), measRecord)

            noiseReplacer.removeSource(measRecord.getId())
            # Log a message if it has been a while since the last log.
            periodicLog.log("Measurement complete for %d sources out of %d", recordIndex + 1, nMeasCat)

        # When done, restore the exposure to its original state
        noiseReplacer.end()

        # Undeblended plugins only fire if we're running everything
        if endOrder is None:
            for sourceIndex, source in enumerate(measCat):
                for plugin in self.undeblendedPlugins.iter():
                    self.doMeasurement(plugin, source, exposure)
                # Log a message if it has been a while since the last log.
                periodicLog.log("Undeblended measurement complete for %d sources out of %d",
                                sourceIndex + 1, nMeasCat)

        # Now we loop over all of the sources one more time to compute the
        # blendedness metrics
        if self.doBlendedness:
            for source in measCat:
                self.blendPlugin.cpp.measureParentPixels(exposure.getMaskedImage(), source)

    def measure(self, measCat, exposure):
        """Backwards-compatibility alias for `run`.
        """
        self.run(measCat, exposure)
