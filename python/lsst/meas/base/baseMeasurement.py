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

"""Base measurement task, which subclassed by the single frame and forced
measurement tasks.
"""

import lsst.pipe.base
import lsst.pex.config

from .pluginRegistry import PluginMap
from .exceptions import FatalAlgorithmError, MeasurementError
from .pluginsBase import BasePluginConfig, BasePlugin
from .noiseReplacer import NoiseReplacerConfig

__all__ = ("BaseMeasurementPluginConfig", "BaseMeasurementPlugin",
           "BaseMeasurementConfig", "BaseMeasurementTask")

# Exceptions that the measurement tasks should always propagate up to their
# callers
FATAL_EXCEPTIONS = (MemoryError, FatalAlgorithmError)


class BaseMeasurementPluginConfig(BasePluginConfig):
    """Base config class for all measurement plugins.

    Notes
    -----
    Most derived classes will want to override `setDefaults` in order to
    customize the default `executionOrder`.

    A derived class whose corresponding Plugin class implements a do `measureN`
    method should additionally add a bool `doMeasureN` field to replace the
    bool class attribute defined here.
    """

    doMeasure = lsst.pex.config.Field(dtype=bool, default=True,
                                      doc="whether to run this plugin in single-object mode")

    doMeasureN = False  # replace this class attribute with a Field if measureN-capable


class BaseMeasurementPlugin(BasePlugin):
    """Base class for all measurement plugins.

    Notes
    -----
    This is class is a placeholder for future behavior which will be shared
    only between measurement plugins and is implemented for symmetry with the
    measurement base plugin configuration class
    """

    pass


class SourceSlotConfig(lsst.pex.config.Config):
    """Assign named plugins to measurement slots.

    Slot configuration which assigns a particular named plugin to each of a set
    of slots.  Each slot allows a type of measurement to be fetched from the
    `lsst.afw.table.SourceTable` without knowing which algorithm was used to
    produced the data.

    Notes
    -----
    The default algorithm for each slot must be registered, even if the default
    is not used.
    """

    Field = lsst.pex.config.Field
    centroid = Field(dtype=str, default="base_SdssCentroid", optional=True,
                     doc="the name of the centroiding algorithm used to set source x,y")
    shape = Field(dtype=str, default="base_SdssShape", optional=True,
                  doc="the name of the algorithm used to set source moments parameters")
    psfShape = Field(dtype=str, default="base_SdssShape_psf", optional=True,
                     doc="the name of the algorithm used to set PSF moments parameters")
    apFlux = Field(dtype=str, default="base_CircularApertureFlux_12_0", optional=True,
                   doc="the name of the algorithm used to set the source aperture instFlux slot")
    modelFlux = Field(dtype=str, default="base_GaussianFlux", optional=True,
                      doc="the name of the algorithm used to set the source model instFlux slot")
    psfFlux = Field(dtype=str, default="base_PsfFlux", optional=True,
                    doc="the name of the algorithm used to set the source psf instFlux slot")
    gaussianFlux = Field(dtype=str, default="base_GaussianFlux", optional=True,
                         doc="the name of the algorithm used to set the source Gaussian instFlux slot")
    calibFlux = Field(dtype=str, default="base_CircularApertureFlux_12_0", optional=True,
                      doc="the name of the instFlux measurement algorithm used for calibration")

    def setupSchema(self, schema):
        """Set up a slots in a schema following configuration directives.

        Parameters
        ----------
        schema : `lsst.afw.table.Schema`
            The schema in which slots will be set up.

        Notes
        -----
        This is defined in this configuration class to support use in unit
        tests without needing to construct an `lsst.pipe.base.Task` object.
        """
        aliases = schema.getAliasMap()
        if self.centroid is not None:
            aliases.set("slot_Centroid", self.centroid)
        if self.shape is not None:
            aliases.set("slot_Shape", self.shape)
        if self.psfShape is not None:
            aliases.set("slot_PsfShape", self.psfShape)
        if self.apFlux is not None:
            aliases.set("slot_ApFlux", self.apFlux)
        if self.modelFlux is not None:
            aliases.set("slot_ModelFlux", self.modelFlux)
        if self.psfFlux is not None:
            aliases.set("slot_PsfFlux", self.psfFlux)
        if self.gaussianFlux is not None:
            aliases.set("slot_GaussianFlux", self.gaussianFlux)
        if self.calibFlux is not None:
            aliases.set("slot_CalibFlux", self.calibFlux)


class BaseMeasurementConfig(lsst.pex.config.Config):
    """Base configuration for all measurement driver tasks.

    Examples
    --------
    Subclasses should define the 'plugins' and 'undeblended' registries, e.g.

    .. code-block:: py

        plugins = PluginBaseClass.registry.makeField(
            multi=True,
            default=[],
            doc="Plugins to be run and their configuration"
        )
        undeblended = PluginBaseClass.registry.makeField(
            multi=True,
            default=[],
            doc="Plugins to run on undeblended image"
        )

    where ``PluginBaseClass`` is the appropriate base class of the plugin
    (e.g., `SingleFramePlugin` or `ForcedPlugin`).
    """

    slots = lsst.pex.config.ConfigField(
        dtype=SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source."
    )

    doReplaceWithNoise = lsst.pex.config.Field(
        dtype=bool, default=True, optional=False,
        doc='When measuring, replace other detected footprints with noise?')

    noiseReplacer = lsst.pex.config.ConfigField(
        dtype=NoiseReplacerConfig,
        doc="configuration that sets how to replace neighboring sources with noise"
    )
    undeblendedPrefix = lsst.pex.config.Field(
        dtype=str, default="undeblended_",
        doc="Prefix to give undeblended plugins"
    )

    def validate(self):
        lsst.pex.config.Config.validate(self)
        if self.slots.centroid is not None and self.slots.centroid not in self.plugins.names:
            raise ValueError("source centroid slot algorithm is not being run.")
        if self.slots.shape is not None and self.slots.shape not in self.plugins.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux,
                     self.slots.gaussianFlux, self.slots.calibFlux):
            if slot is not None:
                for name in self.plugins.names:
                    if len(name) <= len(slot) and name == slot[:len(name)]:
                        break
                else:
                    raise ValueError("source instFlux slot algorithm '%s' is not being run." % slot)


class BaseMeasurementTask(lsst.pipe.base.Task):
    """Ultimate base class for all measurement tasks.

    Parameters
    ----------
    algMetadata : `lsst.daf.base.PropertyList` or `None`
        Will be modified in-place to contain metadata about the plugins being
        run. If `None`, an empty `~lsst.daf.base.PropertyList` will be
        created.
    **kwds
        Additional arguments passed to `lsst.pipe.base.Task.__init__`.

    Notes
    -----
    This base class for `SingleFrameMeasurementTask` and
    `ForcedMeasurementTask` mostly exists to share code between the two, and
    generally should not be used directly.
    """

    ConfigClass = BaseMeasurementConfig
    _DefaultName = "measurement"

    plugins = None
    """Plugins to be invoked (`PluginMap`).

    Initially empty, this will be populated as plugins are initialized. It
    should be considered read-only.
    """

    algMetadata = None
    """Metadata about active plugins (`lsst.daf.base.PropertyList`).

    Contains additional information about active plugins to be saved with
    the output catalog. Will be filled by subclasses.
    """

    def __init__(self, algMetadata=None, **kwds):
        super(BaseMeasurementTask, self).__init__(**kwds)
        self.plugins = PluginMap()
        self.undeblendedPlugins = PluginMap()
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        self.algMetadata = algMetadata

    def initializePlugins(self, **kwds):
        """Initialize plugins (and slots) according to configuration.

        Parameters
        ----------
        **kwds
            Keyword arguments forwarded directly to plugin constructors.

        Notes
        -----
        Derived class constructors should call this method to fill the
        `plugins` attribute and add corresponding output fields and slot
        aliases to the output schema.

        In addition to the attributes added by `BaseMeasurementTask.__init__`,
        a ``schema``` attribute holding the output schema must be present
        before this method is called.

        Keyword arguments are forwarded directly to plugin constructors,
        allowing derived classes to use plugins with different signatures.
        """
        # Make a place at the beginning for the centroid plugin to run first
        # (because it's an OrderedDict, adding an empty element in advance
        # means it will get run first when it's reassigned to the actual
        # Plugin).
        if self.config.slots.centroid is not None:
            self.plugins[self.config.slots.centroid] = None
        # Init the plugins, sorted by execution order.  At the same time add to
        # the schema
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            #   Pass logName to the plugin if the plugin is marked as using it
            #   The task will use this name to log plugin errors, regardless.
            if hasattr(PluginClass, "hasLogName") and PluginClass.hasLogName:
                self.plugins[name] = PluginClass(config, name, metadata=self.algMetadata,
                                                 logName=self.log.getChild(name).name, **kwds)
            else:
                self.plugins[name] = PluginClass(config, name, metadata=self.algMetadata, **kwds)

        # In rare circumstances (usually tests), the centroid slot not be
        # coming from an algorithm, which means we'll have added something we
        # don't want to the plugins map, and we should remove it.
        if self.config.slots.centroid is not None and self.plugins[self.config.slots.centroid] is None:
            del self.plugins[self.config.slots.centroid]
        # Initialize the plugins to run on the undeblended image
        for executionOrder, name, config, PluginClass in sorted(self.config.undeblended.apply()):
            undeblendedName = self.config.undeblendedPrefix + name
            self.undeblendedPlugins[name] = PluginClass(config, undeblendedName, metadata=self.algMetadata,
                                                        **kwds)

    def callMeasure(self, measRecord, *args, **kwds):
        """Call ``measure`` on all plugins and consistently handle exceptions.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            The record corresponding to the object being measured. Will be
            updated in-place with the results of measurement.
        *args
            Positional arguments forwarded to ``plugin.measure``
        **kwds
            Keyword arguments. Two are handled locally:

            beginOrder : `int`
                Beginning execution order (inclusive). Measurements with
                ``executionOrder`` < ``beginOrder`` are not executed. `None`
                for no limit.

            endOrder : `int`
                Ending execution order (exclusive). Measurements with
                ``executionOrder`` >= ``endOrder`` are not executed. `None`
                for no limit.

            Others are forwarded to ``plugin.measure()``.

        Notes
        -----
        This method can be used with plugins that have different signatures;
        the only requirement is that ``measRecord`` be the first argument.
        Subsequent positional arguments and keyword arguments are forwarded
        directly to the plugin.

        This method should be considered "protected": it is intended for use by
        derived classes, not users.
        """
        beginOrder = kwds.pop("beginOrder", None)
        endOrder = kwds.pop("endOrder", None)
        for plugin in self.plugins.iter():
            if beginOrder is not None and plugin.getExecutionOrder() < beginOrder:
                continue
            if endOrder is not None and plugin.getExecutionOrder() >= endOrder:
                break
            self.doMeasurement(plugin, measRecord, *args, **kwds)

    def doMeasurement(self, plugin, measRecord, *args, **kwds):
        """Call ``measure`` on the specified plugin.

        Exceptions are handled in a consistent way.

        Parameters
        ----------
        plugin : subclass of `BasePlugin`
            Plugin that will be executed.
        measRecord : `lsst.afw.table.SourceRecord`
            The record corresponding to the object being measured. Will be
            updated in-place with the results of measurement.
        *args
            Positional arguments forwarded to ``plugin.measure()``.
        **kwds
            Keyword arguments forwarded to ``plugin.measure()``.

        Notes
        -----
        This method can be used with plugins that have different signatures;
        the only requirement is that ``plugin`` and ``measRecord`` be the first
        two arguments.  Subsequent positional arguments and keyword arguments
        are forwarded directly to the plugin.

        This method should be considered "protected": it is intended for use by
        derived classes, not users.
        """
        try:
            plugin.measure(measRecord, *args, **kwds)
        except FATAL_EXCEPTIONS:
            raise
        except MeasurementError as error:
            self.log.getChild(plugin.name).debug(
                "MeasurementError in %s.measure on record %s: %s",
                plugin.name, measRecord.getId(), error)
            plugin.fail(measRecord, error)
        except Exception as error:
            self.log.getChild(plugin.name).debug(
                "Exception in %s.measure on record %s: %s",
                plugin.name, measRecord.getId(), error)
            plugin.fail(measRecord)

    def callMeasureN(self, measCat, *args, **kwds):
        """Call ``measureN`` on all plugins and consistently handle exceptions.

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog containing only the records for the source family to be
            measured, and where outputs should be written.
        *args
            Positional arguments forwarded to ``plugin.measure()``
        **kwds
            Keyword arguments. Two are handled locally:

            beginOrder:
                Beginning execution order (inclusive): Measurements with
                ``executionOrder`` < ``beginOrder`` are not executed. `None`
                for no limit.
            endOrder:
                Ending execution order (exclusive): measurements with
                ``executionOrder`` >= ``endOrder`` are not executed. `None` for
                no ``limit``.

            Others are are forwarded to ``plugin.measure()``.

        Notes
        -----
        This method can be used with plugins that have different signatures;
        the only requirement is that ``measRecord`` be the first argument.
        Subsequent positional arguments and keyword arguments are forwarded
        directly to the plugin.

        This method should be considered "protected": it is intended for use by
        derived classes, not users.
        """
        beginOrder = kwds.pop("beginOrder", None)
        endOrder = kwds.pop("endOrder", None)
        for plugin in self.plugins.iterN():
            if beginOrder is not None and plugin.getExecutionOrder() < beginOrder:
                continue
            if endOrder is not None and plugin.getExecutionOrder() >= endOrder:
                break
            self.doMeasurementN(plugin, measCat, *args, **kwds)

    def doMeasurementN(self, plugin, measCat, *args, **kwds):
        """Call ``measureN`` on the specified plugin.

        Exceptions are handled in a consistent way.

        Parameters
        ----------
        plugin : subclass of `BasePlugin`
            Plugin that will be executed.
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog containing only the records for the source family to be
            measured, and where outputs should be written.
        *args
            Positional arguments forwarded to ``plugin.measureN()``.
        **kwds
            Keyword arguments forwarded to ``plugin.measureN()``.

        Notes
        -----
        This method can be used with plugins that have different signatures;
        the only requirement is that the ``plugin`` and ``measCat`` be the
        first two arguments. Subsequent positional arguments and keyword
        arguments are forwarded directly to the plugin.

        This method should be considered "protected": it is intended for use by
        derived classes, not users.
        """
        try:
            plugin.measureN(measCat, *args, **kwds)
        except FATAL_EXCEPTIONS:
            raise

        except MeasurementError as error:
            for measRecord in measCat:
                self.log.getChild(plugin.name).debug(
                    "MeasurementError in %s.measureN on records %s-%s: %s",
                    plugin.name, measCat[0].getId(), measCat[-1].getId(), error)
                plugin.fail(measRecord, error)
        except Exception as error:
            for measRecord in measCat:
                plugin.fail(measRecord)
                self.log.getChild(plugin.name).debug(
                    "Exception in %s.measureN on records %s-%s: %s",
                    plugin.name, measCat[0].getId(), measCat[-1].getId(), error)
