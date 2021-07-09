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

from collections import namedtuple

import lsst.pipe.base
import lsst.pex.config
import lsst.daf.base

from .pluginsBase import BasePlugin, BasePluginConfig
from .pluginRegistry import PluginRegistry, PluginMap
from . import FatalAlgorithmError, MeasurementError

# Exceptions that the measurement tasks should always propagate up to their
# callers
FATAL_EXCEPTIONS = (MemoryError, FatalAlgorithmError)

__all__ = ("CatalogCalculationPluginConfig", "CatalogCalculationPlugin", "CatalogCalculationConfig",
           "CatalogCalculationTask")


class CatalogCalculationPluginConfig(BasePluginConfig):
    """Default configuration class for catalog calcuation plugins.
    """
    pass


class CatalogCalculationPlugin(BasePlugin):
    """Base class for catalog calculation plugins.

    Parameters
    ----------
    config : `CatalogCalculationPlugin.ConfigClass`
        Plugin configuration.
    name : `str`
        The string the plugin was registered with.
    schema : `lsst.afw.table.Schema`
        The source schema, New fields should be added here to
        hold output produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = CatalogCalculationPluginConfig  # documentation inherited

    registry = PluginRegistry(CatalogCalculationPluginConfig)
    """List of available plugins (`lsst.meas.base.PluginRegistry`).
    """

    plugType = 'single'
    """Does the plugin operate on a single source or the whole catalog (`str`)?

    If the plugin operates on a single source at a time, this should be set to
    ``"single"``; if it expects the whoe catalog, to ``"multi"``.  If the
    plugin is of type ``"multi"``, the `fail` method must be implemented to
    accept the whole catalog. If the plugin is of type ``"single"``, `fail`
    should accept a single source record.
    """

    def __init__(self, config, name, schema, metadata):
        BasePlugin.__init__(self, config, name)

    @classmethod
    def getExecutionOrder(cls):
        r"""Used to set the relative order of plugin execution.

        The values returned by `getExecutionOrder` are compared across all
        plugins, and smaller numbers run first.

        Notes
        -----
        `CatalogCalculationPlugin`\s must run with
        `BasePlugin.DEFAULT_CATALOGCALCULATION` or higher.

        All plugins must implement this method with an appropriate run level
        """
        raise NotImplementedError()

    def calculate(self, cat, **kwargs):
        """Perform the calculation specified by this plugin.

        This method can either be used to operate on a single catalog record
        or a whole catalog, populating it with the output defined by this
        plugin.

        Note that results may be added to catalog records as new columns, or
        may result in changes to existing values.

        Parameters
        ----------
        cat : `lsst.afw.table.SourceCatalog` or `lsst.afw.table.SourceRecord`
            May either be a `~lsst.afw.table.SourceCatalog` or a single
            `~lsst.afw.table.SourceRecord`, depending on the plugin type. Will
            be updated in place to contain the results of plugin execution.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        raise NotImplementedError()


class CCContext:
    """Handle errors that are thrown by catalog calculation plugins.

    This is a context manager.

    Parameters
    ----------
    plugin : `CatalogCalculationPlugin`
        The plugin that is to be run.
    cat : `lsst.afw.table.SourceCatalog` or `lsst.afw.table.SourceRecord`
        May either be a `~lsst.afw.table.SourceCatalog` or a single
        `~lsst.afw.table.SourceRecord`, depending on the plugin type.
    log : `lsst.log.Log`
        A logger. Generally, this should be the logger of the object in which
        the context manager is being used.
    """
    def __init__(self, plugin, cat, log):
        self.plugin = plugin
        self.cat = cat
        self.log = log

    def __enter__(self):
        return

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            return True
        if exc_type in FATAL_EXCEPTIONS:
            raise exc_value
        elif exc_type is MeasurementError:
            self.plugin.fail(self.cat, exc_value)
        else:
            self.log.warning("Error in %s.calculate: %s", self.plugin.name, exc_value)
        return True


class CatalogCalculationConfig(lsst.pex.config.Config):
    """Config class for the catalog calculation driver task.

    Specifies which plugins will execute when the `CatalogCalculationTask`
    associated with this configuration is run.
    """

    plugins = CatalogCalculationPlugin.registry.makeField(
        multi=True,
        default=["base_ClassificationExtendedness",
                 "base_FootprintArea"],
        doc="Plugins to be run and their configuration")


class CatalogCalculationTask(lsst.pipe.base.Task):
    """Run plugins which operate on a catalog of sources.

    This task facilitates running plugins which will operate on a source
    catalog. These plugins may do things such as classifying an object based
    on source record entries inserted during a measurement task.

    Parameters
    ----------
    plugMetaData : `lsst.daf.base.PropertyList` or `None`
        Will be modified in-place to contain metadata about the plugins being
        run. If `None`, an empty `~lsst.daf.base.PropertyList` will be
        created.
    **kwargs
        Additional arguments passed to the superclass constructor.

    Notes
    -----
    Plugins may either take an entire catalog to work on at a time, or work on
    individual records.
    """
    ConfigClass = CatalogCalculationConfig
    _DefaultName = "catalogCalculation"

    def __init__(self, schema, plugMetadata=None, **kwargs):
        lsst.pipe.base.Task.__init__(self, **kwargs)
        self.schema = schema
        if plugMetadata is None:
            plugMetadata = lsst.daf.base.PropertyList()
        self.plugMetadata = plugMetadata
        self.plugins = PluginMap()

        self.initializePlugins()

    def initializePlugins(self):
        """Initialize the plugins according to the configuration.
        """

        pluginType = namedtuple('pluginType', 'single multi')
        self.executionDict = {}
        # Read the properties for each plugin. Allocate a dictionary entry for each run level. Verify that
        # the plugins are above the minimum run level for an catalogCalculation plugin. For each run level,
        # the plugins are sorted into either single record, or multi record groups to later be run
        # appropriately
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            if executionOrder not in self.executionDict:
                self.executionDict[executionOrder] = pluginType(single=[], multi=[])
            if PluginClass.getExecutionOrder() >= BasePlugin.DEFAULT_CATALOGCALCULATION:
                plug = PluginClass(config, name, self.schema, metadata=self.plugMetadata)
                self.plugins[name] = plug
                if plug.plugType == 'single':
                    self.executionDict[executionOrder].single.append(plug)
                elif plug.plugType == 'multi':
                    self.executionDict[executionOrder].multi.append(plug)
            else:
                errorTuple = (PluginClass, PluginClass.getExecutionOrder(),
                              BasePlugin.DEFAULT_CATALOGCALCULATION)
                raise ValueError("{} has an execution order less than the minimum for an catalogCalculation "
                                 "plugin. Value {} : Minimum {}".format(*errorTuple))

    @lsst.pipe.base.timeMethod
    def run(self, measCat):
        """The entry point for the catalog calculation task.

        Parameters
        ----------
        meascat : `lsst.afw.table.SourceCatalog`
            Catalog for measurement.
        """
        self.callCompute(measCat)

    def callCompute(self, catalog):
        """Run each of the plugins on the catalog.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            The catalog on which the plugins will operate.
        """
        for runlevel in sorted(self.executionDict):
            # Run all of the plugins which take a whole catalog first
            for plug in self.executionDict[runlevel].multi:
                with CCContext(plug, catalog, self.log):
                    plug.calculate(catalog)
            # Run all the plugins which take single catalog entries
            for measRecord in catalog:
                for plug in self.executionDict[runlevel].single:
                    with CCContext(plug, measRecord, self.log):
                        plug.calculate(measRecord)
