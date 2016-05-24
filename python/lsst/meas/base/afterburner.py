from collections import namedtuple

import lsst.pipe.base
import lsst.pex.config
import lsst.daf.base

from .pluginsBase import BasePlugin, BasePluginConfig
from .pluginRegistry import PluginRegistry
from .baseLib import FatalAlgorithmError, MeasurementError

# Exceptions that the measurement tasks should always propagate up to their callers
FATAL_EXCEPTIONS = (MemoryError, FatalAlgorithmError)

__all__ = ("AfterburnerPluginConfig", "AfterburnerPlugin", "AfterburnerConfig", "AfterburnerTask")


class AfterburnerPluginConfig(BasePluginConfig):
    '''
    Default configuration class for afterburner plugins
    '''
    pass


class AfterburnerPlugin(BasePlugin):
    '''
    Base class for after burner plugin
    '''
    registry = PluginRegistry(AfterburnerPluginConfig)
    ConfigClass = AfterburnerPluginConfig
    # This defines if the plugin operates on a single source at a time, or expects the whole catalog.
    # The value defaults to single for a single source, set to multi when the plugin expects the whole
    # catalog. If The plugin is of type multi, the fail method should be implemented to accept the whole
    # catalog. If the plugin is of type the fail method should accept a single source record.

    plugType = 'single'

    def __init__(self, config, name, schema, metadata):
        """!
        Initialize the afterburner plugin

        @param[in] config     An instance of afterburner config class.
        @param[in] name       The string the plugin was registered with.
        @param[in,out] schema The source schema, New fields should be added here to
                              hold output produced by this plugin.
        @param[in] metadata   Plugin metadata that will be attached to the output catalog
        """
        BasePlugin.__init__(self, config, name)

    @classmethod
    def getExecutionOrder(cls):
        ''' Sets the relative order of plugins (smaller numbers run first).

        Afterburner plugins must run with BasePlugin.DEFAULT_AFTERBURNER or higher

        All plugins must implement this method with an appropriate run level
        '''
        raise NotImplementedError()

    def burn(self, cat, **kwargs):
        """!
        Process either a single catalog enter or the whole catalog and produce output defined by the plugin

        @param[in,out] cat  Either a lsst source catalog or a catalog entery depending on the plug type
                            specified in the classes configuration. Results may be added to new columns,
                            or existing entries altered.
        @param[in] kwargs   Any additional kwargs that may be passed through the afterburner task.
        """
        raise NotImplementedError()


class AbContext(object):
    '''
    Context manager to handle catching errors that may have been thrown in an afterburner plugin
    @param[in] plugin   The plugin that is to be run
    @param[in] cat      Either a catalog or a source record entry of a catalog, depending of the plugin type,
                        i.e. either working on a whole catalog, or a single record.
    @param[in] log      The log which to write to, most likely will always be the log (self.log) of the object
                        in which the context manager is used.
    '''
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
            self.log.warn("Error in {}.burn: {}".format(self.plugin.name, exc_value))
        return True


class AfterburnerConfig(lsst.pex.config.Config):
    '''
    Default AfterburnerConfig. Currently this is an empty list, meaning that there are no default plugins run.
    The config object for each plugin must use this variable to specify the names of all plugins to be run.
    '''
    plugins = AfterburnerPlugin.registry.makeField(
            multi=True,
            default=["base_ClassificationExtendedness"],
            doc="Plugins to be run and their configuration")


class AfterburnerTask(lsst.pipe.base.Task):
    '''
    This task facilitates running plugins which will operate on a source catalog. These plugins may do things
    such as classifying an object based on source record entries inserted during a measurement task.

    Plugins may either take an entire catalog to work on at a time, or
    '''
    ConfigClass = AfterburnerConfig
    _DefaultName = "afterburner"

    def __init__(self, schema, plugMetadata=None, **kwargs):
        """
        Constructor; only called by derived classes

        @param[in] plugMetaData     An lsst.daf.base.PropertyList that will be filled with metadata
                                    about the plugins being run. If None, an empty empty PropertyList
                                    will be created.
        @param[in] **kwargs         Additional arguments passed to lsst.pipe.base.Task.__init__.
        """
        lsst.pipe.base.Task.__init__(self, **kwargs)
        self.schema = schema
        if plugMetadata is None:
            plugMetadata = lsst.daf.base.PropertyList()
        self.plugMetadata = plugMetadata

        self.initializePlugins()

    def initializePlugins(self):
        '''
        Initialize the plugins according to the configuration.
        '''

        pluginType = namedtuple('pluginType', 'single multi')
        self.executionDict = {}
        # Read the properties for each plugin. Allocate a dictionary entry for each run level. Verify that
        # the plugins are above the minimum run level for an afterburner plugin. For each run level, the
        # plugins are sorted into either single record, or multi record groups to later be run appropriately
        for executionOrder, name, config, PluginClass in self.config.plugins.apply():
            if executionOrder not in self.executionDict:
                self.executionDict[executionOrder] = pluginType(single=[], multi=[])
            if PluginClass.getExecutionOrder() >= BasePlugin.DEFAULT_AFTERBURNER:
                plug = PluginClass(config, name, self.schema, metadata=self.plugMetadata)
                if plug.plugType == 'single':
                    self.executionDict[executionOrder].single.append(plug)
                elif plug.plugType == 'multi':
                    self.executionDict[executionOrder].multi.append(plug)
            else:
                raise ValueError("{} has an execution order less than the minimum for an afterburner plugin."
                                 "Value {} : Minimum {}".format(PluginClass, PluginClass.getExecutionOrder(),
                                                                BasePlugin.DEFAULT_AFTERBURNER))

    def run(self, measCat):
        '''
        The entry point for the afterburner task. This method should be called with a reference to a
        measurement catalog.
        '''
        self.callCompute(measCat)

    def callCompute(self, catalog):
        '''
        Run each of the plugins on the catalog
        @param[in] catalog  The catalog on which the plugins will operate
        '''
        for runlevel in sorted(self.executionDict):
            # Run all of the plugins which take a whole catalog first
            for plug in self.executionDict[runlevel].multi:
                with AbContext(plug, catalog, self.log):
                    plug.burn(catalog)
            # Run all the plugins which take single catalog entries
            for measRecord in catalog:
                for plug in self.executionDict[runlevel].single:
                    with AbContext(plug, measRecord, self.log):
                        plug.burn(measRecord)
