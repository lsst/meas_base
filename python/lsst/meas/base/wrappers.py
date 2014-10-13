import lsst.pex.config

from .base import generateAlgorithmName
from .sfm import SingleFramePlugin, SingleFramePluginConfig
from .forcedMeasurement import ForcedPlugin, ForcedPluginConfig

def wrapAlgorithmControl(Base, Control, hasMeasureN=False, executionOrder=None):
    """!
    Wrap a C++ algorithm's control class into a Python Config class.

    @param[in] Base              Base class for the returned ConfigClass; one of SingleFramePluginConfig or
                                 ForcedPluginConfig
    @param[in] Control           Control class to be wrapped (a Swigged C++ class)
    @param[in] hasMeasureN       Whether the plugin supports fitting multiple objects at once (if so, a
                                 config option to enable/disable this will be added).
    @param[in] executionOrder    If not None, an override for the default executionOrder for this plugin.

    """
    if hasMeasureN:
        # We need to add a Config field to enable multi-object measurement, to replace
        # the simple bool class attribute that's on the base class.  To do that, we
        # create the Config class dynamically here, then call makeControlClass to finish
        # it off by adding fields from the control object.
        cls = type(
            Control.__name__.replace("Control", "Config"),
            (Base,),
            {"doMeasureN": lsst.pex.config.Field(dtype=bool, default=True,
                                                 doc="whether to run this plugin in multi-object mode")}
            )
        ConfigClass = lsst.pex.config.makeControlClass(Control, module=Control.__module__, cls=cls)
    else:
        # If we don't have to add that Config field, we can delegate all of the work to
        # pex_config's makeControlClass
        ConfigClass = lsst.pex.config.makeControlClass(Control, module=Control.__module__, base=base)
    if executionOrder is not None:
        ConfigClass.executionOrder.default = float(executionOrder)
    return ConfigClass

def wrapAlgorithm(Base, AlgClass, name=None, Control=None, ConfigClass=None, needsMetadata=False,
                  hasMeasureN=False, executionOrder=None, doRegister=True):
    if ConfigClass is None:
        if Control is None:
            Control = AlgClass.Control
        ConfigClass = wrapAlgorithmControl(Base.ConfigClass, Control, hasMeasureN=hasMeasureN,
                                           executionOrder=executionOrder)
    PluginClass = type(AlgClass.__name__ + Base.__name__, (Base,),
                       dict(AlgClass=AlgClass, ConfigClass=ConfigClass))
    if doRegister:
        if name is None:
            name = generateAlgorithmName(AlgClass)
        Base.registry.register(name, PluginClass)
    return PluginClass

class WrappedSingleFramePlugin(SingleFramePlugin):

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.cpp = self.AlgClass(config.makeControl(), name, schema, metadata)

    def measure(self, measRecord, exposure):
        self.cpp.measure(measRecord, exposure):

    def measureN(self, measCat, exposure):
        self.cpp.measureN(measRecord, exposure):

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None):

    wrap = classmethod(wrapAlgorithm)


class WrappedForcedPlugin(ForcedPlugin):

    def __init__(self, config, name, schemaMapper, metadata):
        SingleFramePlugin.__init__(self, config, name, schemaMapper, metadata)
        self.cpp = self.AlgClass(config.makeControl(), name, schemaMapper, metadata)

    def measure(self, measRecord, exposure):
        self.cpp.measure(measRecord, exposure):

    def measureN(self, measCat, exposure):
        self.cpp.measureN(measRecord, exposure):

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None):

    wrap = classmethod(wrapAlgorithm)

