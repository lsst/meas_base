import lsst.pex.config

from .base import generateAlgorithmName
from .sfm import SingleFramePlugin, SingleFramePluginConfig
from .forcedMeasurement import ForcedPlugin, ForcedPluginConfig

__all__ = ("wrapSingleFrameAlgorithm", "wrapForcedAlgorithm", "wrapSimpleAlgorithm")


class WrappedSingleFramePlugin(SingleFramePlugin):

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.cpp = self.factory(config, name, schema, metadata)

    def measure(self, measRecord, exposure):
        self.cpp.measure(measRecord, exposure)

    def measureN(self, measCat, exposure):
        self.cpp.measureN(measRecord, exposure)

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None)


class WrappedForcedPlugin(ForcedPlugin):

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        self.cpp = self.factory(config, name, schemaMapper, metadata)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        self.cpp.measure(measRecord, exposure, refRecord, refWcs)

    def measureN(self, measCat, exposure, refCat, refWcs):
        self.cpp.measureN(measRecord, exposure, refCat, refWcs)

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None)


def wrapAlgorithmControl(Base, Control, hasMeasureN=False, executionOrder=None):
    """!
    Wrap a C++ algorithm's control class into a Python Config class.

    @param[in] Base              Base class for the returned ConfigClass; one of SingleFramePluginConfig or
                                 ForcedPluginConfig
    @param[in] Control           Control class to be wrapped (a Swigged C++ class)
    @param[in] hasMeasureN       Whether the plugin supports fitting multiple objects at once (if so, a
                                 config option to enable/disable this will be added).
    @param[in] executionOrder    If not None, an override for the default executionOrder for this plugin.

    @return a new subclass of lsst.pex.config.Config

    This function is generally only called by wrapAlgorithm; it is unlikely users will have to call it
    directly.
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
        ConfigClass = lsst.pex.config.makeConfigClass(Control, module=Control.__module__, cls=cls)
    else:
        # If we don't have to add that Config field, we can delegate all of the work to
        # pex_config's makeControlClass
        ConfigClass = lsst.pex.config.makeConfigClass(Control, module=Control.__module__, base=Base)
    if executionOrder is not None:
        ConfigClass.executionOrder.default = float(executionOrder)
    return ConfigClass


def wrapAlgorithm(Base, AlgClass, factory, name=None, Control=None, ConfigClass=None, doRegister=True,
                  **kwds):
    """!
    Wrap a C++ Algorithm class into a Python Plugin class.

    @param[in] Base            Base class for the returned Plugin; one of SingleFramePlugin or
                               ForcedPlugin
    @param[in] AlgClass        Swigged C++ Algorithm class to convert; must be a subclass of
                               SingleFrameAlgorithm or ForcedAlgorithm (matching the Base argument), or
                               an unrelated class with the same measure() and measureN() signatures as
                               those base classes.
    @param[in] factory         A callable that is used to construct an instance of AlgClass.  It must take
                               four arguments, either (config, name, schema, metadata) or
                               (config, name, schemaMapper, metadata), depending on whether the algorithm is
                               single-frame or forced.
    @param[in] name            String to use when registering the algorithm.  Ignored if doRegistry=False,
                               set to generateAlgorithmName(AlgClass) if None.
    @param[in] Control         Swigged C++ Control class for the algorithm; AlgClass.Control is used if None.
                               Ignored if ConfigClass is not None.
    @param[in] ConfigClass     Python Config class that wraps the C++ Algorithm's swigged Control class.  If
                               None, wrapAlgorithmControl is called to generate a Config class using the
                               Control argument.
    @param[in] doRegister      If True (the default), register the plugin with Base's registry, allowing it
                               to be used by measurement Tasks.
    @param[in] **kwds          Additional keyword arguments passed to generateAlgorithmControl, including:
                               - hasMeasureN:  Whether the plugin supports fitting multiple objects at once
                                 (if so, a config option to enable/disable this will be added).
                               - executionOrder: If not None, an override for the default executionOrder for
                                 this plugin (the default is 2.0, which is usually appropriate for fluxes).

    @return the new Plugin class, a subclass of Base

    This function is generally only called by the public wrapSingleFrameAlgorithm, wrapForcedAlgorithm, and
    wrapSimpleAlgorithm functions; it is unlikely users will have to call it directly.
    """
    if ConfigClass is None:
        if Control is None:
            Control = AlgClass.Control
        ConfigClass = wrapAlgorithmControl(Base.ConfigClass, Control, **kwds)
    PluginClass = type(AlgClass.__name__ + Base.__name__, (Base,),
                       dict(AlgClass=AlgClass, ConfigClass=ConfigClass, factory=factory))
    if doRegister:
        if name is None:
            name = generateAlgorithmName(AlgClass)
        Base.registry.register(name, PluginClass)
    return PluginClass


def wrapSingleFrameAlgorithm(AlgClass, name=None, needsMetadata=False, hasMeasureN=False, **kwds):
    """!
    Wrap a C++ SingleFrameAlgorithm class into a Python SingleFramePlugin class.

    @param[in] AlgClass        Swigged C++ Algorithm class to convert; must be a subclass of
                               SingleFrameAlgorithm, or an unrelated class with the same measure(),
                               measureN(), and fail() signatures.
    @param[in] name            String to use when registering the algorithm.  Ignored if doRegistry=False,
                               set to generateAlgorithmName(AlgClass) if None.
    @param[in] needsMetadata   Sets whether the AlgClass's constructor should be passed a PropertySet
                               metadata argument.
    @param[in] hasMeasureN     Whether the algorithm supports simultaneous measurement of multiple sources.
                               If True, a bool doMeasureN field will be added to the generated Config class,
                               and its value will be passed as the last argument when calling the AlgClass
                               constructor.
    @param[in] **kwds          Additional keyword arguments passed to the lower-level wrapAlgorithm and
                               wrapAlgorithmControl classes.  These include:
                               - Control: Swigged C++ Control class for the algorithm; AlgClass.Control
                                 is used if None. Ignored if ConfigClass is not None.
                               - ConfigClass: Python Config class that wraps the C++ Algorithm's swigged
                                 Control class.  If None, wrapAlgorithmControl is called to generate a
                                 Config class using the Control argument.
                               - doRegister: If True (the default), register the plugin with
                                 SingleFramePlugin.registry, allowing it to be used by
                                 SingleFrameMeasurementTask.
                               - executionOrder: If not None, an override for the default executionOrder for
                                 this plugin (the default is 2.0, which is usually appropriate for fluxes).

    @return the new SingleFramePlugin subclass

    The needsMetadata and hasMeasureN arguments combine to determine the expected constructor signature;
    we always expect the first three arguments to be:
    @verbatim
    Control const & ctrl, std::string const & name, Schema & schema
    @endverbatim
    If needsMetadata, we also append:
    @verbatim
    PropertySet & metadata
    @endverbatim
    If hasMeasureN, we also append:
    @verbatim
    bool doMeasureN
    @endverbatim
    If both are True, the metadata PropertySet precedes the doMeasureN bool.
    """
    if hasMeasureN:
        if needsMetadata:
            def factory(config, name, schema, metadata):
                return AlgClass(config.makeControl(), name, schema, metadata, config.doMeasureN)
        else:
            def factory(config, name, schema, metadata):
                return AlgClass(config.makeControl(), name, schema, config.doMeasureN)
    else:
        if needsMetadata:
            def factory(config, name, schema, metadata):
                return AlgClass(config.makeControl(), name, schema, metadata)
        else:
            def factory(config, name, schema, metadata):
                return AlgClass(config.makeControl(), name, schema)
    return wrapAlgorithm(WrappedSingleFramePlugin, AlgClass, factory=factory,
                         hasMeasureN=hasMeasureN, **kwds)


def wrapForcedAlgorithm(AlgClass, name=None, needsMetadata=False, hasMeasureN=False, needsSchemaOnly=False,
                        **kwds):
    """!
    Wrap a C++ ForcedAlgorithm class into a Python ForcedPlugin class.

    @param[in] AlgClass        Swigged C++ Algorithm class to convert; must be a subclass of
                               ForcedAlgorithm, or an unrelated class with the same measure(), measureN(),
                               and fail() signatures.
    @param[in] name            String to use when registering the algorithm.  Ignored if doRegistry=False,
                               set to generateAlgorithmName(AlgClass) if None.
    @param[in] needsMetadata   Sets whether the AlgClass's constructor should be passed a PropertySet
                               metadata argument.
    @param[in] hasMeasureN     Whether the algorithm supports simultaneous measurement of multiple sources.
                               If True, a bool doMeasureN field will be added to the generated Config class,
                               and its value will be passed as the last argument when calling the AlgClass
                               constructor.
    @param[in] needsSchemaOnly Whether the algorithm constructor expects a Schema argument (representing the
                               output Schema) rather than the full SchemaMapper (which provides access to
                               both the reference Schema and the output Schema).
    @param[in] **kwds          Additional keyword arguments passed to the lower-level wrapAlgorithm and
                               wrapAlgorithmControl classes.  These include:
                               - Control: Swigged C++ Control class for the algorithm; AlgClass.Control
                                 is used if None. Ignored if ConfigClass is not None.
                               - ConfigClass: Python Config class that wraps the C++ Algorithm's swigged
                                 Control class.  If None, wrapAlgorithmControl is called to generate a
                                 Config class using the Control argument.
                               - doRegister: If True (the default), register the plugin with
                                 ForcedPlugin.registry, allowing it to be used by ForcedMeasurementTask.
                               - executionOrder: If not None, an override for the default executionOrder for
                                 this plugin (the default is 2.0, which is usually appropriate for fluxes).

    @return the new ForcedPlugin subclass

    The needsMetadata, hasMeasureN, and needsSchemaOnly arguments combine to determine the expected
    constructor signature; we always expect the first two arguments to be:
    @verbatim
    Control const & ctrl, std::string const & name
    @endverbatim
    If needsSchemaOnly is True, then the third argument will be
    @verbatim
    Schema & schema
    @endverbatim
    otherwise, it will be:
    @verbatim
    SchemaMapper & schemaMapper
    @endverbatim
    If needsMetadata, we also append:
    @verbatim
    PropertySet & metadata
    @endverbatim
    If hasMeasureN, we also append:
    @verbatim
    bool doMeasureN
    @endverbatim
    If both are True, the metadata PropertySet precedes the doMeasureN bool.
    """
    if needsSchemaOnly:
        extractSchemaArg = lambda m: m.editOutputSchema()
    else:
        extractSchemaArg = lambda m: m
    if hasMeasureN:
        if needsMetadata:
            def factory(config, name, schemaMapper, metadata):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                metadata, config.doMeasureN)
        else:
            def factory(config, name, schemaMapper, metadata):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                config.doMeasureN)
    else:
        if needsMetadata:
            def factory(config, name, schemaMapper, metadata):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                metadata)
        else:
            def factory(config, name, schemaMapper, metadata):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper))
    return wrapAlgorithm(WrappedForcedPlugin, AlgClass, factory=factory, **kwds)


def wrapSimpleAlgorithm(AlgClass, name=None, needsMetadata=False, hasMeasureN=False, **kwds):
    """!
    Wrap a C++ SimpleAlgorithm class into both a Python SingleFramePlugin and ForcedPlugin classes

    @param[in] AlgClass        Swigged C++ Algorithm class to convert; must be a subclass of
                               simpleAlgorithm, or an unrelated class with the same measure(), measureN(),
                               and fail() signatures.
    @param[in] name            String to use when registering the algorithm.  Ignored if doRegistry=False,
                               set to generateAlgorithmName(AlgClass) if None.
    @param[in] needsMetadata   Sets whether the AlgClass's constructor should be passed a PropertySet
                               metadata argument.
    @param[in] hasMeasureN     Whether the algorithm supports simultaneous measurement of multiple sources.
                               If True, a bool doMeasureN field will be added to the generated Config class,
                               and its value will be passed as the last argument when calling the AlgClass
                               constructor.
    @param[in] **kwds          Additional keyword arguments passed to the lower-level wrapAlgorithm and
                               wrapAlgorithmControl classes.  These include:
                               - Control: Swigged C++ Control class for the algorithm; AlgClass.Control
                                 is used if None. Ignored if ConfigClass is not None.
                               - ConfigClass: Python Config class that wraps the C++ Algorithm's swigged
                                 Control class.  If None, wrapAlgorithmControl is called to generate a
                                 Config class using the Control argument.
                               - doRegister: If True (the default), register the plugins with Base's
                                 registry, allowing it to be used by measurement Tasks.
                               - executionOrder: If not None, an override for the default executionOrder for
                                 this plugin (the default is 2.0, which is usually appropriate for fluxes).

    @return a two-element tuple, containing the new SingleFramePlugin and ForcedPlugin subclasses

    The needsMetadata and hasMeasureN arguments combine to determine the expected constructor signature;
    we always expect the first three arguments to be:
    @verbatim
    Control const & ctrl, std::string const & name, Schema & schema
    @endverbatim
    If needsMetadata, we also append:
    @verbatim
    PropertySet & metadata
    @endverbatim
    If hasMeasureN, we also append:
    @verbatim
    bool doMeasureN
    @endverbatim
    If both are True, the metadata PropertySet precedes the doMeasureN bool.
    """
    return (wrapSingleFrameAlgorithm(AlgClass, name=name, needsMetadata=needsMetadata, **kwds),
            wrapForcedAlgorithm(AlgClass, name=name, needsMetadata=needsMetadata,
                                needsSchemaOnly=True, **kwds))
