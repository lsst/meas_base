import lsst.pex.config
from .pluginsBase import BasePlugin
from .pluginRegistry import generateAlgorithmName, register
from .apCorrRegistry import addApCorrName
from .sfm import SingleFramePlugin, SingleFramePluginConfig
from .forcedMeasurement import ForcedPlugin, ForcedPluginConfig

__all__ = ("wrapSingleFrameAlgorithm", "wrapForcedAlgorithm", "wrapSimpleAlgorithm",
           "wrapAlgorithm", "wrapAlgorithmControl", "wrapTransform", "GenericPlugin")


class WrappedSingleFramePlugin(SingleFramePlugin):

    def __init__(self, config, name, schema, metadata, logName=None):
        SingleFramePlugin.__init__(self, config, name, schema, metadata, logName=logName)
        if hasattr(self, "hasLogName") and self.hasLogName and logName is not None:
            self.cpp = self.factory(config, name, schema, metadata, logName=logName)
        else:
            self.cpp = self.factory(config, name, schema, metadata)

    def measure(self, measRecord, exposure):
        self.cpp.measure(measRecord, exposure)

    def measureN(self, measCat, exposure):
        self.cpp.measureN(measCat, exposure)

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None)


class WrappedForcedPlugin(ForcedPlugin):

    def __init__(self, config, name, schemaMapper, metadata, logName=None):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata, logName=logName)
        if hasattr(self, "hasLogName") and self.hasLogName and logName is not None:
            self.cpp = self.factory(config, name, schemaMapper, metadata, logName=logName)
        else:
            self.cpp = self.factory(config, name, schemaMapper, metadata)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        self.cpp.measureForced(measRecord, exposure, refRecord, refWcs)

    def measureN(self, measCat, exposure, refCat, refWcs):
        self.cpp.measureNForced(measCat, exposure, refCat, refWcs)

    def fail(self, measRecord, error=None):
        self.cpp.fail(measRecord, error.cpp if error is not None else None)


def wrapAlgorithmControl(Base, Control, module=2, hasMeasureN=False):
    """Wrap a C++ algorithm's control class into a Python config class.

    Parameters
    ----------
    Base : `SingleFramePluginConfig` or `ForcedPluginConfig`
        Base class for the returned config.
    Control : pybind11-wrapped version of a C++ class.
        Control class to be wrapped.
    module : module, `str` or `int`; optional
        Either a module object, a string specifying the name of the module, or
        an integer specifying how far back in the stack to look for the module
        to use: ``0`` is `lsst.pex.config.wrap`, ``1`` is
        `lsst.meas.base.wrappers`, ``2`` is the immediate caller, etc.  This
        will be used to set ``__module__`` for the new config class, and the
        class will also be added to the module.  The default is to use the
        callers' module.
    hasMeasureN : `bool`, optional
        Whether the plugin supports fitting multiple objects at once (if so, a
        config option to enable/disable this will be added).

    Returns
    -------
    ConfigClass : `lsst.pex.config.Config`
        A new subclass of lsst.pex.config.Config.

    Notes
    -----
    This function is generally only called by `wrapAlgorithm`; it is unlikely
    users will have to call it directly.
    """
    if hasMeasureN:
        # We need to add a Config field to enable multi-object measurement, to
        # replace the simple bool class attribute that's on the base class.
        # To do that, we create the Config class dynamically here, then call
        # makeControlClass to finish it off by adding fields from the control
        # object.
        cls = type(
            Control.__name__.replace("Control", "Config"),
            (Base,),
            {"doMeasureN": lsst.pex.config.Field(dtype=bool, default=True,
                                                 doc="whether to run this plugin in multi-object mode")}
        )
        ConfigClass = lsst.pex.config.makeConfigClass(Control, module=module, cls=cls)
    else:
        # If we don't have to add that Config field, we can delegate all of
        # the work to pex_config's makeControlClass
        ConfigClass = lsst.pex.config.makeConfigClass(Control, module=module, base=Base)
    return ConfigClass


def wrapAlgorithm(Base, AlgClass, factory, executionOrder, name=None, Control=None,
                  ConfigClass=None, TransformClass=None, doRegister=True, shouldApCorr=False,
                  apCorrList=(), hasLogName=False, **kwds):
    """Wrap a C++ algorithm class to create a measurement plugin.

    Parameters
    ----------
    Base : `SingleFramePlugin` or `ForcedPlugin`
        Base class for the returned Plugin.
    AlgClass : API compatible with `SingleFrameAlgorithm` or `ForcedAlgorithm`
        C++ algorithm class to convert. May either derive directly from
        `SingleFrameAlgorithm` or `ForcedAlgorithm`, or be an unrelated class
        which has the same ``measure`` and ``measureN`` signatures.
    factory : callable
        A callable that is used to construct an instance of ``AlgClass``.  It
        must take four arguments, either ``(config, name, schema, metadata)``
        or ``(config, name, schemaMapper, metadata)``, depending on whether
        the algorithm is single-frame or forced.
    executionOrder : `float`
        The order this plugin should be run, relative to others
        (see `BasePlugin.getExecutionOrder`).
    name : `str`, optional
        String to use when registering the algorithm. Ignored if
        ``doRegistry=False``, set to ``generateAlgorithmName(AlgClass)`` if
        `None`.
    Control : Pybind11-wrapped version of a C++ class, optional
        Pybind11-wrapped C++ Control class for the algorithm;
        ``AlgClass.Control`` is used if `None`. Ignored if ``ConfigClass``
        is not `None`.
    ConfigClass : subclass of `BaseMeasurementPluginConfig`
        Python config class that wraps the C++ algorithm's pybind11-wrapped
        Control class.  If `None`, `wrapAlgorithmControl` is called to
        generate a Config class using the ``Control`` argument.
    TransformClass : subclass of `MeasurementTransform`, optional
        Transformation which may be used to post-process the results of
        measurement.  If `None`, the default defined by `BasePlugin` is
        used.
    doRegister : `bool`, optional
        If `True` (the default), register the plugin with ``Base``'s
        registry, allowing it to be used by measurement tasks.
    shouldApCorr : `bool`, optional
        Does this algorithm measure an instFlux that can be aperture
        corrected?  This is shorthand for ``apCorrList=[name]`` and is ignored
        if ``apCorrList`` is specified.
    apCorrList : iterable of `str`, optional
        Field name prefixes for instFlux fields to be aperture corrected. If
        an algorithm measures a single instFlux that should be aperture
        corrected, then it is simpler to set ``shouldApCorr=True``. However,
        if an algorithm produces multiple such fields, then specify
        ``apCorrList`` instead.  For example, ``modelfit_CModel`` produces
        three such fields: ``apCorrList= ("modelfit_CModel_exp",
        "modelfit_CModel_exp", "modelfit_CModel_def")`` If ``apCorrList`` is
        not empty then ``shouldApCorr`` is ignored.  If non-empty and
        ``doRegister`` is `True` then the names are added to the set
        retrieved by ``getApCorrNameSet``.
    hasLogName : `bool`, optional
        `True` if the C++ algorithm supports ``logName`` as a constructor
        argument.
    **kwds
        Additional keyword arguments passed to generateAlgorithmControl, which
        may include:

        - ``hasMeasureN``:  Whether the plugin supports fitting multiple
          objects at once ;if so, a config option to enable/disable this will
          be added (`bool`).
        - ``executionOrder``: If not `None`, an override for the default
          execution order for this plugin (the default is ``2.0``, which is
          usually appropriate for fluxes; `bool`).

    Returns
    -------
    PluginClass : subclass of ``Base``
        The new plugin class.
    """
    if ConfigClass is None:
        if Control is None:
            Control = AlgClass.Control
        ConfigClass = wrapAlgorithmControl(Base.ConfigClass, Control, **kwds)

    def getExecutionOrder():
        return executionOrder
    typeDict = dict(AlgClass=AlgClass, ConfigClass=ConfigClass, factory=staticmethod(factory),
                    getExecutionOrder=staticmethod(getExecutionOrder))
    if TransformClass:
        typeDict['getTransformClass'] = staticmethod(lambda: TransformClass)
    PluginClass = type(AlgClass.__name__ + Base.__name__, (Base,), typeDict)
    if doRegister:
        if name is None:
            name = generateAlgorithmName(AlgClass)
        Base.registry.register(name, PluginClass)
        if shouldApCorr:
            addApCorrName(name)
    PluginClass.hasLogName = hasLogName
    return PluginClass


def wrapSingleFrameAlgorithm(AlgClass, executionOrder, name=None, needsMetadata=False, hasMeasureN=False,
                             hasLogName=False, **kwds):
    """Expose a C++ ``SingleFrameAlgorithm`` class as a measurement plugin.

    Parameters
    ----------
    AlgClass : API compatible with `SingleFrameAlgorithm`
        C++ algorithm class to convert. May either derive directly from
        `SingleFrameAlgorithm` or be an unrelated class which has the same
        ``measure``, ``measureN`` and ``fail`` signatures.
    executionOrder : `float`
        The order this plugin should be run, relative to others
        (see `BasePlugin.getExecutionOrder`).
    name : `str`, optional
        Name to use when registering the algorithm. Ignored if
        ``doRegistry=False``; set to ``generateAlgorithmName(AlgClass)`` if
        `None`.
    needsMetadata : `bool`, optional
        Sets whether the ``AlgClass``'s constructor should be passed a
        `~lsst.daf.base.PropertySet` metadata argument.
    hasMeasureN : `bool`, optional
        Does the algorithm support simultaneous measurement of multiple
        sources? If `True`, a `bool` ``doMeasureN`` field will be added to
        the generated config class, and its value will be passed as the last
        argument when calling the ``AlgClass`` constructor.
    hasLogName : `bool`, optional
        `True` if the C++ algorithm supports ``logName`` as a constructor
        argument.
    **kwds
        Additional keyword arguments are passed to the lower-level
        `wrapAlgorithm` and `wrapAlgorithmControl` classes.

    Returns
    -------
    singleFramePlugin : subclass of `SingleFramePlugin`
        The new measurement plugin class.

    Notes
    -----
    The first three arguments to the C++ constructor are expected to be
    ``Control const & ctrl, std::string const & name, Schema & schema``.

    If ``needsMetadata`` is `True`, we also append ``PropertySet & metadata``.

    If ``hasMeasureN`` is `True`, we also append ``bool doMeasureN``.

    If ``hasLogName`` is `True`, we also append ``std::string logName``.

    If more than one of the above is `True`, the metadata ``PropertySet``
    precedes the ``doMeasureN`` ``bool`` and the ``logName`` comes last of the
    three.
    """
    if hasMeasureN:
        if needsMetadata:
            def factory(config, name, schema, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, schema, metadata, config.doMeasureN, **kwargs)
        else:
            def factory(config, name, schema, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, schema, config.doMeasureN, **kwargs)
    else:
        if needsMetadata:
            def factory(config, name, schema, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, schema, metadata, **kwargs)
        else:
            def factory(config, name, schema, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, schema, **kwargs)

    return wrapAlgorithm(WrappedSingleFramePlugin, AlgClass, executionOrder=executionOrder, name=name,
                         factory=factory, hasMeasureN=hasMeasureN, hasLogName=hasLogName, **kwds)


def wrapForcedAlgorithm(AlgClass, executionOrder, name=None, needsMetadata=False,
                        hasMeasureN=False, needsSchemaOnly=False, hasLogName=False, **kwds):
    """Expose a C++ ``ForcedAlgorithm`` class as a measurement plugin.

    Parameters
    ----------
    AlgClass : API compatible with `ForcedAlgorithm`
        C++ algorithm class to convert. May either derive directly from
        `ForcedAlgorithm` or be an unrelated class which has the same
        ``measure``, ``measureN`` and ``fail`` signatures.
    executionOrder : `float`
        The order this plugin should be run, relative to others
        (see `BasePlugin.getExecutionOrder`).
    name : `str`, optional
        Name to use when registering the algorithm. Ignored if
        ``doRegistry=False``; set to ``generateAlgorithmName(AlgClass)`` if
        `None`.
    needsMetadata : `bool`, optional
        Sets whether the ``AlgClass``'s constructor should be passed a
        `~lsst.daf.base.PropertySet` metadata argument.
    hasMeasureN : `bool`, optional
        Does the algorithm support simultaneous measurement of multiple
        sources? If `True`, a `bool` ``doMeasureN`` field will be added to
        the generated config class, and its value will be passed as the last
        argument when calling the ``AlgClass`` constructor.
    hasLogName : `bool`, optional
        `True` if the C++ algorithm supports ``logName`` as a constructor
        argument.
    needsSchemaOnly : `bool`, optional
        Whether the algorithm constructor expects a Schema argument
        (representing the output `~lsst.afw.table.Schema`) rather than the
        full `~lsst.afw.table.SchemaMapper` (which provides access to both the
        reference schema and the output schema).
    **kwds
        Additional keyword arguments are passed to the lower-level
        `wrapAlgorithm` and `wrapAlgorithmControl` classes.

    Returns
    -------
    forcedPlugin : subclass of `ForcedPlugin`
        The new measurement plugin class.

    Notes
    -----
    The first two arguments to the C++ constructor are expected to be
    ``Control const & ctrl, std::string const & name``

    If ``needsSchemaOnly`` is `True`, then the third argument will be
    ``Schema & schema``; otherwise, it will be ``SchemaMapper &
    schemaMapper``.

    If ``needsMetadata`` is `True`, we also append ``PropertySet &
    metadata``.

    If ``hasMeasureN`` is `True`, we also append ``bool doMeasureN``.

    If ``hasLogName`` is `True`, we also append ``std::string logName``.

    If more than one of the above is `True`, the metadata ``PropertySet``
    precedes the ``doMeasureN`` ``bool`` and the ``logName`` comes last of the
    three.
    """
    if needsSchemaOnly:
        def extractSchemaArg(m):
            return m.editOutputSchema()
    else:
        def extractSchemaArg(m):
            return m
    if hasMeasureN:
        if needsMetadata:
            def factory(config, name, schemaMapper, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                metadata, config.doMeasureN, **kwargs)
        else:
            def factory(config, name, schemaMapper, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                config.doMeasureN, **kwargs)
    else:
        if needsMetadata:
            def factory(config, name, schemaMapper, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper),
                                metadata, **kwargs)
        else:
            def factory(config, name, schemaMapper, metadata, **kwargs):
                return AlgClass(config.makeControl(), name, extractSchemaArg(schemaMapper), **kwargs)

    return wrapAlgorithm(WrappedForcedPlugin, AlgClass, executionOrder=executionOrder, name=name,
                         factory=factory, hasLogName=hasLogName, **kwds)


def wrapSimpleAlgorithm(AlgClass, executionOrder, name=None, needsMetadata=False, hasMeasureN=False,
                        hasLogName=False, **kwds):
    r"""Expose a C++ ``SimpleAlgorithm`` class as a measurement plugin.

    ``SimpleAlgorithm``\ s are made available as both `SingleFramePlugin`\ s
    and `ForcedPlugin`\ s.

    Parameters
    ----------
    AlgClass : Subclass of C++ ``SimpleAlgorithm``, or API compatible
        Algorithm class to convert. The C++ class should be wrapped with
        Pybind11, and must provide ``measure()``, ``measureN()`` and ``fail()`
        signatures equivalent to ``SimpleAlgorithm``.
    executionOrder : `float`
        The order this plugin should be run, relative to others
        (see `~BasePlugin.getExecutionOrder`).
    name : `str`, optional
        Name to use when registering the algorithm. Ignored if
        ``doRegistry=False``; set to ``generateAlgorithmName(AlgClass)`` if
        `None`.
    needsMetadata : `bool`, optional
        Sets whether the ``AlgClass``'s constructor should be passed a
        `~lsst.daf.base.PropertySet` metadata argument.
    hasMeasureN : `bool`, optional
        Does the algorithm support simultaneous measurement of multiple
        sources? If `True`, a `bool` ``doMeasureN`` field will be added to
        the generated config class, and its value will be passed as the last
        argument when calling the ``AlgClass`` constructor.
    hasLogName : `bool`, optional
        `True` if the C++ algorithm supports ``logName`` as a constructor
        argument.
    **kwds
        Additional keyword arguments are passed to the lower-level
        `wrapAlgorithm` and `wrapAlgorithmControl` classes.

    Returns
    -------
    singleFramePlugin : subclass of `SingleFramePlugin`
        The new single frame measurement plugin class.
    forcedPlugin : subclass of `ForcedPlugin`
        The new forced measurement plugin class.

    Notes
    -----
    The first three arguments to the C++ constructor are expected to be
    ``Control const & ctrl, std::string const & name, Schema & schema``.

    If ``needsMetadata`` is `True`, we also append ``PropertySet &
    metadata``.

    If ``hasMeasureN`` is `True`, we also append ``bool doMeasureN``.

    If ``hasLogName`` is `True`, we also append ``std::string logName``.

    If more than one of the above is `True`, the metadata ``PropertySet``
    precedes the ``doMeasureN`` ``bool`` and the ``logName`` comes last of the
    three.
    """
    return (wrapSingleFrameAlgorithm(AlgClass, executionOrder=executionOrder, name=name,
                                     needsMetadata=needsMetadata, hasLogName=hasLogName, **kwds),
            wrapForcedAlgorithm(AlgClass, executionOrder=executionOrder, name=name,
                                needsMetadata=needsMetadata, hasLogName=hasLogName,
                                needsSchemaOnly=True, **kwds))


def wrapTransform(transformClass, hasLogName=False):
    """Modify a C++ transform to accept either a ``Config`` or a ``Control``.

    That is, the configuration may either be provided as a (C++) ``Control``
    object or an instance of a Python class derived from
    `~lsst.meas.base.BasePluginConfig`.

    Parameters
    ----------
    transformClass : Subclass of C++ ``BaseTransform``
        A C++ transform class, wrapped with pybind11. Its constructor must
        take a ``Control`` object, a ``std::string``, and a
        `~lsst.afw.table.SchemaMapper`, in that order.
    hasLogName : `bool`, optional
        Unused.
    """
    oldInit = transformClass.__init__

    def _init(self, ctrl, name, mapper, logName=None):
        if hasattr(ctrl, "makeControl"):
            ctrl = ctrl.makeControl()
        # logName signature needs to be on this Class __init__, but is not
        # needed by the C++ plugin.
        oldInit(self, ctrl, name, mapper)

    transformClass.__init__ = _init


class GenericPlugin(BasePlugin):
    """Abstract base class for a generic plugin.

    Parameters
    ----------
    config : `lsst.pex.config.Config`
        An instance of this class' ``ConfigClass``.
    name : `str`
        Name of this measurement plguin, for registering.
    schema : `lsst.afw.table.Schema`
        The catalog schema. New fields should be added here to
        hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Metadata that will be attached to the output catalog.
    logName : `str`, optional
        Name of log component.

    Notes
    -----
    A generic plugin can be used with the `singleFramePluginFromGeneric`
    and/or `forcedPluginFromGeneric` wrappers to create classes that can be
    used for single frame measurement and/or forced measurement (as
    appropriate). The only real difference between `SingleFramePlugin` and
    `ForcedPlugin` is the ``measure`` method; this class introduces a shared
    signature for `measure` that, in combination with the aforementioned
    wrappers, allows both plugin styles to share a single implementation.

    This doesn't use `abc.ABCMeta` because I couldn't get it to work
    with a superclass.

    Sub-classes should set `ConfigClass` and implement the `measure` and
    `measureN` methods. They may optionally provide alternative
    implementations for the `__init__`, `fail` and `getExecutionOrder`
    methods.

    This default implementation simply adds a field for recording
    a fatal failure of the measurement plugin.
    """
    ConfigClass = None

    @classmethod
    def getExecutionOrder(cls):
        return 0

    def __init__(self, config, name, schema, metadata, logName=None):
        BasePlugin.__init__(self, config, name, logName=logName)
        self._failKey = schema.addField(name + '_flag', type="Flag", doc="Set for any fatal failure")

    def measure(self, measRecord, exposure, center):
        """Measure a single source.

        It is the responsibility of this method to perform the desired
        measurement and record the result in the `measRecord`.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Catalog record for the source being measured.
        exposure : `lsst.afw.image.Exposure`
            Exposure on which the source is being measured.
        center : `lsst.geom.Point2D`
            Pixel coordinates of the object.

        Raises
        ------
        MeasurementError
            Raised if the measurement fails for a known/justifiable reason.
        """
        raise NotImplementedError()

    def measureN(self, measCat, exposure, refCat, refWcs):
        """Measure multiple sources.

        It is the responsibility of this method to perform the desired
        measurement and record the result in the `measCat`.

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog for the sources being measured.
        exposure : `lsst.afw.image.Exposure`
            Exposure on which the source is being measured.
        refCat : `lsst.afw.table.SourceCatalog`
            Reference catalog.
        refWcs : `lsst.afw.image.Wcs`
            Astrometric solution for the reference image.

        Raises
        ------
        MeasurementError
            Raised if the measurement fails for a known/justifiable reason.
        """
        raise NotImplementedError()

    def fail(self, measRecord, error=None):
        """Record a measurement failure.

        This default implementation simply records the failure in the source
        record.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Catalog record for the source being measured.
        error : `Exception`
            Error causing failure, or `None`.
        """
        measRecord.set(self._failKey, True)

    @classmethod
    def makeSingleFramePlugin(cls, name):
        """Produce a SingleFramePlugin subclass from this GenericPlugin class.

        The class is also registered.

        Parameters
        ----------
        name : `str`
            Name of plugin to register.
        """
        class SingleFrameFromGenericConfig(cls.ConfigClass, SingleFramePluginConfig):
            pass

        @register(name)
        class SingleFrameFromGenericPlugin(SingleFramePlugin):
            ConfigClass = SingleFrameFromGenericConfig

            def __init__(self, config, name, schema, metadata, logName=None):
                SingleFramePlugin.__init__(self, config, name, schema, metadata, logName=logName)
                self._generic = cls(config, name, schema, metadata)

            def measure(self, measRecord, exposure):
                center = measRecord.getCentroid()
                return self._generic.measure(measRecord, exposure, center)

            def measureN(self, measCat, exposure, refCat, refWcs):
                return self._generic.measureN(measCat, exposure, refCat, refWcs)

            def fail(self, measRecord, error=None):
                self._generic.fail(measRecord, error if error is not None else None)

            @staticmethod
            def getExecutionOrder():
                return cls.getExecutionOrder()

            def getTransformClass(self):
                return self._generic.getTransformClass()

        return SingleFrameFromGenericPlugin

    @classmethod
    def makeForcedPlugin(cls, name):
        """Produce a ForcedPlugin subclass from this GenericPlugin class.

        The class is also registered.

        Parameters
        ----------
        name : `str`
            Name of plugin to register.
        """
        class ForcedFromGenericConfig(cls.ConfigClass, ForcedPluginConfig):
            pass

        @register(name)
        class ForcedFromGenericPlugin(ForcedPlugin):
            ConfigClass = ForcedFromGenericConfig

            def __init__(self, config, name, schemaMapper, metadata, logName=None):
                ForcedPlugin.__init__(self, config, name, schemaMapper, metadata, logName=logName)
                schema = schemaMapper.editOutputSchema()
                self._generic = cls(config, name, schema, metadata)

            def measure(self, measRecord, exposure, refRecord, refWcs):
                center = exposure.getWcs().skyToPixel(refWcs.pixelToSky(refRecord.getCentroid()))
                return self._generic.measure(measRecord, exposure, center)

            def measureN(self, measCat, exposure, refCat, refWcs):
                return self._generic.measureN(measCat, exposure, refCat, refWcs)

            def fail(self, measRecord, error=None):
                self._generic.fail(measRecord, error if error is not None else None)

            @staticmethod
            def getExecutionOrder():
                return cls.getExecutionOrder()

            def getTransformClass(self):
                return self._generic.getTransformClass()

        return ForcedFromGenericPlugin
