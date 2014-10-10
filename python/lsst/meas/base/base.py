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
"""Base and utility classes for measurement frameworks

This includes base classes for plugins and tasks and utility classes such as PluginMap and
MeasurementDataFlags that are shared by the single-frame measurement framework and the forced
measurement framework.
"""

import traceback
import collections

import lsst.pipe.base
import lsst.pex.config

from .baseLib import *
from .noiseReplacer import *

FATAL_EXCEPTIONS = (MemoryError, FatalAlgorithmError)  # Exceptions that the framework should always propagate up

def generateAlgorithmName(AlgClass):
    """Generate a string name for an algorithm class that strips away terms that are generally redundant
    while (hopefully) remaining easy to trace to the code.

    The returned name will cobmine the package name, with any "lsst" and/or "meas" prefix removed,
    with the class name, with any "Algorithm" suffix removed.  For instance,
    lsst.meas.base.SdssShapeAlgorithm becomes "base_SdssShape".
    """
    name = AlgClass.__name__
    pkg = AlgClass.__module__
    name = name.replace("Algorithm", "")
    terms = pkg.split(".")
    if terms[-1].endswith("Lib"):
        terms = terms[:-1]
    if terms[0] == "lsst":
        terms = terms[1:]
    if terms[0] == "meas":
        terms = terms[1:]
    if name.lower().startswith(terms[-1].lower()):
        terms = terms[:-1]
    return "%s_%s" % ("_".join(terms), name)

# Translation map from new PixelFlags to old ones defined in meas_algorithms

_flagMap = {
    "base_PixelFlags_flag_bad":"flags.pixel.bad",
    "base_PixelFlags_flag_edge":"flags.pixel.edge",
    "base_PixelFlags_flag_interpolated":"flags.pixel.interpolated.any",
    "base_PixelFlags_flag_saturated":"flags.pixel.saturated.any",
    "base_PixelFlags_flag_cr":"flags.pixel.cr.any",
    "base_PixelFlags_flag_interpolatedCenter":"flags.pixel.interpolated.center",
    "base_PixelFlags_flag_saturatedCenter":"flags.pixel.saturated.center",
    "base_PixelFlags_flag_crCenter":"flags.pixel.cr.center",
    }

def Version0FlagMapper(flags):
    _flags = []
    for name in flags:
        if name in _flagMap.keys():
            _flags.append(_flagMap[name])
        else:
            _flags.append(name)
    return _flags

def addDependencyFlagAliases(AlgClass, name, schema):
    """!
    Add aliases to flag fields that an algorithm depends on.

    When an algorithm relies on the slot centroid or shape as an input, it can fail (or be unreliable)
    simply because the algorithm it depends on failed.  In that case, we already have a flag field
    that indicates why the algorithm failed (the slot failure flag), but it's not obviously connected
    to the failure of the dependent algorithm.  This function adds aliases to the appropriate slot flag
    for C++ algorithms that depend on that slot, in the namespace of the dependent algorithm.

    @param[in]  AlgClass    Swig-wrappped C++ algorithm class; must have a class attribute Input that
                            constains the type of Input object it requires.
    @param[in]  name        Name of the dependent algorithm
    @param[out] schema      Schema to add the aliases to

    """
    if issubclass(AlgClass.Input, FootprintCentroidInput):
        schema.getAliasMap().set("%s_flag_badCentroid" % name, "slot_Centroid_flag")
    if issubclass(AlgClass.Input, FootprintCentroidShapeInput):
        schema.getAliasMap().set("%s_flag_badShape" % name, "slot_Shape_flag")


class PluginRegistry(lsst.pex.config.Registry):
    """!
    Base class for plugin registries

    The Plugin class allowed in the registry is defined on the ctor of the registry
    Single-frame and forced plugins have different registries.
    """

    class Configurable(object):
        """!
        Class used as the actual element in the registry

        Rather than constructing a Plugin instance, it returns a tuple
        of (runlevel, name, config, PluginClass), which can then
        be sorted before the plugins are instantiated.
        """

        __slots__ = "PluginClass", "name"

        def __init__(self, name, PluginClass):
            """!
            Initialize registry with Plugin Class
            """
            self.name = name
            self.PluginClass = PluginClass

        @property
        def ConfigClass(self): return self.PluginClass.ConfigClass

        def __call__(self, config):
            return (config.executionOrder, self.name, config, self.PluginClass)

    def register(self, name, PluginClass):
        """!
        Register a Plugin class with the given name.

        The same Plugin may be registered multiple times with different names; this can
        be useful if we often want to run it multiple times with different configuration.

        The name will be used as a prefix for all fields produced by the Plugin, and it
        should generally contain the name of the Plugin or Algorithm class itself
        as well as enough of the namespace to make it clear where to find the code.
        """
        lsst.pex.config.Registry.register(self, name, self.Configurable(name, PluginClass))

    def makeField(self, doc, default=None, optional=False, multi=False):
        return lsst.pex.config.RegistryField(doc, self, default, optional, multi)


def register(name):
    """!
    A Python decorator that registers a class, using the given name, in its base class's PluginRegistry.
    For example,
    @code
    @register("base_TransformedCentroid")
    class ForcedTransformedCentroidPlugin(ForcedPlugin):
        ...
    @endcode
    is equivalent to:
    @code
    class ForcedTransformedCentroidPlugin(ForcedPlugin):
        ...
    @ForcedPlugin.registry.register("base_TransformedCentroid", ForcedTransformedCentroidPlugin)
    @enedcode
    """
    def decorate(PluginClass):
        PluginClass.registry.register(name, PluginClass)
        return PluginClass
    return decorate


class PluginMap(collections.OrderedDict):
    """!
    Map of plugins to be run for a task

    We assume Plugins are added to the PluginMap according to their executionOrder, so this
    class doesn't actually do any of the sorting (though it does have to maintain that order).
    """

    def iter(self):
        """Call each plugin in the map which has a measure() method
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasure:
                yield plugin

    def iterN(self):
        """Call each plugin in the map which has a measureN() method
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasureN:
                yield plugin

class BasePluginConfig(lsst.pex.config.Config):
    """!
    Base class measurement Plugin config classes.

    Most derived classes will want to override setDefaults() in order to customize
    the default exceutionOrder.

    A derived class whose corresponding Plugin class implements measureN() should
    additionally add a bool doMeasureN field to replace the bool class attribute
    defined here.
    """

    executionOrder = lsst.pex.config.Field(
        dtype=float, default=2.0,
        doc="""Sets the relative order of plugins (smaller numbers run first).

In general, the following values should be used (intermediate values
are also allowed, but should be avoided unless they are needed):
   0.0 ------ centroids and other algorithms that require only a Footprint and
              its Peaks as input
   1.0 ------ shape measurements and other algorithms that require
              getCentroid() to return a good centroid in addition to a
              Footprint and its Peaks.
   2.0 ------ flux algorithms that require both getShape() and getCentroid()
              in addition to the Footprint and its Peaks
   3.0 ------ Corrections applied to fluxes (i.e. aperture corrections, tying
              model to PSF fluxes). All flux measurements should have an
              executionOrder < 3.0, while all algorithms that rely on corrected
              fluxes (i.e. classification) should have executionOrder > 3.0.
"""
        )


    doMeasure = lsst.pex.config.Field(dtype=bool, default=True,
                                      doc="whether to run this plugin in single-object mode")

    doMeasureN = False  # replace this class attribute with a Field if measureN-capable

class BasePlugin(object):
    """!
    Base class for measurement plugins.

    This is the base class for SingleFramePlugin and ForcedPlugin; derived classes should inherit
    from one of those.
    """

    def fail(self, measRecord, error=None):
        """!
        Record a failure of the measure or measureN() method.

        When measure() raises an exception, the measurement framework
        will call fail() to allow the plugin to set its failure flag
        field(s).  When measureN() raises an exception, fail() will be
        called repeatedly with all the records that were being
        measured.

        If the exception is a MeasurementError, it will be passed as
        the error argument; in all other cases the error argument will
        be None, and the failure will be logged by the measurement
        framework as a warning.
        """
        traceback.print_exc()
        message = ("The algorithm '%s' thinks it cannot fail, but it did; "
                   "please report this as a bug (the full traceback is above)."
                   % self.__class__.__name__)
        raise NotImplementedError(message)

class MeasurementDataFlags(object):
    """!
    Flags that describe data to be measured, allowing plugins with the same signature but
    different requirements to assert their appropriateness for the data they will be run on.

    This is just a placeholder for now; see DM-1026.
    """

    PRECONVOLVED = 0x01  # the image has already been convolved with its PSF;
                         # set for preconvolved difference images and Kaiser-style coadds

    DIFFERENCE = 0x02    # the image is a difference, and hence may have negative surface brightnesses

    COADD = 0x04  # the image is a coadd

    NO_PSF = 0x08 # the image has no Psf object attached

    NO_WCS = 0x10 # the image has no Wcs object attached

class SourceSlotConfig(lsst.pex.config.Config):
    """!
    Slot configuration which assigns a particular named plugin to each of a set of
    slots.  Each slot allows a type of measurement to be fetched from the SourceTable
    without knowing which algorithm was used to produced the data.

    NOTE: the default algorithm for each slot must be registered, even if the default is not used.
    """

    centroid = lsst.pex.config.Field(dtype=str, default="base_SdssCentroid", optional=True,
                                     doc="the name of the centroiding algorithm used to set source x,y")
    shape = lsst.pex.config.Field(dtype=str, default="base_SdssShape", optional=True,
                                  doc="the name of the algorithm used to set source moments parameters")
    apFlux = lsst.pex.config.Field(dtype=str, default="base_SincFlux", optional=True,
                                   doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = lsst.pex.config.Field(dtype=str, default="base_GaussianFlux", optional=True,
                                      doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = lsst.pex.config.Field(dtype=str, default="base_PsfFlux", optional=True,
                                    doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = lsst.pex.config.Field(dtype=str, default="base_GaussianFlux", optional=True,
                                     doc="the name of the algorithm used to set the source inst flux slot")

    def setupTable(self, table, prefix=None):
        """Convenience method to setup a table's slots according to the config definition.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        if prefix is None: prefix = ""
        if self.centroid is not None: table.defineCentroid(prefix + self.centroid)
        if self.shape is not None: table.defineShape(prefix + self.shape)
        if self.apFlux is not None: table.defineApFlux(prefix + self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(prefix + self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(prefix + self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(prefix + self.instFlux)


class BaseMeasurementConfig(lsst.pex.config.Config):
    """!
    Base config class for all measurement driver tasks.
    """

    slots = lsst.pex.config.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source."
        )

    doReplaceWithNoise = lsst.pex.config.Field(dtype=bool, default=True, optional=False,
        doc='When measuring, replace other detected footprints with noise?')

    noiseReplacer = lsst.pex.config.ConfigField(
        dtype=NoiseReplacerConfig,
        doc="configuration that sets how to replace neighboring sources with noise"
        )

    def validate(self):
        lsst.pex.config.Config.validate(self)
        if self.slots.centroid is not None and self.slots.centroid not in self.plugins.names:
            raise ValueError("source centroid slot algorithm is not being run.")
        if self.slots.shape is not None and self.slots.shape not in self.plugins.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None:
                for name in self.plugins.names:
                    if len(name) <= len(slot) and name == slot[:len(name)]:
                        break
                else:
                    raise ValueError("source flux slot algorithm '%s' is not being run." % slot)


## @addtogroup LSST_task_documentation
## @{
## @page baseMeasurementTask
## BaseMeasurementTask @copybrief BaseMeasurementTask
## @}

class BaseMeasurementTask(lsst.pipe.base.Task):
    """!
    Ultimate base class for all measurement tasks.

    This base class for SingleFrameMeasurementTask and ForcedMeasurementTask mostly exists to share
    code between the two, and generally should not be used directly.
    """

    ConfigClass = BaseMeasurementConfig
    _DefaultName = "measurement"
    tableVersion = 1

    def __init__(self, algMetadata=None, **kwds):
        """!
        Constructor; only called by derived classes.

        @param[in]  algMetadata     An lsst.daf.base.PropertyList that will be filled with metadata
                                    about the plugins being run.  If None, an empty PropertyList will
                                    be created.
        @param[in]  **kwds          Additional arguments passed to lsst.pipe.base.Task.__init__.

        This attaches two public attributes to the class for use by derived classes and parent tasks:
         - plugins: an empty PluginMap, which will eventually contain all active plugins that will by
           invoked by the run() method (to be filled by subclasses).  This should be considered read-only.
         - algMetadata: a lsst.daf.base.PropertyList that will contain additional information about the
           active plugins to be saved with the output catalog (to be filled by subclasses).
        """
        super(BaseMeasurementTask, self).__init__(**kwds)
        self.plugins = PluginMap()
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        self.algMetadata = algMetadata

    def callMeasure(self, measRecord, *args, **kwds):
        """!
        Call the measure() method on all plugins, handling exceptions in a consistent way.

        @param[in,out]  measRecord     lsst.afw.table.SourceRecord that corresponds to the object being
                                       measured, and where outputs should be written.
        @param[in]      *args          Positional arguments forwarded to Plugin.measure()
        @param[in]      **kwds         Keyword arguments forwarded to Plugin.measure()

        This method can be used with plugins that have different signatures; the only requirement is that
        'measRecord' be the first argument.  Subsequent positional arguments and keyword arguments are
        forwarded directly to the plugin.

        This method should be considered "protected"; it is intended for use by derived classes, not users.
        """
        for plugin in self.plugins.iter():
            try:
                plugin.measure(measRecord, *args, **kwds)
            except FATAL_EXCEPTIONS:
                raise
            except MeasurementError as error:
                plugin.fail(measRecord, error)
            except Exception as error:
                self.log.warn("Error in %s.measure on record %s: %s"
                              % (plugin.name, measRecord.getId(), error))
                plugin.fail(measRecord)

    def callMeasureN(self, measCat, *args, **kwds):
        """!
        Call the measureN() method on all plugins, handling exceptions in a consistent way.

        @param[in,out]  measCat        lsst.afw.table.SourceCatalog containing records for just
                                       the source family to be measured, and where outputs should
                                       be written.
        @param[in]      *args          Positional arguments forwarded to Plugin.measure()
        @param[in]      **kwds         Keyword arguments forwarded to Plugin.measure()

        This method can be used with plugins that have different signatures; the only requirement is that
        'measRecord' be the first argument.  Subsequent positional arguments and keyword arguments are
        forwarded directly to the plugin.

        This method should be considered "protected"; it is intended for use by derived classes, not users.
        """
        for plugin in self.plugins.iterN():
            try:
                plugin.measureN(measCat, *args, **kwds)
            except FATAL_EXCEPTIONS:
                raise
            except MeasurementError as error:
                for measRecord in measCat:
                    plugin.fail(measRecord, error)
            except Exception as error:
                for measRecord in measCat:
                    plugin.fail(measRecord)
                self.log.warn("Error in %s.measureN on records %s-%s: %s"
                              % (plugin.name, measCat[0].getId(), measCat[-1].getId(), error))
