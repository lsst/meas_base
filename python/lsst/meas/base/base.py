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

This includes base classes for plugins and tasks and utility classes such as PluginMap
that are shared by the single-frame measurement framework and the forced measurement framework.
"""

import re
import traceback
import collections

import lsst.pipe.base
import lsst.pex.config

from .baseLib import *
from .noiseReplacer import *

# Exceptions that the measurement tasks should always propagate up to their callers
FATAL_EXCEPTIONS = (MemoryError, FatalAlgorithmError)

#
# Constant dict of the predefined slots and their types.
#
# This essentially duplicates information in afw::table::SourceTable, but we don't have a way to pull
# it out of there right now, and it's best to wait to address that until we can remove the
# meas_algorithms measurement system and simplify SourceTable at the same time.
SLOT_TYPES = {
    "Centroid": "Centroid",
    "Shape": "Shape",
    "PsfFlux", "Flux",
    "ModelFlux", "Flux",
    "ApFlux", "Flux",
    "InstFlux", "Flux",
}

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

    # A dict of index values that can be used to fill slots, and the types these
    # correspond to.
    SLOT_CHOICES = {}

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

    


class BaseMeasurementConfig(lsst.pex.config.Config):
    """!
    Base config class for all measurement driver tasks.
    """

    # Regular expression used to match slot definitions
    _SLOT_REGEXP = re.compile("(?P<plugin>\w+)((\[(?P<n>\d+)\])|(\['(?P<s1>\w+)'\])|(\[\"(?P<s2>\w+)\"\]))?")

    def _parseSlot(self, name):
        """Parse a re.Match object corresponding to _SLOT_REGEXP, returning a tuple of
        (<plugin-name>, <index>)
        """
        value = self.slots[name]
        # n.b. regexp enforces that we have at most one of the groups below non-None
        match = self._SLOT_REGEXP.match(value)
        if not match:
            raise lsst.pex.config.FieldValidationError(
                type(self).slots,
                self,
                "Cannot parse %r for %s slot" % (value, name)
                )
        index = (int(match.group("n")) if match.group("n") else None)
        if index is None: index = match.group("s1")
        if index is None: index = match.group("s2")
        return (match.get("plugin"), index)

    slots = lsst.pex.config.DictField(
        keytype=str,
        itemtype=str,
        doc=("Defines a mapping from algorithms to special aliases in Source.\n"
             "The allowed Keys are: Centroid, Shape, PsfFlux, ModelFlux, ApFlux, and\n"
             "InstFlux.  Values are the name of a measurement plugin of the appropriate\n"
             "type, optionally followed by an algorithm-dependent numeric or string\n"
             "index in square brackets, e.g.:\n"
             "\n"
             "measurement.slots['Centroid'] = 'base_SdssCentroid'\n"
             "measurement.slots['ApFlux'] = 'base_CircularApertureFlux[2]\n"
             "measurement.slots['ModelFlux'] = 'modelfit_CModel['exp']\n"
             "\n"
             "Valid indices for a particular plugin are listed in the SLOT_CHOICES\n"
             "sequence in the plugin's config class."),
        default={
            "Centroid": "base_SdssCentroid",
            "Shape": "base_SdssShape",
            "PsfFlux": "base_PsfFlux",
            "ModelFlux": "base_GaussianFlux",
            "ApFlux": "base_SincFlux",
            "InstFlux": "base_GaussianFlux",
            }
        )

    doReplaceWithNoise = lsst.pex.config.Field(dtype=bool, default=True, optional=False,
        doc='When measuring, replace other detected footprints with noise?')

    noiseReplacer = lsst.pex.config.ConfigField(
        dtype=NoiseReplacerConfig,
        doc="configuration that sets how to replace neighboring sources with noise"
        )

    def validate(self):
        lsst.pex.config.Config.validate(self)
        for slotName in self.slots:
            if slotName not in self._SLOT_TYPES:
                raise lsst.pex.config.FieldValidationError(
                    type(self).slots,
                    self,
                    "Unsupported slot: %r" % slotName
                    )
            plugin, index = self._parseSlot(slotName)
            if plugin not in self.plugins.names:
                raise lsst.pex.config.FieldValidationError(
                    type(self).slots,
                    self,
                    "Plugin %r for %s slot is not being run" % (plugin, slotName)
                    )
            if index not in self.plugins[plugin].SLOT_CHOICES:
                raise lsst.pex.config.FieldValidationError(
                    type(self).slots,,
                    self,
                    "%r is not a valid index for %r in slot %s" % (index, plugin, slotName)
                    )


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

    def initializePlugins(self, **kwds):
        # Transform slots config dict values by parsing strings into Structs of plugin name and index,
        # then reversing the mapping so the plugins are keys, and each value is a list of Structs of
        # (slot, index).  We do this so we can initialize the slots in the order we initialize
        # plugins, so dependent plugins can expect to see the slots they need already setup by the
        # time they're initialized.
        slotsByPlugin = {}
        for slotName in self.config.slots:
            plugin, index = self._parseSlot(slotName)
            slotsByPlugin.setdefault(plugin, []).append((slotName, index))
        # Make a place at the beginning for the centroid plugin to run first (because it's an OrderedDict,
        # adding an empty element in advance means it will get run first when it's reassigned to the
        # actual Plugin).
        if "Centroid" in self.config.slots:
            self.plugins[self.config.slots["Centroid"].plugin] = None
        # Initialize the plugins, sorted by execution order.  At the same time add fields to the schema
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            plugin = PluginClass(config, name, metadata=self.algMetadata, **kwds)
            self.plugins[name] = plugin
            for slotName, index in slotsByPlugin.get(name, []):
                plugin.setupSlot(slotName, 


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
