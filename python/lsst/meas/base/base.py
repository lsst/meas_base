#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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

import collections

import lsst.pipe.base
import lsst.pex.config

# Set of names of algorithms that measure fluxes that can be aperture corrected
_ApCorrNameSet = set()

def addApCorrName(name):
    """!Add to the set of field name prefixes for fluxes that should be aperture corrected

    @param[in] name  field name prefix for a flux that should be aperture corrected.
        The corresponding field names are {name}_flux, {name}_fluxSigma and {name}_flag.
        For example name "base_PsfFlux" corresponds to fields base_PsfFlux_flux,
        base_PsfFlux_fluxSigma and base_PsfFlux_flag.
    """
    global _ApCorrNameSet
    _ApCorrNameSet.add(str(name))

def getApCorrNameSet():
    """!Return a copy of the set of field name prefixes for fluxes that should be aperture corrected

    For example the returned set will likely include "base_PsfFlux" and "base_GaussianFlux".
    """
    global _ApCorrNameSet
    return _ApCorrNameSet.copy()

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

class PluginRegistry(lsst.pex.config.Registry):
    """!
    Base class for plugin registries

    The Plugin class allowed in the registry is defined in the ctor of the registry.

    Single-frame and forced plugins have different registries.
    """

    class Configurable(object):
        """!
        Class used as the actual element in the registry

        Rather than constructing a Plugin instance, its __call__ method
        (invoked by RegistryField.apply) returns a tuple
        of (executionOrder, name, config, PluginClass), which can then
        be sorted before the plugins are instantiated.
        """

        __slots__ = "PluginClass", "name"

        def __init__(self, name, PluginClass):
            """!
            Create a Configurable object for the given PluginClass and name
            """
            self.name = name
            self.PluginClass = PluginClass

        @property
        def ConfigClass(self): return self.PluginClass.ConfigClass

        def __call__(self, config):
            return (self.PluginClass.getExecutionOrder(), self.name, config, self.PluginClass)

    def register(self, name, PluginClass, shouldApCorr=False, apCorrList=()):
        """!
        Register a Plugin class with the given name.

        The same Plugin may be registered multiple times with different names; this can
        be useful if we often want to run it multiple times with different configuration.

        @param[in] name  name of plugin class. This is used as a prefix for all fields produced by the Plugin,
            and it should generally contain the name of the Plugin or Algorithm class itself
            as well as enough of the namespace to make it clear where to find the code.
            For example "base_GaussianFlux" indicates an algorithm in meas_base
            that measures Gaussian Flux and produces fields such as "base_GaussianFlux_flux",
            "base_GaussianFlux_fluxSigma" and "base_GaussianFlux_flag".
        @param[in] shouldApCorr  if True then this algorithm measures a flux that should be aperture
            corrected. This is shorthand for apCorrList=[name] and is ignored if apCorrList is specified.
        @param[in] apCorrList  list of field name prefixes for flux fields that should be aperture corrected.
            If an algorithm produces a single flux that should be aperture corrected then it is simpler
            to set shouldApCorr=True. But if an algorithm produces multiple such fields then it must
            specify apCorrList, instead. For example modelfit_CModel produces 3 such fields:
                apCorrList=("modelfit_CModel_exp", "modelfit_CModel_exp", "modelfit_CModel_def")
            If apCorrList is non-empty then shouldApCorr is ignored.
        """
        lsst.pex.config.Registry.register(self, name, self.Configurable(name, PluginClass))
        if shouldApCorr and not apCorrList:
            apCorrList = [name]
        for prefix in apCorrList:
            addApCorrName(prefix)

    def makeField(self, doc, default=None, optional=False, multi=False):
        return lsst.pex.config.RegistryField(doc, self, default, optional, multi)


def register(name, shouldApCorr=False):
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
    @endcode
    """
    def decorate(PluginClass):
        PluginClass.registry.register(name, PluginClass, shouldApCorr=shouldApCorr)
        return PluginClass
    return decorate


class PluginMap(collections.OrderedDict):
    """!
    Map of plugins (instances of subclasses of BasePlugin) to be run for a task

    We assume plugins are added to the PluginMap according to their "Execution Order", so this
    class doesn't actually do any of the sorting (though it does have to maintain that order,
    which it does by inheriting from OrderedDict).
    """

    def iter(self):
        """!Return an iterator over plugins for which plugin.config.doMeasure is true

        @note plugin.config.doMeasure is usually a simple boolean class attribute, not a normal Config field.
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasure:
                yield plugin

    def iterN(self):
        """!Return an iterator over plugins for which plugin.config.doMeasureN is true

        @note plugin.config.doMeasureN is usually a simple boolean class attribute, not a normal Config field.
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasureN:
                yield plugin
