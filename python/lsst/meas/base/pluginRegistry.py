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

"""Registry for measurement plugins and utilities for plugin management.
"""

import collections

import lsst.pipe.base
import lsst.pex.config
from .apCorrRegistry import addApCorrName

__all__ = ("generateAlgorithmName", "PluginRegistry", "register", "PluginMap")


def generateAlgorithmName(AlgClass):
    """Generate a name for an algorithm.

    This generates a short name for an algorithmic class that strips away
    terms that are generally redundant while remaining easy to trace to the
    code.

    Parameters
    ----------
    AlgClass : subclass of `BaseAlgorithm`
        The class to generate a name for.

    Returns
    -------
    name : `str`
        A short name for the algorithm.

    Notes
    -----
    The returned name will cobmine the package name, with any ``lsst`` and/or
    ``meas`` prefix removed, with the class name, with any ``Algorithm``
    suffix removed.  For instance, ``lsst.meas.base.SdssShapeAlgorithm``
    becomes ``base_SdssShape``.
    """
    name = AlgClass.__name__
    pkg = AlgClass.__module__
    name = name.replace("Algorithm", "")
    terms = pkg.split(".")
    # Hide private module name only if it's part of a public package
    if len(terms) > 1 and terms[-1].startswith("_"):
        terms = terms[:-1]
    if len(terms) > 1 and terms[-1].endswith("Lib"):
        terms = terms[:-1]
    if terms[0] == "lsst":
        terms = terms[1:]
    if terms[0] == "meas":
        terms = terms[1:]
    if name.lower().startswith(terms[-1].lower()):
        terms = terms[:-1]
    return "%s_%s" % ("_".join(terms), name)


class PluginRegistry(lsst.pex.config.Registry):
    """Base class for plugin registries.

    Notes
    -----
    The class of plugins allowed in the registry is defined in the constructor
    of the registry.

    Single-frame and forced plugins have different registries.
    """

    class Configurable:
        """Class used as the element in the plugin registry.

        Parameters
        ----------
        name : `str`
            Name under which the plugin is registerd.
        PluginClass : subclass of `BasePlugin`
            The class of plugin which can be stored in the registry.

        Notes
        -----
        Rather than constructing a Plugin instance, its __call__ method
        (invoked by RegistryField.apply) returns a tuple
        of ``(executionOrder, name, config, PluginClass)``, which can then
        be sorted before the plugins are instantiated.
        """

        __slots__ = "PluginClass", "name"

        def __init__(self, name, PluginClass):
            self.name = name
            self.PluginClass = PluginClass

        @property
        def ConfigClass(self):
            return self.PluginClass.ConfigClass

        def __call__(self, config):
            return (self.PluginClass.getExecutionOrder(), self.name, config, self.PluginClass)

    def register(self, name, PluginClass, shouldApCorr=False, apCorrList=()):
        """Register a plugin class with the given name.

        Parameters
        ----------
        name : `str`
            The name of the plugin. This is used as a prefix for all fields
            produced by the plugin, and it should generally contain the name
            of the plugin or algorithm class itself as well as enough of the
            namespace to make it clear where to find the code.  For example
            ``base_GaussianFlux`` indicates an algorithm in `lsst.meas.base`
            that measures Gaussian Flux and produces fields such as
            ``base_GaussianFlux_instFlux``, ``base_GaussianFlux_instFluxErr``
            and ``base_GaussianFlux_flag``.
        shouldApCorr : `bool`
            If `True`, then this algorithm measures an instFlux that should
            be aperture corrected. This is shorthand for ``apCorrList=[name]``
            and is ignored if ``apCorrList`` is specified.
        apCorrList : `list` of `str`
            List of field name prefixes for instFlux fields to be aperture
            corrected.  If an algorithm produces a single instFlux that should
            be aperture corrected then it is simpler to set
            ``shouldApCorr=True``. But if an algorithm produces multiple such
            fields then it must specify ``apCorrList`` instead. For example,
            ``modelfit_CModel`` produces three such fields:
            ``apCorrList=("modelfit_CModel_exp", "modelfit_CModel_exp",
            "modelfit_CModel_def")``. If ``apCorrList`` is not empty then
            shouldApCorr is ignored.

        Notes
        -----
        The same plugin may be registered multiple times with different names;
        this can be useful if we often want to run it multiple times with
        different configuration.
        """
        lsst.pex.config.Registry.register(self, name, self.Configurable(name, PluginClass))
        if shouldApCorr and not apCorrList:
            apCorrList = [name]
        for prefix in apCorrList:
            addApCorrName(prefix)

    def makeField(self, doc, default=None, optional=False, multi=False):
        return lsst.pex.config.RegistryField(doc, self, default, optional, multi)


def register(name, shouldApCorr=False, apCorrList=()):
    """A decorator to register a plugin class in its base class's registry.

    Parameters
    ----------
    shouldApCorr : `bool`
        If `True`, then this algorithm measures an instFlux that should be
        aperture corrected. This is shorthand for ``apCorrList=[name]`` and is
        ignored if ``apCorrList`` is specified.
    apCorrList : `list` of `str`
        List of field name prefixes for instFlux fields to be aperture
        corrected.  If an algorithm produces a single instFlux that should be
        aperture corrected then it is simpler to set ``shouldApCorr=True``.
        But if an algorithm produces multiple such fields then it must specify
        ``apCorrList`` instead. For example, ``modelfit_CModel`` produces
        three such fields: ``apCorrList=("modelfit_CModel_exp",
        "modelfit_CModel_exp", "modelfit_CModel_def")``. If ``apCorrList`` is
        not empty then shouldApCorr is ignored.

    """

    def decorate(PluginClass):
        PluginClass.registry.register(name, PluginClass, shouldApCorr=shouldApCorr, apCorrList=apCorrList)
        return PluginClass
    return decorate


class PluginMap(collections.OrderedDict):
    """Map of plugins to be run for a given task.

    Notes
    -----
    Plugins are classes derived from `BasePlugin`.

    We assume plugins are added to the plugin map according to their
    "Execution Order", so this class doesn't actually do any of the sorting
    (though it does have to maintain that order, which it does by inheriting
    from `collections.OrderedDict`).
    """

    def iter(self):
        """Return an iterator over plugins for use in single-object mode.

        Notes
        -----
        Plugins which should be used in single-object mode are identified by
        having the `doMeasure` config attribute evaluate to `True`. This is
        usually a simple boolean class attribute.
        """
        for plugin in self.values():
            if plugin.config.doMeasure:
                yield plugin

    def iterN(self):
        """Return an iterator over plugins for use in multi-object mode.

        Notes
        -----
        Plugins which should be used in multi-object mode are identified by
        having the `doMeasureN` config attribute evaluate to `True`.
        This is usually a simple boolean class attribute.
        """
        for plugin in self.values():
            if plugin.config.doMeasureN:
                yield plugin
