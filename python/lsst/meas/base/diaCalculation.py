# This file is part of ap_association.
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
import numpy as np
import pandas as pd

from .catalogCalculation import (CatalogCalculationPluginConfig,
                                 CatalogCalculationPlugin,
                                 CatalogCalculationConfig,
                                 CatalogCalculationTask,
                                 CCContext)
from .pluginsBase import BasePlugin
from .pluginRegistry import (PluginRegistry, PluginMap)
import lsst.pipe.base

# Enforce an error for unsafe column/array value setting in pandas.
pd.options.mode.chained_assignment = 'raise'

__all__ = ("DiaObjectCalculationPlugin", "DiaObjectCalculationPluginConfig",
           "DiaObjectCalculationTask", "DiaObjectCalculationConfig")


class DiaObjectCalculationPluginConfig(CatalogCalculationPluginConfig):
    """Default configuration class for DIA catalog calculation plugins.
    """
    pass


class DiaObjectCalculationPlugin(CatalogCalculationPlugin):
    """Base class for DIA catalog calculation plugins.

    Task follows CatalogCalculationPlugin with modifications for use in AP.

    Parameters
    ----------
    config : `DiaObjectCalculationPlugin.ConfigClass`
        Plugin configuration.
    name : `str`
        The string the plugin was registered with.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = DiaObjectCalculationPluginConfig

    registry = PluginRegistry(DiaObjectCalculationPluginConfig)
    """List of available plugins (`lsst.meas.base.PluginRegistry`).
    """

    FLUX_MOMENTS_CALCULATED = 5.0
    """Add order after flux means and stds are calculated.
    """

    plugType = 'single'
    """Does the plugin operate on a single source or the whole catalog (`str`)?
    If the plugin operates on a single source at a time, this should be set to
    ``"single"``; if it expects the whoe catalog, to ``"multi"``.  If the
    plugin is of type ``"multi"``, the `fail` method must be implemented to
    accept the whole catalog. If the plugin is of type ``"single"``, `fail`
    should accept a single source record.
    """

    inputCols = []
    """DiaObject column names required by the plugin in order to run and
    complete its calculation. DiaCalculationTask should raise an error is a
    plugin is instantiated without the needed column available. Input columns
    should be defined in the DPDD/cat/Apdb schema. Filter dependent columns
    should be specified without the filter name perpended to them. eg
    ``PSFluxMean`` instead of ``uPSFluxMean``.
    """
    outputCols = []
    """DiaObject column names output by the plugin. DiaCalculationTask should
    raise an error if another pluging is run output to the same column.
    Output columns should be defined in the DPDD/cat/Apdb schema. Filter
    dependent columns should be specified without the filter name perpended to
    them. eg ``PSFluxMean`` instead of ``uPSFluxMean``.
    """

    needsFilter = True
    """This plugin requires a filter to be specified. Plugin's using filter
    names usually deal with fluxes and only a sub-set of the DiaSource
    catalog. Plugins that to not use the filter name usually run over a value
    common across all observations/detections such as position.
    """

    def __init__(self, config, name, metadata):
        BasePlugin.__init__(self, config, name)

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes=None,
                  filterName=None,
                  **kwargs):
        """Perform the calculation specified by this plugin.

        This method can either be used to operate on a single catalog record
        or a whole catalog, populating it with the output defined by this
        plugin.

        Note that results may be added to catalog records as new columns, or
        may result in changes to existing values.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaFluxes : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        raise NotImplementedError()

    def fail(self, diaObject, columns, error=None):
        """Set diaObject position values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        columns : `list` of `str`
            List of string names of columns to write a the failed value.
        error : `BaseException` or `None`
            Error to pass. Kept for consistency with CatologCalculationPlugin.
            Unused.
        """
        for colName in columns:
            diaObject[colName] = np.nan


class DiaObjectCalculationConfig(CatalogCalculationConfig):
    """Config class for the catalog calculation driver task.

    Specifies which plugins will execute when the `CatalogCalculationTask`
    associated with this configuration is run.
    """

    plugins = DiaObjectCalculationPlugin.registry.makeField(
        multi=True,
        default=["ap_meanPosition",
                 "ap_meanFlux"],
        doc="Plugins to be run and their configuration")


class DiaObjectCalculationTask(CatalogCalculationTask):
    """Run plugins which operate on a catalog of DIA sources.

    This task facilitates running plugins which will operate on a source
    catalog. These plugins may do things such as classifying an object based
    on source record entries inserted during a measurement task.

    This task differs from CatalogCaculationTask in the following ways:

    -No multi mode is available for plugins. All plugins are assumed to run
     in single mode.

    -Input and output catalog types are assumed to be `pandas.DataFrames` with
     columns following those used in the Apdb.

    -No schema argument is passed to the plugins. Each plugin specifies
     output columns and required inputs.

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
    ConfigClass = DiaObjectCalculationConfig
    _DefaultName = "diaObjectCalculation"

    def __init__(self, plugMetadata=None, **kwargs):
        lsst.pipe.base.Task.__init__(self, **kwargs)
        if plugMetadata is None:
            plugMetadata = lsst.daf.base.PropertyList()
        self.plugMetadata = plugMetadata
        self.plugins = PluginMap()
        self.outputCols = []

        self.initializePlugins()

    def initializePlugins(self):
        """Initialize the plugins according to the configuration.
        """

        pluginType = namedtuple('pluginType', 'single multi')
        self.executionDict = {}
        # Read the properties for each plugin. Allocate a dictionary entry for
        # each run level. Verify that the plugins are above the minimum run
        # level for an catalogCalculation plugin. For each run level, the
        # plugins are sorted into either single record, or multi record groups
        # to later be run appropriately
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            if executionOrder not in self.executionDict:
                self.executionDict[executionOrder] = pluginType(single=[], multi=[])
            if PluginClass.getExecutionOrder() >= BasePlugin.DEFAULT_CATALOGCALCULATION:
                plug = PluginClass(config, name, metadata=self.plugMetadata)

                self._validatePluginCols(plug)

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

    def _validatePluginCols(self, plug):
        """Assert that output columns are not duplicated and input columns
        exist for dependent plugins.

        Parameters
        ----------
        plug : `lsst.ap.association.DiaCalculationPlugin`
            Plugin to test for output collisions and input needs.
        """
        for inputName in plug.inputCols:
            if inputName not in self.outputCols:
                errorTuple = (plug.name, plug.getExecutionOrder(),
                              inputName)
                raise ValueError(
                    "Plugin, {} with execution order {} requires DiaObject "
                    "column {} to exist. Check the execution order of the "
                    "plugin and make sure it runs after a plugin creating "
                    "the column is run.".format(*errorTuple))
        for outputName in plug.outputCols:
            if outputName in self.outputCols:
                errorTuple = (plug.name, plug.getExecutionOrder(),
                              outputName)
                raise ValueError(
                    "Plugin, {} with execution order {} is attempting to "
                    "output a column {}, however the column is already being "
                    "produced by another plugin. Check other plugins for "
                    "collisions with this one.".format(*errorTuple))
            else:
                self.outputCols.append(outputName)

    @lsst.pipe.base.timeMethod
    def run(self,
            diaObjectCat,
            diaSourceCat,
            updatedDiaObjectIds,
            filterNames):
        """The entry point for the DIA catalog calculation task.

        Run method both updates the values in the diaObjectCat and appends
        newly created DiaObjects to the catalog. For catalog column names
        see the lsst.cat schema definitions for the DiaObject and DiaSource
        tables (http://github.com/lsst/cat).

        Parameters
        ----------
        diaObjectCat : `pandas.DataFrame`
            DiaObjects to update values of and append new objects to. DataFrame
            should be indexed on "diaObjectId"
        diaSourceCat : `pandas.DataFrame`
            DiaSources associated with the DiaObjects in diaObjectCat.
            DataFrame should be indexed on
            `["diaObjectId", "filterName", "diaSourceId"]`
        updatedDiaObjectIds : `numpy.ndarray`
            Integer ids of the DiaObjects to update and create.
        filterNames : `list` of `str`
            List of string names of filters to be being processed.

        Returns
        -------
        returnStruct : `lsst.pipe.base.Struct`
            Struct containing:

            ``diaObjectCat``
                Full set of DiaObjects including both un-updated and
                updated/new DiaObjects (`pandas.DataFrame`).
            ``updatedDiaObjects``
                Catalog of DiaObjects  that were updated or created by this
                task (`pandas.DataFrame`).
        """
        if diaObjectCat.index.name is None:
            diaObjectCat.set_index("diaObjectId", inplace=True, drop=False)
        elif diaObjectCat.index.name != "diaObjectId":
            self.log.warning(
                "Input diaObjectCat is indexed on column(s) incompatible with "
                "this task. Should be indexed on 'diaObjectId'. Trying to set "
                "index regardless")
            diaObjectCat.set_index("diaObjectId", inplace=True, drop=False)

        # ``names`` by default is FrozenList([None]) hence we access the first
        # element and test for None.
        if diaSourceCat.index.names[0] is None:
            diaSourceCat.set_index(
                ["diaObjectId", "filterName", "diaSourceId"],
                inplace=True,
                drop=False)
        elif (diaSourceCat.index.names
              != ["diaObjectId", "filterName", "diaSourceId"]):
            self.log.warning(
                "Input diaSourceCat is indexed on column(s) incompatible with "
                "this task. Should be indexed on 'multi-index, "
                "['diaObjectId', 'filterName', 'diaSourceId']. Trying to set "
                "index regardless.")
            diaSourceCat.set_index(
                ["diaObjectId", "filterName", "diaSourceId"],
                inplace=True,
                drop=False)

        return self.callCompute(diaObjectCat,
                                diaSourceCat,
                                updatedDiaObjectIds,
                                filterNames)

    @lsst.pipe.base.timeMethod
    def callCompute(self,
                    diaObjectCat,
                    diaSourceCat,
                    updatedDiaObjectIds,
                    filterNames):
        """Run each of the plugins on the catalog.

        For catalog column names see the lsst.cat schema definitions for the
        DiaObject and DiaSource tables (http://github.com/lsst/cat).

        Parameters
        ----------
        diaObjectCat : `pandas.DataFrame`
            DiaObjects to update values of and append new objects to. DataFrame
            should be indexed on "diaObjectId"
        diaSourceCat : `pandas.DataFrame`
            DiaSources associated with the DiaObjects in diaObjectCat.
            DataFrame must be indexed on
            ["diaObjectId", "filterName", "diaSourceId"]`
        updatedDiaObjectIds : `numpy.ndarray`
            Integer ids of the DiaObjects to update and create.
        filterNames : `list` of `str`
            List of string names of filters to be being processed.

        Returns
        -------
        returnStruct : `lsst.pipe.base.Struct`
            Struct containing:

            ``diaObjectCat``
                Full set of DiaObjects including both un-updated and
                updated/new DiaObjects (`pandas.DataFrame`).
            ``updatedDiaObjects``
                Catalog of DiaObjects  that were updated or created by this
                task (`pandas.DataFrame`).

        Raises
        ------
        KeyError
            Raises if `pandas.DataFrame` indexing is not properly set.
        """
        # DiaObjects will be updated in place.
        diaObjectsToUpdate = diaObjectCat.loc[updatedDiaObjectIds, :]
        self.log.info("Calculating summary stats for %i DiaObjects" %
                      len(diaObjectsToUpdate))

        updatingDiaSources = diaSourceCat.loc[updatedDiaObjectIds, :]
        diaSourcesGB = updatingDiaSources.groupby(level=0)
        for runlevel in sorted(self.executionDict):
            for plug in self.executionDict[runlevel].single:
                if plug.needsFilter:
                    continue
                for updatedDiaObjectId in updatedDiaObjectIds:

                    # Sub-select diaSources associated with this diaObject.
                    objDiaSources = updatingDiaSources.loc[updatedDiaObjectId]

                    # Sub-select on diaSources observed in the current filter.
                    with CCContext(plug, updatedDiaObjectId, self.log):
                        # We feed the catalog we need to update and the id
                        # so as to get a few into the catalog and not a copy.
                        # This updates the values in the catalog.
                        plug.calculate(diaObjects=diaObjectsToUpdate,
                                       diaObjectId=updatedDiaObjectId,
                                       diaSources=objDiaSources,
                                       filterDiaSources=None,
                                       filterName=None)
            for plug in self.executionDict[runlevel].multi:
                if plug.needsFilter:
                    continue
                with CCContext(plug, diaObjectsToUpdate, self.log):
                    plug.calculate(diaObjects=diaObjectsToUpdate,
                                   diaSources=diaSourcesGB,
                                   filterDiaSources=None,
                                   filterName=None)

        for filterName in filterNames:
            try:
                updatingFilterDiaSources = updatingDiaSources.loc[
                    (slice(None), filterName), :
                ]
            except KeyError:
                self.log.warning(f"No DiaSource data with fitler={filterName}. "
                                 "Continuing...")
                continue
            # Level=0 here groups by diaObjectId.
            filterDiaSourcesGB = updatingFilterDiaSources.groupby(level=0)

            for runlevel in sorted(self.executionDict):
                for plug in self.executionDict[runlevel].single:
                    if not plug.needsFilter:
                        continue
                    for updatedDiaObjectId in updatedDiaObjectIds:

                        # Sub-select diaSources associated with this diaObject.
                        objDiaSources = updatingDiaSources.loc[updatedDiaObjectId]

                        # Sub-select on diaSources observed in the current filter.
                        try:
                            filterObjDiaSources = objDiaSources.loc[filterName]
                        except KeyError:
                            self.log.warning(
                                "DiaObjectId={updatedDiaObjectId} has no "
                                "DiaSources for filter={filterName}. "
                                "Continuing...")
                        with CCContext(plug, updatedDiaObjectId, self.log):
                            # We feed the catalog we need to update and the id
                            # so as to get a few into the catalog and not a copy.
                            # This updates the values in the catalog.
                            plug.calculate(diaObjects=diaObjectsToUpdate,
                                           diaObjectId=updatedDiaObjectId,
                                           diaSources=objDiaSources,
                                           filterDiaSources=filterObjDiaSources,
                                           filterName=filterName)
                for plug in self.executionDict[runlevel].multi:
                    if not plug.needsFilter:
                        continue
                    with CCContext(plug, diaObjectsToUpdate, self.log):
                        plug.calculate(diaObjects=diaObjectsToUpdate,
                                       diaSources=diaSourcesGB,
                                       filterDiaSources=filterDiaSourcesGB,
                                       filterName=filterName)
        # Need to store the newly updated diaObjects directly as the editing
        # a view into diaObjectsToUpdate does not update the values of
        # diaObjectCat.
        diaObjectCat.loc[updatedDiaObjectIds, :] = diaObjectsToUpdate
        return lsst.pipe.base.Struct(
            diaObjectCat=diaObjectCat,
            updatedDiaObjects=diaObjectsToUpdate)

    def _initialize_dia_object(self, objId):
        """Create a new DiaObject with values required to be initialized by the
        Apdb.

        Parameters
        ----------
        objid : `int`
            ``diaObjectId`` value for the of the new DiaObject.

        Returns
        -------
        diaObject : `dict`
            Newly created DiaObject with keys:

            ``diaObjectId``
                Unique DiaObjectId (`int`).
            ``pmParallaxNdata``
                Number of data points used for parallax calculation (`int`).
            ``nearbyObj1``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj2``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj3``
                Id of the a nearbyObject in the Object table (`int`).
            ``?PSFluxData``
                Number of data points used to calculate point source flux
                summary statistics in each bandpass (`int`).
        """
        new_dia_object = {"diaObjectId": objId,
                          "pmParallaxNdata": 0,
                          "nearbyObj1": 0,
                          "nearbyObj2": 0,
                          "nearbyObj3": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            new_dia_object["%sPSFluxNdata" % f] = 0
        return new_dia_object
