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
"""Base classes for single-frame measurement plugins and the driver task for these.

In single-frame measurement, we assume that detection and probably deblending have already been run on
the same frame, so a SourceCatalog has already been created with Footprints (which may be HeavyFootprints).
Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).
"""

import lsst.pex.config
import lsst.pipe.base
import lsst.daf.base

from .base import *
from .noiseReplacer import *

__all__ = ("SingleFramePluginConfig", "SingleFramePlugin", "WrappedSingleFramePlugin",
           "SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")

class SingleFramePluginConfig(BasePluginConfig):
    """!
    Base class for configs of single-frame plugin algorithms.
    """
    pass

class SingleFramePlugin(BasePlugin):
    """!
    Base class for single-frame plugin algorithms.

    New Plugins can be created in Python by inheriting directly from this class
    and implementing measure(), fail() (from BasePlugin), and optionally __init__
    and measureN().  Plugins can also be defined in C++ via the WrappedSingleFramePlugin
    class.
    """

    # All subclasses of SingleFramePlugin should be registered here
    registry = PluginRegistry(SingleFramePluginConfig)
    ConfigClass = SingleFramePluginConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        """!
        Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        @param[in,out]  schema   The Source schema.  New fields should be added here to
                                 hold measurements produced by this plugin.
        @param[in]  flags        A set of bitflags describing the data that the plugin
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       A PluginMap of previously-initialized plugins
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        super(SingleFramePlugin, self).__init__()
        self.config = config
        self.name = name

    def measure(self, measRecord, exposure):
        """!
        Measure the properties of a source on a single image (single-epoch image or coadd).

        @param[in,out] measRecord  lsst.afw.table.SourceRecord to be filled with outputs,
                                   and from which previously-measured quantities can be
                                   retreived.

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  All
                                 other sources in the image will have been replaced by
                                 noise according to deblender outputs.

        """
        raise NotImplementedError()

    def measureN(self, measCat, exposure):
        """!
        Measure the properties of a group of blended sources on a single image
        (single-epoch image or coadd).

        @param[in,out] measCat   lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be
                                 retrieved, containing only the sources that should be
                                 measured together in this call.

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  Sources
                                 not in the blended hierarchy to be measured will have
                                 been replaced with noise using deblender outputs.

        Derived classes that do not implement measureN() should just inherit this
        disabled version.  Derived classes that do implement measureN() should additionally
        add a bool doMeasureN config field to their config class to signal that measureN-mode
        is available.
        """
        raise NotImplementedError()

class WrappedSingleFramePlugin(SingleFramePlugin):
    """!
    A base class for SingleFramePlugins that delegates the algorithmic work to a C++
    Algorithm class.

    Derived classes of WrappedSingleFramePlugin must set the AlgClass class attribute
    to the C++ class being wrapped, which must meet the requirements defined in the
    "Implementing New Plugins and Algorithms" section of the meas_base documentation.  This is usually done
    by calling the generate() class method.
    """

    AlgClass = None

    def __init__(self, config, name, schema, flags, others, metadata):
        super(WrappedSingleFramePlugin, self).__init__(config, name, schema, flags, others, metadata)
        self.resultMapper = self.AlgClass.makeResultMapper(schema, name, config.makeControl())
        # TODO: check flags

    def measure(self, measRecord, exposure):
        inputs = self.AlgClass.Input(measRecord)
        result = self.AlgClass.apply(exposure, inputs, self.config.makeControl())
        self.resultMapper.apply(measRecord, result)

    def measureN(self, measCat, exposure):
        assert hasattr(AlgClass, "applyN")  # would be better if we could delete this method somehow
        inputs = self.AlgClass.Input.Vector(measCat)
        results = self.AlgClass.applyN(exposure, inputs, self.config.makeControl())
        for result, measRecord in zip(results, measCat):
            self.resultMapper.apply(measRecord, result)

    def fail(self, measRecord, error=None):
        # The ResultMapper will set detailed flag bits describing the error if error is not None,
        # and set a general failure bit otherwise.
        if self.name == measRecord.getTable().getCentroidDefinition():
            if len(measRecord.getFootprint().getPeaks()) > 0:
                measRecord[self.name+"_x"] = measRecord.getFootprint().getPeaks()[0].getCentroid().getX()
                measRecord[self.name+"_y"] = measRecord.getFootprint().getPeaks()[0].getCentroid().getY()
                measRecord[self.name+"_xSigma"] = 0
                measRecord[self.name+"_ySigma"] = 0
        self.resultMapper.fail(measRecord, error)

    @classmethod
    def generate(Base, AlgClass, name=None, doRegister=True, ConfigClass=None, executionOrder=None):
        """!
        Create a new derived class of WrappedSingleFramePlugin from a C++ Algorithm class.

        @param[in]   AlgClass   The name of the (Swigged) C++ Algorithm class this Plugin will delegate to.
        @param[in]   name       The name to use when registering the Plugin (ignored if doRegister=False).
                                Defaults to the result of generateAlgorithmName(AlgClass).
        @param[in]   doRegister   If True (default), register the new Plugin so it can be configured to be
                                  run by SingleFrameMeasurementTask.
        @param[in]   ConfigClass  The ConfigClass associated with the new Plugin.  This should have a
                                  makeControl() method that returns the Control object used by the C++
                                  Algorithm class.
        @param[in]   executionOrder   If not None, a float that overrides the default executionOrder
                                      for this algorithm (see BasePluginConfig.executionOrder).

        For more information, please see the "Adding New Algorithms" section of the main meas_base
        documentation.
        """
        if ConfigClass is None:
            ConfigClass = lsst.pex.config.makeConfigClass(AlgClass.Control, base=Base.ConfigClass,
                                                          module=AlgClass.__module__)
            if hasattr(AlgClass, "applyN"):
                ConfigClass.doMeasureN = lsst.pex.config.Field(
                    dtype=bool, default=True,
                    doc="whether to run this plugin multi-object mode"
                    )
        if executionOrder is not None:
            ConfigClass.executionOrder.default = float(executionOrder)
        PluginClass = type(AlgClass.__name__ + "SingleFramePlugin", (Base,),
                           dict(AlgClass=AlgClass, ConfigClass=ConfigClass))
        if doRegister:
            if name is None:
                name = generateAlgorithmName(AlgClass)
            Base.registry.register(name, PluginClass)
        return PluginClass

class SingleFrameMeasurementConfig(BaseMeasurementConfig):
    """!
    Config class for single frame measurement driver task.
    """

    plugins = SingleFramePlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_SdssCentroid",
                 "base_GaussianCentroid",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_NaiveFlux",
                 "base_PsfFlux",
                 "base_SincFlux",
                 "base_ClassificationExtendedness",
                 "base_SkyCoord",
                 ],
        doc="Plugins to be run and their configuration"
        )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")

## @addtogroup LSST_task_documentation
## @{
## @page singleFrameMeasurementTask
## SingleFrameMeasurementTask
## @copybrief SingleFrameMeasurementTask
## @}

class SingleFrameMeasurementTask(BaseMeasurementTask):
    """!
    A subtask for measuring the properties of sources on a single exposure.

    The task is configured with a list of "plugins": each plugin defines the values it
    measures (i.e. the columns in a table it will fill) and conducts that measurement
    on each detected source (see SingleFramePlugin).  The job of the
    measurement task is to initialize the set of plugins (which includes setting up the
    catalog schema) from their configuration, and then invoke each plugin on each
    source.

    When run after the deblender (see lsst.meas.deblender.SourceDeblendTask),
    SingleFrameMeasurementTask also replaces each source's neighbors with noise before
    measuring each source, utilizing the HeavyFootprints created by the deblender (see
    NoiseReplacer).

    SingleFrameMeasurementTask has only two methods: __init__() and run().  For configuration
    options, see SingleFrameMeasurementConfig.

    @section meas_base_sfm_Example	A complete example of using SingleFrameMeasurementTask

    The code below is in examples/runSingleFrameTask.py

    @dontinclude runSingleFrameTask.py

    See meas_algorithms_detection_Example for more information on SourceDetectionTask.

    Import the tasks (there are some other standard imports; read the file if you're confused)
    @skip SourceDetectionTask
    @until SingleFrameMeasurementTask

    We need to create our tasks before processing any data as the task constructors
    can add extra columns to the schema.  The most important argument we pass these to these
    is a lsst.afw.table.Schema object, which contains information about the fields (i.e. columns) of the
    measurement catalog we'll create, including names, types, and additional documentation.
    Tasks that operate on a catalog are typically passed a Schema upon construction, to which
    they'll add the fields they'll fill later when run.  We construct a mostly empty Schema that
    contains just the fields required for a SourceCatalog like this:
    @skipline schema
    Now we can configure and create the SourceDetectionTask:
    @until detectionTask
    (tableVersion is a temporary parameter that will be removed after the transition from the meas_algorithms
    Tasks to the meas_base Tasks is complete; for now, setting tableVersion=1 is necessary when using
    meas_base Tasks via pipe_tasks drivers).

    We then move on to configuring the measurement task:
    @until config
    While a reasonable set of plugins is configured by default, we'll customize the list.
    We also need to unset one of the slots at the same time, because we're
    not running the algorithm that it' set to by default, and that would cause problems later:
    @until config.slot
    MeasurementDataFlags are just a placeholder for now; eventually they'll be used to tell plugins about
    certain features of the data, such as whether the data being processed is a difference image, a regular
    image, or a coadd, or whether it has been "preconvolved" by the PSF.
    @until flags

    Now, finally, we can construct the measurement task:
    @until measureTask

    After constructing all the tasks, we can inspect the Schema we've created:
    @skipline print schema
    All of the fields in the
    schema can be accessed via the get() method on a record object.  See afwTable for more
    information.

    We're now ready to process the data (we could loop over multiple exposures/catalogs using the same
    task objects).  First create the output table and process the image to find sources:
    @skipline afwTable
    @skip result
    @until sources

    Then measure them:
    @skipline measure

    We then might plot the results (@em e.g. if you set @c --ds9 on the command line)
    @skip display
    @until RED
    and end up with something like
    @image html runSingleFrameTask-ds9.png
    """

    ConfigClass = SingleFrameMeasurementConfig

    def __init__(self, schema, algMetadata=None, flags=None, **kwds):
        """!
        Initialize the task. Set up the execution order of the plugins and initialize
        the plugins, giving each plugin an opportunity to add its measurement fields to
        the output schema and to record information in the task metadata.

        @param[in,out] schema      lsst.afw.table.Schema, to be initialized to include the
                                   measurement fields from the plugins already
        @param[in,out] algMetadata lsst.daf.base.PropertyList used to record information about
                                   each algorithm.  An empty PropertyList will be created if None.
        @param[in]     flags       lsst.meas.base.MeasurementDataFlags
        @param[in]     **kwds      Keyword arguments passed from lsst.pipe.base.task.Task
        """
        super(SingleFrameMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.schema = schema
        # Init the plugins, sorted by execution order.  At the same time add to the schema
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, schema=schema, flags=flags,
                others=self.plugins, metadata=self.algMetadata)

    def run(self, measCat, exposure):
        """!
        Run single frame measurement over an exposure and source catalog

        @param[in,out]  measCat  lsst.afw.table.SourceCatalog to be filled with outputs.  Must
                                 contain all the SourceRecords to be measured (with Footprints
                                 attached), and have a schema that is a superset of self.schema.

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.
        """
        # Temporary workaround for change in order of arguments; will be removed when transition
        # from meas_algorithms to meas_base is complete.
        if exposure.__class__.__name__ == "SourceCatalog":
            temp = exposure
            exposure = measCat
            measCat = temp
        assert measCat.getSchema().contains(self.schema)
        self.config.slots.setupTable(measCat.table)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
            for measRecord in measCat}

        # noiseReplacer is used to fill the footprints with noise and save heavy footprints
        # of the source pixels so that they can be restored one at a time for measurement.
        # After the NoiseReplacer is constructed, all pixels in the exposure.getMaskedImage()
        # which belong to objects in measCat will be replaced with noise
        noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure, footprints, log=self.log)

        # First, create a catalog of all parentless sources
        # Loop through all the parent sources, first processing the children, then the parent
        measParentCat = measCat.getChildren(0)

        self.log.info("Measuring %d sources (%d parents, %d children) "
                      % (len(measCat), len(measParentCat), len(measCat) - len(measParentCat)))

        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                self.callMeasure(measChildRecord, exposure)
                noiseReplacer.removeSource(measChildRecord.getId())
            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            self.callMeasure(measParentRecord, exposure)
            # Finally, process both the parent and the child set through measureN
            self.callMeasureN(measParentCat[parentIdx:parentIdx+1], exposure)
            self.callMeasureN(measChildCat, exposure)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()

    def measure(self, measCat, exposure):
        """!
        Backwards-compatibility alias for run()
        """
        self.run(measCat, exposure)
