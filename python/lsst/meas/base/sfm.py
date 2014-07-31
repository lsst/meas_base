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

__all__ = ("SingleFramePluginConfig", "SingleFramePlugin", "WrappedSingleFramePlugin",
           "SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")

class SingleFramePluginConfig(BasePluginConfig):
    """Base class for configs of single-frame plugin algorithms."""
    pass

class SingleFramePlugin(BasePlugin):
    """Base class for single-frame plugin algorithms."""

    # All subclasses of SingleFramePlugin should be registered here
    registry = PluginRegistry(SingleFramePluginConfig)
    ConfigClass = SingleFramePluginConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        @param[in,out]  schema   The Source schema.  New fields should be added here to
                                 hold measurements produced by this plugin.
        @param[in]  flags        A set of bitflags describing the data that the plugin
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       A PluginMap of previously-initialized plugins
        @param[in]  metadata     Plugin metadata that will be attached to the output catalog
        """
        self.config = config
        self.name = name

    def measure(self, measRecord, exposure):
        """Measure the properties of a source on a single image
        (single-epoch image or coadd).

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
        """Measure the properties of a group of blended sources on a single image
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
    """A base class for SingleFramePlugins that delegates the algorithmic work to a C++
    Algorithm class.

    Derived classes of WrappedSingleFramePlugin must set the AlgClass class attribute
    to the C++ class being wrapped, which must meet the requirements defined in the
    "Implementing New Plugins and Algorithms" section of the meas_base documentation.  This is usually done
    by calling the generate() class method.
    """

    AlgClass = None

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
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
        """Create a new derived class of WrappedSingleFramePlugin from a C++ Algorithm class.

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
    """Config class for single frame measurement driver task."""

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
                 #"correctfluxes",
                 "base_ClassificationExtendedness",
                 "base_SkyCoord",
                 ],
        doc="Plugin plugins to be run and their configuration"
        )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")

## \addtogroup LSST_task_documentation
## \{
## \page singleFrameMeasurementTask
## \ref SingleFrameMeasurementTask_ "SingleFrameMeasurementTask"
## \copybrief SingleFrameMeasurementTask
## \}

class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """!
\anchor SingleFrameMeasurementTask_
\brief The SingleFrameMeasurementTask is used to measure the properties of sources on a single exposure.


\section meas_base_sfm_Contents Contents

 - \ref meas_base_sfm_Purpose
 - \ref meas_base_sfm_Initialize
 - \ref meas_base_sfm_IO
 - \ref meas_base_sfm_Config
 - \ref meas_base_sfm_Example

\section meas_base_sfm_Purpose	Description

\copybrief SingleFrameMeasurementTask

The task is configured with a list of "plugins": each plugin defines the values it
is measuring and conducts that measurement on each detected source. The job of the
measurement task is to call each plugin at the appropriate time for initialization
and measurement of each detected source, and to save the results of the
measurement to a Source Record.

\section meas_base_sfm_Initialize	Task initialisation

\copydoc init

\section meas_base_sfm_IO		Inputs/Outputs to the run method

\copydoc run 

\section meas_base_sfm_Config       Configuration parameters

See \ref SingleFrameMeasurementConfig

\section meas_base_sfm_Example	A complete example of using SingleFrameMeasurementTask

This code is in \link runSFMTask.py\endlink in the examples directory, and can be run as \em e.g.
\code
examples/runSFMTask.py --ds9
\endcode
\dontinclude runSFMTask.py

See \ref meas_algorithms_detection_Example for a few more details on the DetectionTask.

Import the tasks (there are some other standard imports; read the file if you're confused)
\skip SourceDetectionTask
\until SingleFrameMeasurementTask

We need to create our tasks before processing any data as the task constructors
can add extra columns to the schema.  First the detection task
\skipline makeMinimalSchema
\skip SourceDetectionTask.ConfigClass
\until detectionTask
and then the measurement task using the default plugins (as set by SingleFrameMeasurementConfig.plugins):
\skipline SingleFrameMeasurementTask.ConfigClass
\until measureTask

We're now ready to process the data (we could loop over multiple exposures/catalogues using the same
task objects).  First create the output table and process the image to find sources:
\skipline afwTable
\skip result
\until sources

Then measure them:
\skipline measure

We then might plot the results (\em e.g. if you set \c --ds9 on the command line)
\skip display
\until RED

\dontinclude runSFMTask.py
Rather than accept a default set you can select which plugins should be run.
First create the Config object:
\skipline SingleFrameMeasurementTask.ConfigClass
Then specify which plugins we're interested in and set any needed parameters:
\until radii = radii

Unfortunately that won't quite work as there are still "slots" (mappings between measurements like PSF fluxes
and the plugins that calculate them) pointing to some of the discarded plugins (see SourceSlotConfig):

\skip instFlux
\until psfFlux
and create the task as before:
\skipline measureTask
and create the task as before:
\skipline measureTask
We can find out what aperture radii were chosen with
\skipline radii
and add them to the display code:
\skip s in sources
\until YELLOW

and end up with something like
\image html measAlgTasks-ds9.png

    """
class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """Single-frame measurement driver task"""

    ConfigClass = SingleFrameMeasurementConfig
    _DefaultName = "measurement"
    tableVersion = 1

    #   The algMetadata parameter is currently required by the pipe_tasks running mechanism
    #   This is a temporary state until pipe_tasks is converted to the new plugin framework.
    def init(self, schema, algMetadata=None, flags=None, **kwds):
        """!Initialize the task. Set up the execution order of the plugins and initialize
        the plugins, giving each plugin an opportunity to add its measurement fields to
        the output schema and to record information in the task metadata.

        \param[in,out] schema      lsst.afw.table.Schema, to be initialized to include the
                                   measurement fields from the plugins already
        \param[in,out] algMetaData lsst.daf.base.PropertyList used to record information about
                                   each algorithm.
        \param[in]     flags       lsst.meas.base.MeasurementDataFlags
        \param[in]     **kwds      Keyword arguments passed from lsst.pipe.base.task.Task
        """
        __init__(self, schema, algMetadata=None, flags=None, **kwds)

    #   The algMetadata parameter is currently required by the pipe_tasks running mechanism
    #   This is a temporary state until pipe_tasks is converted to the new plugin framework.
    def __init__(self, schema, algMetadata=None, flags=None, **kwds):
        """!\copydoc init
        """
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.schema = schema
        self.algMetadata = lsst.daf.base.PropertyList()
        self.plugins = PluginMap()
        # Init the plugins, sorted by execution order.  At the same time add to the schema
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            self.plugins[name] = PluginClass(config, name, schema=schema, flags=flags,
                others=self.plugins, metadata=self.algMetadata)

    def run(self, measCat, exposure):
        """!Run single frame measurement over an exposure and source catalog

        \param[in,out]  measCat  lsst.afw.table.SourceCatalog to be filled with outputs.  Must
                                 contain all the SourceRecords to be measured (with Footprints
                                 attached), and have a schema that is a superset of self.schema.

        \param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.
        """
        if exposure.__class__.__name__ == "SourceCatalog":
            temp = exposure
            exposure = measCat
            measCat = temp
        assert measCat.getSchema().contains(self.schema)
        self.config.slots.setupTable(measCat.table)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
            for measRecord in measCat}

        # noiseReplacer is used to fill the footprints with noise and save off heavy footprints
        # of the source pixels so that they can re restored one at a time for measurement.
        # After the NoiseReplacer is constructed, all pixels in the exposure.getMaskedImage()
        # which belong to objects in measCat will be replaced with noise
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)

        # First, create a catalog of all parentless sources
        # Loop through all the parent sources, first processing the children, then the parent
        measParentCat = measCat.getChildren(0)

        self.log.info("Measuring %d sources (%d parents, %d children) "
                      % (len(measCat), len(measParentCat), len(measCat) - len(measParentCat)))

        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                callMeasure(self, measChildRecord, exposure)
                noiseReplacer.removeSource(measChildRecord.getId())
            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            callMeasure(self, measParentRecord, exposure)
            # Finally, process both the parent and the child set through measureN
            callMeasureN(self, measParentCat[parentIdx:parentIdx+1], exposure)
            callMeasureN(self, measChildCat, exposure)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()

    def measure(self, measCat, exposure):
        self.run(measCat, exposure)

