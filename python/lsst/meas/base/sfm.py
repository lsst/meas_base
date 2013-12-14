"""Base classes for single-frame measurement plugin algorithms and the driver task for these.

In single-frame measurement, we assumes that detection and probably deblending have already been run on
the same frame, so a SourceCatalog has already been created with Footprints (which may be HeavyFootprints).
Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).
"""

import lsst.pex.config
import lsst.pipe.base
import lsst.daf.base
from .base import *
__all__ = ("SingleFrameAlgorithmConfig", "SingleFrameAlgorithm", "SingleFrameMeasurementConfig",
           "SingleFrameMeasurementTask")
class SingleFrameAlgorithmConfig(BaseAlgorithmConfig):
    """Base class for configs of single-frame plugin algorithms."""
    pass

class SingleFrameAlgorithm(BaseAlgorithm):
    """Base class for single-frame plugin algorithms."""


    # All subclasses of SingleFrameAlgorithm should be registered here
    registry = AlgorithmRegistry(SingleFrameAlgorithmConfig)
    ConfigClass = SingleFrameAlgorithmConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the algorithm was registered with.
        @param[in,out]  schema   The Source schema.  New fields should be added here to
                                 hold measurements produced by this algorithm.
        @param[in]  flags        A set of bitflags describing the data that the algorithm
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       An AlgorithmMap of previously-initialized algorithms
        @param[in]  metadata     Algorithm metadata that will be attached to the output catalog
        """
        self.config = config

    def measureSingle(self, exposure, source):
        """Measure the properties of a source on a single image
        (single-epoch image or coadd).

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  All
                                 other sources in the image will have been replaced by
                                 noise according to deblender outputs.

        @param[in,out] measRecord  lsst.afw.table.SourceRecord to be filled with outputs,
                                   and from which previously-measured quantities can be
                                   retreived.

        """
        raise NotImplementedError()

    def measureMulti(self, exposure, sources):
        """Measure the properties of a group of blended sources on a single image
        (single-epoch image or coadd).

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  Sources
                                 not in the blended hierarchy to be measured will have
                                 been replaced with noise using deblender outputs.

        @param[in,out] measCat   lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be
                                 retrieved, containing only the sources that should be
                                 measured together in this call.

        """
        raise NotImplementedError()

class SingleFrameMeasurementConfig(lsst.pex.config.Config):
    """Config class for single-frame measurement driver task."""

    algorithms = SingleFrameAlgorithm.registry.makeField(
        multi=True,
        default=[], #TODO
        doc="Plugin algorithms to be run and their configuration"
        )

    doReplaceWithNoise = lsst.pex.config.Field(dtype=bool, default=True, optional=False,
        doc='When measuring, replace other detected footprints with noise?')

    #  These parameters are to allow the internal NoiseReplacer to be parameterized.
    noiseSource = lsst.pex.config.ChoiceField(
        doc='How to choose mean and variance of the Gaussian noise we generate?',
        dtype=str, allowed={'measure': 'Measure clipped mean and variance from the whole image',
        'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
        'variance': "Mean = 0, variance = the image's variance",},
        default='measure', optional=False)
    noiseOffset = lsst.pex.config.Field(dtype=float, optional=False, default=0.,
                                  doc='Add ann offset to the generated noise.')
    noiseSeed = lsst.pex.config.Field(dtype=int, default=0,
        doc='The seed value to use for random number generation.')

class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """Single-frame measurement driver task"""

    ConfigClass = SingleFrameMeasurementConfig
    _DefaultName = "measurement"

    """ Initialize the task, including setting up the execution order of the algorithms
        and providing the task with the metadata and schema objects

        @param[in] schema      lsst.afw.table.Schema, which should have been initialized
                               to include the measurement fields from the algorithms already
    """
    #   The algMetadata parameter is currently required by the pipe_base running mechanism
    #   This is a termporary state untile the plugins are converted.
    def __init__(self, schema, algMetadata=None, **kwds):
        flags = None
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.schema = schema
        self.algMetadata = lsst.daf.base.PropertyList()
        self.algorithms = AlgorithmMap()
        for executionOrder, name, config, AlgorithmClass in sorted(self.config.algorithms.apply()):
            self.algorithms[name] = AlgorithmClass(config, name, schema=schema, flags=flags,
                others=self.algorithms, metadata=self.algMetadata)

    """ Run single frame measurement over an exposure and source catalog

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.

        @param[in, out] measCat  lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be retreived.

    """
    def run(self, exposure, measCat):
        assert measCat.getSchema().contains(self.schema)
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
            for measRecord in measCat}
        # noiseReplacer is used to fill the footprints with noise and save off heavy footprints
        # of what was in the exposure beforehand.
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)

        # loop through all the parent sources, processing the children, then the parent
        measParentCat = measCat.getChildren(0)
        for parentIdx, measParentRecord in enumerate(measParentCat):
            # first get all the children of this parent, insert footprint in turn, and measure
            measChildCat = measCat.getChildren(measParentRecord.getId())
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                for algorithm in self.algorithms.iterSingle():
                    algorithm.measureSingle(exposure, measChildRecord)
                noiseReplacer.removeSource(measChildRecord.getId())
            # Then insert the parent footprint, and measure that
            noiseReplacer.insertSource(measParentRecord.getId())
            for algorithm in self.algorithms.iterSingle():
                algorithm.measureSingle(exposure, measParentRecord)
            for algorithm in self.algorithms.iterMulti():
                algorithm.measureMulti(exposure, measParentCat[parentIndex:parentIndex+1])
                algorithm.measureMulti(exposure, measChildCat)
            noiseReplacer.removeSource(measParentRecord.getId())
        # when done, restore the exposure to its original state
        noiseReplacer.end()


