"""Base classes for single-frame measurement plugin algorithms and the driver task for these.

In single-frame measurement, we assumes that detection and probably deblending have already been run on
the same frame, so a SourceCatalog has already been created with Footprints (which may be HeavyFootprints).
Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).
"""

import pdb
import lsst.pex.config
import lsst.pipe.base
import lsst.daf.base
from .base import *
__all__ = ("SingleFrameAlgorithmConfig", "SingleFrameAlgorithm", "SingleFrameMeasurementConfig", "SingleFrameMeasurementTask")
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
        default=[
              "test.null", #TODO
            ],
        doc="Plugin algorithms to be run and their configuration"
        )
    doReplaceWithNoise = lsst.pex.config.Field(dtype=bool, default=True, optional=False,
                                         doc='When measuring, replace other detected footprints with noise?')

    #  This is to allow the internal NoiseReplacer to be parameterized.
    noiseSource = lsst.pex.config.ChoiceField(doc='How do we choose the mean and variance of the Gaussian noise we generate?',
                                      dtype=str, allowed={
                                          'measure': 'Measure clipped mean and variance from the whole image',
                                          'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
                                          'variance': "Mean = 0, variance = the image's variance",
                                          },
                                      default='measure',
                                      optional=False)

    noiseOffset = lsst.pex.config.Field(dtype=float, optional=False, default=0.,
                                  doc='Add ann offset to the generated noise.')

    noiseSeed = lsst.pex.config.Field(dtype=int, default=0, doc='The seed value to use for random number generation.')



class SingleFrameMeasurementTask(lsst.pipe.base.Task):
    """Single-frame measurement driver task"""

    ConfigClass = SingleFrameMeasurementConfig
    _DefaultName = "measurement"

    # FIX:  The algMetadata parameter is currently required by the pipe_base running mechanism
    #       Problem should be resolved when the plugins are converted.
    def __init__(self, schema, algMetadata=None, **kwds):
        flags = None
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.schema = schema
        self.algMetadata = lsst.daf.base.PropertyList()
        self.algorithms = AlgorithmMap()
        for executionOrder, name, config, AlgorithmClass in sorted(self.config.algorithms.apply()):
            self.algorithms[name] = AlgorithmClass(config, name, schema=schema, flags=flags,
                                                   others=self.algorithms, metadata=self.algMetadata)
    def getChildren(self, cat, id):
        retcat = lsst.afw.table.SourceCatalog(cat.getSchema())
        for rec in cat:
            if rec.getParent() == id:
                retcat.append(rec)
        #view = cat.getColumnView()
        #mask = view["parent"] == id
        #return cat.subset(mask)
        return retcat

    def run(self, exposure, measCat):
        assert measCat.getSchema().contains(self.schema)
        #measCat.writeFits("measCat.fits")
        exposure.getMaskedImage().getImage().writeFits("image.fits")
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in measCat}
        # noiseReplacer is used to fill the footprints with noise and save off heavy footprints
        # of what was in the exposure beforehand
        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource, self.config.noiseOffset, self.config.noiseSeed)

        # loop through all the parent sources, processing the children then the parent
        measParentCat = measCat.getChildren(0)
        self.log.info("There are %d parent sources"%len(measParentCat))
        for parentIdx, measParentRecord in enumerate(measParentCat):
            measChildCat = measCat.getChildren(measParentRecord.getId())
            for measChildRecord in measChildCat:
                noiseReplacer.insertSource(measChildRecord.getId())
                for algorithm in self.algorithms.iterSingle():
                    algorithm.measureSingle(exposure, measChildRecord)
                noiseReplacer.removeSource(measChildRecord.getId())
            noiseReplacer.insertSource(measParentRecord.getId())
            for algorithm in self.algorithms.iterSingle():
                algorithm.measureSingle(exposure, measParentRecord)
            for algorithm in self.algorithms.iterMulti():
                algorithm.measureMulti(exposure, measParentCat[parentIndex:parentIndex+1])
                algorithm.measureMulti(exposure, measChildCat)
            noiseReplacer.removeSource(measParentRecord.getId())
        noiseReplacer.end()


