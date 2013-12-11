from lsst.meas.base.sfm import *
import lsst.meas.base.base
import lsst.afw.detection
import numpy
import lsst.meas.algorithms

__all__ = ('TestNull', 'SdssCentroid')

class TestNullConfig(SingleFrameAlgorithmConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=-1.0, doc="sets relative order of algorithms")

#  Algoritm which does nothing This is a placeholder for now
class TestNull(SingleFrameAlgorithm):

    ConfigClass = TestNullConfig

    def measureSingle(self, exposure, source):
        return

    def measureMulti(self, exposure, sources):
        return

class TestFluxConfig(SingleFrameAlgorithmConfig):
    pass

SingleFrameAlgorithm.registry.register("test.null", TestNull)

class SdssCentroidConfig(SingleFrameAlgorithmConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=-1.0, doc="sets relative order of algorithms")

class SdssCentroid(SingleFrameAlgorithmCpp):

    ConfigClass = SdssCentroidConfig

    def measureSingle(self, exposure, source):
        try: 
            foot = source.getFootprint() 
            peak = foot.getPeaks()[0]
            center = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
            self.alg.apply(source, exposure, center)
        except Exception as e:
            pass #print "err in sdss.centroid: %s, id = %d"%(e.message,source.getId())

SingleFrameAlgorithm.registry.register("centroid.sdss", SdssCentroid)

class SdssShapeConfig(SingleFrameAlgorithmConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=1.0, doc="sets relative order of algorithms")

class SdssShape(SingleFrameAlgorithmCpp):

    ConfigClass = SdssShapeConfig

SingleFrameAlgorithm.registry.register("shape.sdss", SdssShape)

class FluxGaussianConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxGaussian(SingleFrameAlgorithmCpp):

    ConfigClass = FluxGaussianConfig


SingleFrameAlgorithm.registry.register("flux.gaussian", FluxGaussian)


class FlagsPixelConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class FlagsPixel(SingleFrameAlgorithmCpp):

    ConfigClass = FlagsPixelConfig


SingleFrameAlgorithm.registry.register("flags.pixel", FlagsPixel)

class CentroidGaussianConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class CentroidGaussian(SingleFrameAlgorithmCpp):

    ConfigClass = CentroidGaussianConfig


SingleFrameAlgorithm.registry.register("centroid.gaussian", CentroidGaussian)

class CentroidNaiveConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class CentroidNaive(SingleFrameAlgorithmCpp):

    ConfigClass = CentroidNaiveConfig


SingleFrameAlgorithm.registry.register("centroid.naive", CentroidNaive)


class FluxNaiveConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxNaive(SingleFrameAlgorithmCpp):

    ConfigClass = FluxNaiveConfig


SingleFrameAlgorithm.registry.register("flux.naive", FluxNaive)

class FluxPsfConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxPsf(SingleFrameAlgorithmCpp):

    ConfigClass = FluxPsfConfig


SingleFrameAlgorithm.registry.register("flux.psf", FluxPsf)

class FluxSincConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxSinc(SingleFrameAlgorithmCpp):

    ConfigClass = FluxSincConfig


SingleFrameAlgorithm.registry.register("flux.sinc", FluxSinc)

class ClassificationExtendednessConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class ClassificationExtendedness(SingleFrameAlgorithmCpp):

    ConfigClass = ClassificationExtendednessConfig


SingleFrameAlgorithm.registry.register("classification.extendedness", ClassificationExtendedness)

class SkyCoordConfig(SingleFrameAlgorithmConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class SkyCoord(SingleFrameAlgorithmCpp):

    ConfigClass = SkyCoordConfig

SingleFrameAlgorithm.registry.register("skycoord", SkyCoord)
