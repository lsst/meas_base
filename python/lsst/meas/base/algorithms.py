from lsst.meas.base.sfm import *
from lsst.meas.base.forced import *
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

# ------------------------------------------------

class ForcedFluxConfig(ForcedAlgorithmConfig):
    pass

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class ForcedFlux(ForcedAlgorithm):
    ConfigClass = ForcedFluxConfig
    doMeasureSingle = True
    doMeasureMulti = False

    def __init__(self, config, name, schemaMapper, flags=None, others=None, metadata=None):
        schemaMapper.addOutputField(lsst.afw.table.Field_D("test.flux", "sum of flux in object footprint", ""))
        schemaMapper.addOutputField(lsst.afw.table.Field_I("test.fluxcount", "number of pixels in object footprint", ""))
        schemaMapper.addOutputField(lsst.afw.table.Field_D("test.back", "avg of flux in background", ""))
        schemaMapper.addOutputField(lsst.afw.table.Field_I("test.backcount", "number of pixels in surrounding background", ""))
        self.config = config

    def measureSingle(self, exposure, source, ref, refWcs):

        schema = source.getSchema()
        fluxkey = schema.find("test.flux").key
        fluxcountkey = schema.find("test.fluxcount").key
        backkey = schema.find("test.back").key
        backcountkey = schema.find("test.backcount").key
        id = source.getId()

        expregion = lsst.afw.geom.Box2I(exposure.getXY0(),
            lsst.afw.geom.Extent2I(exposure.getWidth(), exposure.getHeight()))
        targetWcs = exposure.getWcs()
        foot = source.getFootprint().transform(refWcs, targetWcs, expregion, True)
 
        image = exposure.getMaskedImage().getImage()
        array = image.getArray()

        # sum the footprint area for this source
        sumarray = numpy.ndarray((foot.getArea()), array.dtype)
        lsst.afw.detection.flattenArray(foot, image.getArray(), sumarray, image.getXY0())
        flux = sumarray.sum(dtype=numpy.float64)
        area = foot.getArea()
        source.set(fluxkey, flux)
        source.set(fluxcountkey, area)

        # Now find an area which is 100 pixels larger in all directions than the foot.getBBox()
        fbbox = foot.getBBox()
        border = 100
        xmin = fbbox.getMinX() - border
        ymin = fbbox.getMinY() - border
        xmax = fbbox.getMaxX() + border
        ymax = fbbox.getMaxY() + border
        x0 = image.getX0()
        y0 = image.getY0()
        if xmin < x0: xmin = x0
        if ymin < y0: ymin = y0
        if xmax > (x0 + exposure.getWidth()): xmax = x0+exposure.getWidth()
        if ymax > (y0 + exposure.getHeight()): ymax = y0+exposure.getHeight()
        bigarraysub = array[ymin-y0:ymax-y0, xmin-x0:xmax-x0]
        bigflux = bigarraysub.sum(dtype=numpy.float64)
        bigarea = (ymax-ymin)*(xmax-xmin)
        source.set(backkey, bigflux - flux)
        source.set(backcountkey, bigarea - area)


    def measureMulti(self, exposure, sources):
        return


ForcedAlgorithm.registry.register("forced.flux", ForcedFlux)

