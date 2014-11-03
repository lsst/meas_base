#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python measureSources.py
or
   python
   >>> import measureSources; measureSources.run()
"""
import itertools
import math
import unittest
import numpy
import lsst.utils.tests as tests
import lsst.pex.logging as pexLogging
import lsst.pex.exceptions
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.meas.base as measBase
import lsst.meas.algorithms as measAlg

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9

FwhmPerSigma = 2*math.sqrt(2*math.log(2)) # FWHM for an N(0, 1) Gaussian

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MeasureSourcesTestCase(unittest.TestCase):
    """A test case for Measure"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testNaiveMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #
        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))
        alg = measBase.NaiveFluxAlgorithm
        control = alg.Control()
        control.radius = 10.0
        result = alg.Result()
        alg.apply(exp, afwGeom.Point2D(30 + x0, 50 + y0),result,control)
        flux = 3170.0
        self.assertEqual(result.flux, flux)

    def testCircularApertureMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #

        radii =  ( 1.0,   5.0,   10.0)  # radii to use

        control = measBase.CircularApertureFluxAlgorithm.Control()
        control.radii = radii

        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=float)
        schema.addField("centroid_y", type=float)

        alg = measBase.CircularApertureFluxAlgorithm(control, "myplugin", schema)
        table = afwTable.SourceTable.make(schema)
        table.defineCentroid("centroid")
        source = table.makeRecord()
        source.set("centroid_x", 30+x0)
        source.set("centroid_y", 50+y0)
        alg.measure(source, exp)
        for i, r in enumerate(radii):
            currentFlux = source.get("myplugin_flux_%d" % i)
            self.assertAlmostEqual(10.0*math.pi*r*r/currentFlux, 1.0, places=4)

    def testPeakLikelihoodFlux(self):
        """Test measurement with PeakLikelihoodFlux
        """
        # make mp: a flux measurer
        alg = measBase.PeakLikelihoodFluxAlgorithm

        # make and measure a series of exposures containing just one star, approximately centered
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(100, 101))
        kernelWidth = 35
        var = 100
        fwhm = 3.0
        sigma = fwhm/FwhmPerSigma
        convolutionControl = afwMath.ConvolutionControl()
        psf = measAlg.SingleGaussianPsf(kernelWidth, kernelWidth, sigma)
        psfKernel = psf.getLocalKernel()
        psfImage = psf.computeKernelImage()
        sumPsfSq = numpy.sum(psfImage.getArray()**2)
        psfSqArr = psfImage.getArray()**2

        for flux in (1000, 10000):
            ctrInd = afwGeom.Point2I(50, 51)
            ctrPos = afwGeom.Point2D(ctrInd)

            kernelBBox = psfImage.getBBox()
            kernelBBox.shift(afwGeom.Extent2I(ctrInd))

            # compute predicted flux error
            unshMImage = makeFakeImage(bbox, [ctrPos], [flux], fwhm, var)

            # filter image by PSF
            unshFiltMImage = afwImage.MaskedImageF(unshMImage.getBBox())
            afwMath.convolve(unshFiltMImage, unshMImage, psfKernel, convolutionControl)

            # compute predicted flux = value of image at peak / sum(PSF^2)
            # this is a sanity check of the algorithm, as much as anything
            predFlux = unshFiltMImage.getImage().get(ctrInd[0], ctrInd[1]) / sumPsfSq
            self.assertLess(abs(flux - predFlux), flux * 0.01)

            # compute predicted flux error based on filtered pixels
            # = sqrt(value of filtered variance at peak / sum(PSF^2)^2)
            predFluxErr = math.sqrt(unshFiltMImage.getVariance().get(ctrInd[0], ctrInd[1])) / sumPsfSq

            # compute predicted flux error based on unfiltered pixels
            # = sqrt(sum(unfiltered variance * PSF^2)) / sum(PSF^2)
            # and compare to that derived from filtered pixels;
            # again, this is a test of the algorithm
            varView = afwImage.ImageF(unshMImage.getVariance(), kernelBBox)
            varArr = varView.getArray()
            unfiltPredFluxErr = math.sqrt(numpy.sum(varArr*psfSqArr)) / sumPsfSq
            self.assertLess(abs(unfiltPredFluxErr - predFluxErr), predFluxErr * 0.01)

            for fracOffset in (afwGeom.Extent2D(0, 0), afwGeom.Extent2D(0.2, -0.3)):
                adjCenter = ctrPos + fracOffset
                if fracOffset == (0, 0):
                    maskedImage = unshMImage
                    filteredImage = unshFiltMImage
                else:
                    maskedImage = makeFakeImage(bbox, [adjCenter], [flux], fwhm, var)
                    # filter image by PSF
                    filteredImage = afwImage.MaskedImageF(maskedImage.getBBox())
                    afwMath.convolve(filteredImage, maskedImage, psfKernel, convolutionControl)

                exposure = afwImage.makeExposure(filteredImage)
                exposure.setPsf(psf)

                control = alg.Control()
                result = alg.Result()
                alg.apply(exposure, afwGeom.Point2D(*adjCenter), result, control)
                measFlux = result.flux
                measFluxErr = result.fluxSigma
                self.assertLess(abs(measFlux - flux), flux * 0.003)

                self.assertLess(abs(measFluxErr - predFluxErr), predFluxErr * 0.2)

                # try nearby points and verify that the flux is smaller;
                # this checks that the sub-pixel shift is performed in the correct direction
                for dx in (-0.2, 0, 0.2):
                    for dy in (-0.2, 0, 0.2):
                        if dx == dy == 0:
                            continue
                        offsetCtr = afwGeom.Point2D(adjCenter[0] + dx, adjCenter[1] + dy)
                        result = alg.Result()
                        alg.apply(exposure, offsetCtr, result)
                        self.assertLess(result.flux, measFlux)

        # source so near edge of image that PSF does not overlap exposure should result in failure

        for edgePos in (
            (1, 50),
            (50, 1),
            (50, bbox.getHeight() - 1),
            (bbox.getWidth() - 1, 50),
        ):
            self.assertRaises(
                lsst.pex.exceptions.RangeError,
                alg.apply,
                exposure,
                afwGeom.Point2D(*edgePos),
                alg.Result(),
            )

        # no PSF should result in failure: flags set
        noPsfExposure = afwImage.ExposureF(filteredImage)
        self.assertRaises(
            lsst.pex.exceptions.InvalidParameterError,
            alg.apply,
            noPsfExposure,
            afwGeom.Point2D(*adjCenter),
            alg.Result(),
        )

    def testPixelFlags(self):
        width, height = 100, 100
        mi = afwImage.MaskedImageF(width, height)
        exp = afwImage.makeExposure(mi)
        mi.getImage().set(0)
        mask = mi.getMask()
        sat = mask.getPlaneBitMask('SAT')
        interp = mask.getPlaneBitMask('INTRP')
        edge = mask.getPlaneBitMask('EDGE')
        bad = mask.getPlaneBitMask('BAD')
        mask.set(0)
        mask.set(20, 20, sat)
        mask.set(60, 60, interp)
        mask.set(40, 20, bad)
        mask.Factory(mask, afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(3, height))).set(edge)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))
        alg = measBase.PixelFlagsAlgorithm
        allFlags = [measBase.PixelFlagsAlgorithm.EDGE,
                    measBase.PixelFlagsAlgorithm.BAD,
                    measBase.PixelFlagsAlgorithm.SATURATED_CENTER,
                    measBase.PixelFlagsAlgorithm.SATURATED,
                    measBase.PixelFlagsAlgorithm.INTERPOLATED_CENTER,
                    measBase.PixelFlagsAlgorithm.INTERPOLATED,
                    ]
        for x, y, setFlags in [(1, 50, [measBase.PixelFlagsAlgorithm.EDGE]),
                               (40, 20, [measBase.PixelFlagsAlgorithm.BAD]),
                               (20, 20, [measBase.PixelFlagsAlgorithm.SATURATED_CENTER, measBase.PixelFlagsAlgorithm.SATURATED]),
                               (20, 22, [measBase.PixelFlagsAlgorithm.SATURATED]),
                               (60, 60, [measBase.PixelFlagsAlgorithm.INTERPOLATED_CENTER, measBase.PixelFlagsAlgorithm.INTERPOLATED]),
                               (60, 62, [measBase.PixelFlagsAlgorithm.INTERPOLATED]),
                    ]:
            foot = afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5)
            result = alg.Result()
            alg.apply(mi, afwGeom.Point2D(x + x0, y + y0), foot, result)
            for flag in allFlags:
                value = result.getFlag(flag)
                if flag in setFlags:
                    self.assertTrue(value, "Flag %s should be set for %f,%f" % (flag, x, y))
                else:
                    self.assertFalse(value, "Flag %s should not be set for %f,%f" % (flag, x, y))

        result = alg.Result()
        self.assertRaises(
            lsst.pex.exceptions.InvalidParameterError,
            alg.apply,
            mi,
            afwGeom.Point2D(float("NAN"), 40),
            afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5),
            alg.Result(),
            alg.Control()
        )

def addStar(image, center, flux, fwhm):
    """Add a perfect single Gaussian star to an image

    @warning uses Python to iterate over all pixels (because there is no C++
    function that computes a Gaussian offset by a non-integral amount).

    @param[in,out] image: Image to which to add star
    @param[in] center: position of center of star on image (pair of float)
    @param[in] flux: flux of Gaussian star, in counts
    @param[in] fwhm: FWHM of Gaussian star, in pixels
    """
    sigma = fwhm/FwhmPerSigma
    func = afwMath.GaussianFunction2D(sigma, sigma, 0)
    starImage = afwImage.ImageF(image.getBBox())
    # The flux in the region of the image will not be exactly the desired flux because the Gaussian
    # does not extend to infinity, so keep track of the actual flux and correct for it
    actFlux = 0
    # No function exists that has a fractional x and y offset, so set the image the slow way
    for i in range(image.getWidth()):
        x = center[0] - i
        for j in range(image.getHeight()):
            y = center[1] - j
            pixVal = flux * func(x, y)
            actFlux += pixVal
            starImage[i, j] += pixVal
    starImage *= flux / actFlux

    image += starImage

def makeFakeImage(bbox, centerList, fluxList, fwhm, var):
    """Make a fake image containing a set of stars variance = image + var

    (It is trivial to add Poisson noise, which would be more accurate,
    but hard to make a unit test  that can reliably determine whether such an image passes a test)

    @param[in] bbox: bounding box for image
    @param[in] centerList: list of positions of center of star on image (pairs of float)
    @param[in] fluxList: flux of each star, in counts
    @param[in] fwhm: FWHM of Gaussian star, in pixels
    @param[in] var: value of variance plane (counts)
    """
    if len(centerList) != len(fluxList):
        raise RuntimeError("len(centerList) != len(fluxList)")
    maskedImage = afwImage.MaskedImageF(bbox)
    image = maskedImage.getImage()
    for center, flux in itertools.izip(centerList, fluxList):
        addStar(image, center=center, flux=flux, fwhm=fwhm)
    variance = maskedImage.getVariance()
    variance[:] = image
    variance += var
    return maskedImage


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureSourcesTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
