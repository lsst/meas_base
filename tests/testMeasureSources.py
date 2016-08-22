#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""
Tests for measuring sources using the meas_base framework

"""
from __future__ import absolute_import, division, print_function
import itertools
import math
import unittest

import numpy as np

import lsst.pex.logging as pexLogging
import lsst.pex.exceptions
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.utils.tests

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("afwDetection.Measure", verbose)

try:
    type(display)
except NameError:
    display = False

FwhmPerSigma = 2*math.sqrt(2*math.log(2))  # FWHM for an N(0, 1) Gaussian


def makePluginAndCat(alg, name, control, metadata=False, centroid=None):
    schema = afwTable.SourceTable.makeMinimalSchema()
    if centroid:
        schema.addField(centroid + "_x", type=float)
        schema.addField(centroid + "_y", type=float)
        schema.addField(centroid + "_flag", type='Flag')
        schema.getAliasMap().set("slot_Centroid", centroid)
    if metadata:
        plugin = alg(control, name, schema, dafBase.PropertySet())
    else:
        plugin = alg(control, name, schema)
    cat = afwTable.SourceCatalog(schema)
    return plugin, cat

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class MeasureSourcesTestCase(lsst.utils.tests.TestCase):
    """A test case for Measure"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testCircularApertureMeasure(self):
        mi = afwImage.MaskedImageF(afwGeom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #

        radii = (1.0,   5.0,   10.0)  # radii to use

        control = measBase.ApertureFluxControl()
        control.radii = radii

        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))

        plugin, cat = makePluginAndCat(measBase.CircularApertureFluxAlgorithm,
                                       "test", control, True, centroid="centroid")
        source = cat.makeRecord()
        source.set("centroid_x", 30+x0)
        source.set("centroid_y", 50+y0)
        plugin.measure(source, exp)

        for r in radii:
            currentFlux = source.get("%s_flux" %
                                     measBase.CircularApertureFluxAlgorithm.makeFieldPrefix("test", r))
            self.assertAlmostEqual(10.0*math.pi*r*r/currentFlux, 1.0, places=4)

    def testPeakLikelihoodFlux(self):
        """Test measurement with PeakLikelihoodFlux."""
        # make and measure a series of exposures containing just one star, approximately centered
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(100, 101))
        kernelWidth = 35
        var = 100
        fwhm = 3.0
        sigma = fwhm/FwhmPerSigma
        convolutionControl = afwMath.ConvolutionControl()
        psf = afwDetection.GaussianPsf(kernelWidth, kernelWidth, sigma)
        psfKernel = psf.getLocalKernel()
        psfImage = psf.computeKernelImage()
        sumPsfSq = np.sum(psfImage.getArray()**2)
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
            unfiltPredFluxErr = math.sqrt(np.sum(varArr*psfSqArr)) / sumPsfSq
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

                exp = afwImage.makeExposure(filteredImage)
                exp.setPsf(psf)
                control = measBase.PeakLikelihoodFluxControl()
                plugin, cat = makePluginAndCat(measBase.PeakLikelihoodFluxAlgorithm, "test",
                                               control, centroid="centroid")
                source = cat.makeRecord()
                source.set("centroid_x", adjCenter.getX())
                source.set("centroid_y", adjCenter.getY())
                plugin.measure(source, exp)
                measFlux = source.get("test_flux")
                measFluxErr = source.get("test_fluxSigma")
                self.assertLess(abs(measFlux - flux), flux * 0.003)

                self.assertLess(abs(measFluxErr - predFluxErr), predFluxErr * 0.2)

                # try nearby points and verify that the flux is smaller;
                # this checks that the sub-pixel shift is performed in the correct direction
                for dx in (-0.2, 0, 0.2):
                    for dy in (-0.2, 0, 0.2):
                        if dx == dy == 0:
                            continue
                        offsetCtr = afwGeom.Point2D(adjCenter[0] + dx, adjCenter[1] + dy)
                        source = cat.makeRecord()
                        source.set("centroid_x", offsetCtr.getX())
                        source.set("centroid_y", offsetCtr.getY())
                        plugin.measure(source, exp)
                        self.assertLess(source.get("test_flux"), measFlux)

        # source so near edge of image that PSF does not overlap exposure should result in failure
        for edgePos in (
            (1, 50),
            (50, 1),
            (50, bbox.getHeight() - 1),
            (bbox.getWidth() - 1, 50),
        ):
            source = cat.makeRecord()
            source.set("centroid_x", edgePos[0])
            source.set("centroid_y", edgePos[1])
            with self.assertRaises(lsst.pex.exceptions.RangeError):
                plugin.measure(source,exp)

        # no PSF should result in failure: flags set
        noPsfExposure = afwImage.ExposureF(filteredImage)
        source = cat.makeRecord()
        source.set("centroid_x", edgePos[0])
        source.set("centroid_y", edgePos[1])
        with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
            plugin.measure(source, noPsfExposure)

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
        nodata = mask.getPlaneBitMask('NO_DATA')
        mask.addMaskPlane('CLIPPED')
        clipped = mask.getPlaneBitMask('CLIPPED')
        mask.set(0)
        mask.set(20, 20, sat)
        mask.set(60, 60, interp)
        mask.set(40, 20, bad)
        mask.set(20, 80, nodata)
        mask.set(30, 30, clipped)
        mask.Factory(mask, afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(3, height))).set(edge)
        x0, y0 = 1234, 5678
        exp.setXY0(afwGeom.Point2I(x0, y0))
        control = measBase.PixelFlagsControl()
        # Change the configuration of control to test for clipped mask
        control.masksFpAnywhere = ['CLIPPED']
        plugin, cat = makePluginAndCat(measBase.PixelFlagsAlgorithm, "test", control, centroid="centroid")
        allFlags = [
            "",
            "edge",
            "interpolated",
            "interpolatedCenter",
            "saturated",
            "saturatedCenter",
            "cr",
            "crCenter",
            "bad",
            "clipped",
        ]
        for x, y, setFlags in [(1, 50, ['edge']),
                               (40, 20, ['bad']),
                               (20, 20, ['saturatedCenter',
                                         'saturated']),
                               (20, 22, ['saturated']),
                               (60, 60, ['interpolatedCenter',
                                         'interpolated']),
                               (60, 62, ['interpolated']),
                               (20, 80, ['edge']),
                               (30, 30, ['clipped']),
                               ]:
            foot = afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5)
            source = cat.makeRecord()
            source.setFootprint(foot)
            source.set("centroid_x", x+x0)
            source.set("centroid_y", y+y0)
            plugin.measure(source, exp)
            for flag in allFlags[1:]:
                value = source.get("test_flag_" + flag)
                if flag in setFlags:
                    self.assertTrue(value, "Flag %s should be set for %f,%f" % (flag, x, y))
                else:
                    self.assertFalse(value, "Flag %s should not be set for %f,%f" % (flag, x, y))

        # the new code which grabs the center of a record throws when a Nan is set in the
        # centroid slot and the algorithm attempts to get the default center position
        source = cat.makeRecord()
        source.set("centroid_x", float("NAN"))
        source.set("centroid_y", 40)
        source.set("centroid_flag", True)
        source.setFootprint(afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5))
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source,exp)
        # Test that if there is no center and centroider that the object should look at the footprint
        plugin, cat = makePluginAndCat(measBase.PixelFlagsAlgorithm, "test", control)
        # The first test should raise exception because there is no footprint
        source = cat.makeRecord()
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source, exp)
        # The second test will raise an error because no peaks are present
        source.setFootprint(afwDetection.Footprint(afwGeom.Point2I(afwGeom.Point2D(x + x0, y + y0)), 5))
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source, exp)
        # The final test should pass because it detects a peak, we are reusing the location of the
        # clipped bit in the mask plane, so we will check first that it is false, then true
        source.getFootprint().addPeak(x+x0, y+y0, 100)
        self.assertFalse(source.get("test_flag_clipped"), "The clipped flag should be set False")
        plugin.measure(source, exp)
        self.assertTrue(source.get("test_flag_clipped"), "The clipped flag should be set True")


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

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
