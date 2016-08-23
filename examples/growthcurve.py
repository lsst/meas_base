#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 Aura/LSST
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

#
# Original filename: examples/growthcurve.py
#
# Author: Steve Bickerton
# Email: bick@astro.princeton.edu
# Date: Mon 2009-10-26 13:42:37
#
# Summary:
#
# python version growthcurve.cc example
#
"""
%prog [options] arg
"""
from __future__ import print_function
import sys
import optparse
import math

from builtins import map
from builtins import range
from builtins import object
import numpy

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.meas.base as measBase

# =====================================================================
# a functor for the PSF
#


class Gaussian(object):  # public std::binary_function<double, double, double> {

    def __init__(self, xcen, ycen, sigma, a):
        self.xcen = xcen
        self.ycen = ycen
        self.sigma = sigma
        self.a = a

    def __call__(self, x, y):
        xx = x - self.xcen
        yy = y - self.ycen
        ss = self.sigma*self.sigma
        coeff = self.a * (1.0/(2.0*numpy.pi*ss))
        expon = numpy.exp(-(xx*xx + yy*yy) / (2.0*ss))
        return coeff*expon


# =====================================================================
# a radial functor for the PSF
#
# This functor isn't currently used in the routine
# I'll leave it here in case I (someday) figure out how to integrate a python functor
class RGaussian(object):  # public std::unary_function<double, double> {

    def __init__(self, sigma, a, apradius, aptaper):
        self.sigma = sigma
        self.a = a
        self.apradius = apradius
        self.aptaper = aptaper

    def __call__(self, r):
        ss = self.sigma*self.sigma
        gauss = self.a * (1.0/(2.0*numpy.pi*ss)) * numpy.exp(-(r*r)/(2.0*ss))
        aperture = 0.0
        if (r <= apradius):
            aperture = 1.0
        elif (r > apradius and r < apradius + aptaper):
            aperture = 0.5*(1.0 + cos(numpy.pi*(r - apradius)/aptaper))
        return aperture*gauss*(r*2.0*numpy.pi)


#############################################################
#
# Main body of code
#
#############################################################

def main():

    ########################################################################
    # command line arguments and options
    ########################################################################

    parser = optparse.OptionParser(usage=__doc__)
    # parser.add_option("-a", "--aa", dest="aa", type=float,
    #                  default=1.0, help="default=%default")

    opts, args = parser.parse_args()

    if len(args) == 0:
        r1, r2, dr = 3.0, 3.0, 0.5
    elif len(args) == 3:
        r1, r2, dr = list(map(float, args))
    else:
        parser.print_help()
        sys.exit(1)

    # make a list of radii to compute the growthcurve points
    radius = []
    nR = int((r2 - r1)/dr + 1)
    for iR in range(nR):
        radius.append(r1 + iR*dr)

    # make an image big enough to hold the largest requested aperture
    xwidth = 2*(0 + 128)
    ywidth = xwidth

    # initializations
    sigmas = [1.5, 2.5]  # the Gaussian widths of the psfs we'll use
    nS = len(sigmas)
    alg = measBase.PsfFluxAlgorithm
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("centroid_x", type=float)
    schema.addField("centroid_y", type=float)
    plugin = alg(measBase.PsfFluxControl(), "test", schema)
    cat = afwTable.SourceCatalog(schema)
    cat.defineCentroid("centroid")
    print("# sig rad  Naive Sinc Psf")
    for iS in range(nS):
        sigma = sigmas[iS]

        # make a Gaussian star to measure
        gauss = afwMath.GaussianFunction2D(sigma, sigma)
        kernel = afwMath.AnalyticKernel(xwidth, ywidth, gauss)
        kimg = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimg, False)
        kimg *= 100.0
        mimg = afwImage.MaskedImageF(kimg.convertFloat(),
                                     afwImage.MaskU(kimg.getDimensions(), 0x0),
                                     afwImage.ImageF(kimg.getDimensions(), 0.0))
        exposure = afwImage.ExposureF(mimg)

        # loop over all the radii in the growthcurve
        for iR in range(nR):

            psfH = int(2.0*(r2 + 2.0)) + 1
            psfW = int(2.0*(r2 + 2.0)) + 1

            # this test uses a single Gaussian instead of the original double gaussian
            psf = afwDet.GaussianPsf(psfW, psfH, sigma)

            # get the aperture fluxes for Naive and Sinc methods

            axes = afwGeom.ellipses.Axes(radius[iR], radius[iR], math.radians(0))
            center = afwGeom.Point2D(0, 0)
            ellipse = afwGeom.ellipses.Ellipse(axes, center)
            resultSinc = measBase.ApertureFluxAlgorithm.computeSincFlux(mimg.getImage(), ellipse)
            resultNaive = measBase.ApertureFluxAlgorithm.computeNaiveFlux(mimg.getImage(), ellipse)
            source = cat.addNew()
            source["centroid_x"] = 0
            source["centroid_y"] = 0
            exposure.setPsf(psf)
            plugin.measure(source, exposure)
            fluxNaive = resultNaive.flux
            fluxSinc = resultSinc.flux
            fluxPsf = source["test_flux"]

            # get the exact flux for the theoretical smooth PSF
            # rpsf = RGaussian(sigma, a, radius[iR], aptaper)
            # *** not sure how to integrate a python functor ***

            print("%.2f %.2f  %.3f %.3f %.3f" % (sigma, radius[iR], fluxNaive, fluxSinc, fluxPsf))


#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()
