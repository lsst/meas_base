#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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

import numpy

import lsst.afw.table
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.geom.ellipses
import lsst.afw.coord
import lsst.pex.exceptions

from .sfm import SingleFrameMeasurementTask
from .forcedMeasurement import ForcedMeasurementTask
from .baseLib import CentroidResultKey

__all__ = ("BlendContext", "TestDataset", "AlgorithmTestCase", "TransformTestCase",
           "SingleFramePluginTransformSetupHelper", "ForcedPluginTransformSetupHelper",
           "FluxTransformTestCase", "CentroidTransformTestCase")


class BlendContext(object):
    """!
    A Python context manager used to add multiple overlapping sources along with a parent source
    that represents all of them together.

    This is used as the return value for TestDataset.addBlend(), and this is the only way it should
    be used.  The only public method is addChild().
    """

    def __init__(self, owner):
        self.owner = owner
        self.parentRecord = self.owner.catalog.addNew()
        self.parentImage = lsst.afw.image.ImageF(self.owner.exposure.getBBox())
        self.children = []

    def __enter__(self):
        # BlendContext is its own context manager, so we just return self.
        return self

    def addChild(self, flux, centroid, shape=None):
        """!
        Add a child source to the blend, and return the truth catalog record that corresponds to it.

        @param[in]  flux      Total flux of the source to be added.
        @param[in]  centroid  Position of the source to be added (lsst.afw.geom.Point2D).
        @param[in]  shape     2nd moments of the source before PSF convolution
                              (lsst.afw.geom.ellipses.Quadrupole).  Note that the truth catalog
                              records post-convolution moments)
        """
        record, image = self.owner.addSource(flux, centroid, shape)
        record.set(self.owner.keys["parent"], self.parentRecord.getId())
        self.parentImage += image
        self.children.append((record, image))
        return record

    def __exit__(self, type_, value, tb):
        # We're not using the context manager for any kind of exception safety or guarantees;
        # we just want the nice "with" statement syntax.
        if type_ is not None:  # exception was raised; just skip all this and let it propagate
            return
        # On exit, we need to compute and set the truth values for the parent object.
        self.parentRecord.set(self.owner.keys["nChild"], len(self.children))
        # Compute flux from sum of component fluxes
        flux = 0.0
        for record, image in self.children:
            flux += record.get(self.owner.keys["flux"])
        self.parentRecord.set(self.owner.keys["flux"], flux)
        # Compute centroid from flux-weighted mean of component centroids
        x = 0.0
        y = 0.0
        for record, image in self.children:
            w = record.get(self.owner.keys["flux"])/flux
            x += record.get(self.owner.keys["centroid"].getX())*w
            y += record.get(self.owner.keys["centroid"].getY())*w
        self.parentRecord.set(self.owner.keys["centroid"], lsst.afw.geom.Point2D(x, y))
        # Compute shape from flux-weighted mean of offset component shapes
        xx = 0.0
        yy = 0.0
        xy = 0.0
        for record, image in self.children:
            w = record.get(self.owner.keys["flux"])/flux
            dx = record.get(self.owner.keys["centroid"].getX()) - x
            dy = record.get(self.owner.keys["centroid"].getY()) - y
            xx += (record.get(self.owner.keys["shape"].getIxx()) + dx**2)*w
            yy += (record.get(self.owner.keys["shape"].getIyy()) + dy**2)*w
            xy += (record.get(self.owner.keys["shape"].getIxy()) + dx*dy)*w
        self.parentRecord.set(self.owner.keys["shape"], lsst.afw.geom.ellipses.Quadrupole(xx, yy, xy))
        # Run detection on the parent image to get the parent Footprint.
        self.owner._installFootprint(self.parentRecord, self.parentImage)
        # Create perfect HeavyFootprints for all children; these will need to be modified later to account
        # for the noise we'll add to the image.
        deblend = lsst.afw.image.MaskedImageF(self.owner.exposure.getMaskedImage(), True)
        for record, image in self.children:
            deblend.getImage().getArray()[:, :] = image.getArray()
            heavyFootprint = lsst.afw.detection.HeavyFootprintF(self.parentRecord.getFootprint(), deblend)
            record.setFootprint(heavyFootprint)


class TestDataset(object):
    """!
    A simulated dataset consisting of a test image and an associated truth catalog.

    TestDataset creates an idealized image made of pure Gaussians (including a Gaussian PSF),
    with simple noise and idealized Footprints/HeavyFootprints that simulated the outputs
    of detection and deblending.  Multiple noise realizations can be created from the same
    underlying sources, allowing uncertainty estimates to be verified via Monte Carlo.

    Typical usage:
    @code
    bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0,0), lsst.afw.geom.Point2I(100, 100))
    dataset = TestDataset(bbox)
    dataset.addSource(flux=1E5, centroid=lsst.afw.geom.Point2D(25, 26))
    dataset.addSource(flux=2E5, centroid=lsst.afw.geom.Point2D(75, 24),
                      shape=lsst.afw.geom.ellipses.Quadrupole(8, 7, 2))
    with dataset.addBlend() as family:
        family.addChild(flux=2E5, centroid=lsst.afw.geom.Point2D(50, 72))
        family.addChild(flux=1.5E5, centroid=lsst.afw.geom.Point2D(51, 74))
    exposure, catalog = dataset.realize(noise=100.0, schema=TestDataset.makeMinimalSchema())
    @endcode
    """

    @classmethod
    def makeMinimalSchema(cls):
        """Return the minimal schema needed to hold truth catalog fields.

        When TestDataset.realize() is called, the schema must include at least these fields.
        Usually it will include additional fields for measurement algorithm outputs, allowing
        the same catalog to be used for both truth values (the fields from the minimal schema)
        and the measurements.
        """
        if not hasattr(cls, "_schema"):
            schema = lsst.afw.table.SourceTable.makeMinimalSchema()
            cls.keys = {}
            cls.keys["parent"] = schema.find("parent").key
            cls.keys["nChild"] = schema.addField("deblend_nChild", type=int)
            cls.keys["flux"] = schema.addField("truth_flux", type=float, doc="true flux", units="count")
            cls.keys["centroid"] = lsst.afw.table.Point2DKey.addFields(
                schema, "truth", "true simulated centroid", "pixel"
            )
            cls.keys["centroid_flag"] = schema.addField("truth_flag", type="Flag",
                                                        doc="set if the object is a star")
            cls.keys["shape"] = lsst.afw.table.QuadrupoleKey.addFields(
                schema, "truth", "true shape after PSF convolution", lsst.afw.table.CoordinateType_PIXEL
            )
            cls.keys["isStar"] = schema.addField("truth_isStar", type="Flag",
                                                 doc="set if the object is a star")
            schema.getAliasMap().set("slot_Shape", "truth")
            schema.getAliasMap().set("slot_Centroid", "truth")
            schema.getAliasMap().set("slot_ModelFlux", "truth")
            schema.getCitizen().markPersistent()
            cls._schema = schema
        schema = lsst.afw.table.Schema(cls._schema)
        schema.disconnectAliases()
        return schema

    @staticmethod
    def makePerturbedWcs(oldWcs, minScaleFactor=1.2, maxScaleFactor=1.5,
                         minRotation=None, maxRotation=None,
                         minRefShift=None, maxRefShift=None,
                         minPixShift=2.0, maxPixShift=4.0):
        """!
        Create a new undistorted TanWcs that is similar but not identical to another, with random
        scaling, rotation, and offset (in both pixel position and reference position).

        The maximum and minimum arguments are interpreted as absolute values for a split
        range that covers both positive and negative values (as this method is used
        in testing, it is typically most important to avoid perturbations near zero).
        Scale factors are treated somewhat differently: the actual scale factor is chosen between
        minScaleFactor and maxScaleFactor OR (1/maxScaleFactor) and (1/minScaleFactor).

        The default range for rotation is 30-60 degrees, and the default range for reference shift
        is 0.5-1.0 arcseconds (these cannot be safely included directly as default values because Angle
        objects are mutable).
        """
        if minRotation is None:
            minRotation = 30.0*lsst.afw.geom.degrees
        if maxRotation is None:
            maxRotation = 60.0*lsst.afw.geom.degrees
        if minRefShift is None:
            minRefShift = 0.5*lsst.afw.geom.arcseconds
        if maxRefShift is None:
            maxRefShift = 1.0*lsst.afw.geom.arcseconds

        def splitRandom(min1, max1, min2=None, max2=None):
            if min2 is None:
                min2 = -max1
            if max2 is None:
                max2 = -min1
            if numpy.random.uniform() > 0.5:
                return float(numpy.random.uniform(min1, max1))
            else:
                return float(numpy.random.uniform(min2, max2))
        # Generate random perturbations
        scaleFactor = splitRandom(minScaleFactor, maxScaleFactor, 1.0/maxScaleFactor, 1.0/minScaleFactor)
        rotation = splitRandom(minRotation.asRadians(), maxRotation.asRadians())*lsst.afw.geom.radians
        refShiftRa = splitRandom(minRefShift.asRadians(), maxRefShift.asRadians())*lsst.afw.geom.radians
        refShiftDec = splitRandom(minRefShift.asRadians(), maxRefShift.asRadians())*lsst.afw.geom.radians
        pixShiftX = splitRandom(minPixShift, maxPixShift)
        pixShiftY = splitRandom(minPixShift, maxPixShift)
        # Compute new CD matrix
        oldTransform = lsst.afw.geom.LinearTransform(oldWcs.getCDMatrix())
        rTransform = lsst.afw.geom.LinearTransform.makeRotation(rotation)
        sTransform = lsst.afw.geom.LinearTransform.makeScaling(scaleFactor)
        newTransform = oldTransform*rTransform*sTransform
        matrix = newTransform.getMatrix()
        # Compute new coordinate reference pixel (CRVAL)
        oldSkyOrigin = oldWcs.getSkyOrigin().toIcrs()
        newSkyOrigin = lsst.afw.coord.IcrsCoord(oldSkyOrigin.getRa() + refShiftRa,
                                                oldSkyOrigin.getDec() + refShiftDec)
        # Compute new pixel reference pixel (CRPIX)
        oldPixOrigin = oldWcs.getPixelOrigin()
        newPixOrigin = lsst.afw.geom.Point2D(oldPixOrigin.getX() + pixShiftX,
                                             oldPixOrigin.getY() + pixShiftY)
        return lsst.afw.image.makeWcs(newSkyOrigin, newPixOrigin,
                                      matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1])

    @staticmethod
    def makeEmptyExposure(bbox, wcs=None, crval=None, cdelt=None, psfSigma=2.0, psfDim=17, fluxMag0=1E12):
        """!
        Create an Exposure, with a Calib, Wcs, and Psf, but no pixel values set.

        @param[in]    bbox        Bounding box of the image (image coordinates) as returned by makeCatalog.
        @param[in]    wcs         New Wcs for the exposure (created from crval and cdelt if None).
        @param[in]    crval       afw.coord.Coord: center of the TAN WCS attached to the image.
        @param[in]    cdelt       afw.geom.Angle: pixel scale of the image
        @param[in]    psfSigma    Radius (sigma) of the Gaussian PSF attached to the image
        @param[in]    psfDim      Width and height of the image's Gaussian PSF attached to the image
        @param[in]    fluxMag0    Flux at magnitude zero (in e-) used to set the Calib of the exposure.
        """
        if wcs is None:
            if crval is None:
                crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
            if cdelt is None:
                cdelt = 0.2*lsst.afw.geom.arcseconds
            crpix = lsst.afw.geom.Box2D(bbox).getCenter()
            wcs = lsst.afw.image.makeWcs(crval, crpix, cdelt.asDegrees(), 0.0, 0.0, cdelt.asDegrees())
        exposure = lsst.afw.image.ExposureF(bbox)
        psf = lsst.afw.detection.GaussianPsf(psfDim, psfDim, psfSigma)
        calib = lsst.afw.image.Calib()
        calib.setFluxMag0(fluxMag0)
        exposure.setWcs(wcs)
        exposure.setPsf(psf)
        exposure.setCalib(calib)
        return exposure

    @staticmethod
    def drawGaussian(bbox, flux, ellipse):
        """!
        Create an image of an elliptical Gaussian.

        @param[in,out] bbox        Bounding box of image to create.
        @param[in]     flux        Total flux of the Gaussian (normalized analytically, not using pixel
                                   values)
        @param[in]     ellipse     lsst.afw.geom.ellipses.Ellipse holding the centroid and shape.
        """
        x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX()),
                              numpy.arange(bbox.getBeginY(), bbox.getEndY()))
        t = ellipse.getGridTransform()
        xt = t[t.XX] * x + t[t.XY] * y + t[t.X]
        yt = t[t.YX] * x + t[t.YY] * y + t[t.Y]
        image = lsst.afw.image.ImageF(bbox)
        image.getArray()[:, :] = numpy.exp(-0.5*(xt**2 + yt**2))*flux/(2.0*ellipse.getCore().getArea())
        return image

    def __init__(self, bbox, threshold=10.0, exposure=None, **kwds):
        """!
        Initialize the dataset.

        @param[in]   bbox       Bounding box of the test image.
        @param[in]   threshold  Threshold absolute value used to determine footprints for
                                simulated sources.  This thresholding will be applied before noise is
                                actually added to images (or before the noise level is even known), so
                                this will necessarily produce somewhat artificial footprints.
        @param[in]   exposure   lsst.afw.image.ExposureF test sources should be added to.  Ownership should
                                be considered transferred from the caller to the TestDataset.
                                Must have a Gaussian Psf for truth catalog shapes to be exact.
        @param[in]   **kwds     Keyword arguments forwarded to makeEmptyExposure if exposure is None.
        """
        if exposure is None:
            exposure = self.makeEmptyExposure(bbox, **kwds)
        self.threshold = lsst.afw.detection.Threshold(threshold, lsst.afw.detection.Threshold.VALUE)
        self.exposure = exposure
        self.psfShape = self.exposure.getPsf().computeShape()
        self.schema = self.makeMinimalSchema()
        self.catalog = lsst.afw.table.SourceCatalog(self.schema)

    def _installFootprint(self, record, image):
        """Create a Footprint for a simulated source and add it to its truth catalog record.
        """
        # Run detection on the single-source image
        fpSet = lsst.afw.detection.FootprintSet(image, self.threshold)
        # the call below to the FootprintSet ctor is actually a grow operation
        fpSet = lsst.afw.detection.FootprintSet(fpSet, int(self.psfShape.getDeterminantRadius() + 1.0), True)
        # Update the full exposure's mask plane to indicate the detection
        fpSet.setMask(self.exposure.getMaskedImage().getMask(), "DETECTED")
        # Attach the new footprint to the exposure
        if len(fpSet.getFootprints()) > 1:
            raise RuntimeError("Threshold value results in multiple Footprints for a single object")
        if len(fpSet.getFootprints()) == 0:
            raise RuntimeError("Threshold value results in zero Footprints for object")
        record.setFootprint(fpSet.getFootprints()[0])

    def addSource(self, flux, centroid, shape=None):
        """!
        Add a source to the simulation

        @param[in]  flux      Total flux of the source to be added.
        @param[in]  centroid  Position of the source to be added (lsst.afw.geom.Point2D).
        @param[in]  shape     2nd moments of the source before PSF convolution
                              (lsst.afw.geom.ellipses.Quadrupole).  Note that the truth catalog
                              records post-convolution moments).  If None, a point source
                              will be added.

        @return a truth catalog record and single-source image corresponding to the new source.
        """
        # Create and set the truth catalog fields
        record = self.catalog.addNew()
        record.set(self.keys["flux"], flux)
        record.set(self.keys["centroid"], centroid)
        if shape is None:
            record.set(self.keys["isStar"], True)
            fullShape = self.psfShape
        else:
            record.set(self.keys["isStar"], False)
            fullShape = shape.convolve(self.psfShape)
        record.set(self.keys["shape"], fullShape)
        # Create an image containing just this source
        image = self.drawGaussian(self.exposure.getBBox(), flux,
                                  lsst.afw.geom.ellipses.Ellipse(fullShape, centroid))
        # Generate a footprint for this source
        self._installFootprint(record, image)
        # Actually add the source to the full exposure
        self.exposure.getMaskedImage().getImage().getArray()[:, :] += image.getArray()
        return record, image

    def addBlend(self):
        """!
        Return a context manager that allows a blend of multiple sources to be added.

        Example:
        @code
        d = TestDataset(...)
        with d.addBlend() as b:
            b.addChild(flux1, centroid1)
            b.addChild(flux2, centroid2, shape2)
        @endcode

        Note that nothing stops you from creating overlapping sources just using the addSource() method,
        but addBlend() is necesssary to create a parent object and deblended HeavyFootprints of the type
        produced by the detection and deblending pipelines.
        """
        return BlendContext(self)

    def transform(self, wcs, **kwds):
        """!
        Create a copy of the dataset transformed to a new WCS, with new Psf and Calib.

        @param[in]  wcs      Wcs for the new dataset.
        @param[in]  **kwds   Additional keyword arguments passed on to makeEmptyExposure.  If not
                             specified, these revert to the defaults for makeEmptyExposure, not the
                             values in the current dataset.
        """
        bboxD = lsst.afw.geom.Box2D()
        xyt = lsst.afw.image.XYTransformFromWcsPair(wcs, self.exposure.getWcs())
        for corner in lsst.afw.geom.Box2D(self.exposure.getBBox()).getCorners():
            bboxD.include(xyt.forwardTransform(lsst.afw.geom.Point2D(corner)))
        bboxI = lsst.afw.geom.Box2I(bboxD)
        result = TestDataset(bbox=bboxI, wcs=wcs, **kwds)
        oldCalib = self.exposure.getCalib()
        newCalib = result.exposure.getCalib()
        oldPsfShape = self.exposure.getPsf().computeShape()
        for record in self.catalog:
            if record.get(self.keys["nChild"]):
                raise NotImplementedError("Transforming blended sources in TestDatasets is not supported")
            magnitude = oldCalib.getMagnitude(record.get(self.keys["flux"]))
            newFlux = newCalib.getFlux(magnitude)
            oldCentroid = record.get(self.keys["centroid"])
            newCentroid = xyt.forwardTransform(oldCentroid)
            if record.get(self.keys["isStar"]):
                newDeconvolvedShape = None
            else:
                affine = xyt.linearizeForwardTransform(oldCentroid)
                oldFullShape = record.get(self.keys["shape"])
                oldDeconvolvedShape = lsst.afw.geom.ellipses.Quadrupole(
                    oldFullShape.getIxx() - oldPsfShape.getIxx(),
                    oldFullShape.getIyy() - oldPsfShape.getIyy(),
                    oldFullShape.getIxy() - oldPsfShape.getIxy(),
                    False
                )
                newDeconvolvedShape = oldDeconvolvedShape.transform(affine.getLinear())
            result.addSource(newFlux, newCentroid, newDeconvolvedShape)
        return result

    def realize(self, noise, schema):
        """!
        Create a simulated with noise and a simulated post-detection catalog with (Heavy)Footprints.

        @param[in]   noise      Standard deviation of noise to be added to the exposure.  The noise will be
                                Gaussian and constant, appropriate for the sky-limited regime.
        @param[in]   schema     Schema of the new catalog to be created.  Must start with self.schema (i.e.
                                schema.contains(self.schema) must be True), but typically contains fields for
                                already-configured measurement algorithms as well.

        @return a tuple of (exposure, catalog)
        """
        assert schema.contains(self.schema)
        mapper = lsst.afw.table.SchemaMapper(self.schema)
        mapper.addMinimalSchema(self.schema, True)
        exposure = self.exposure.clone()
        exposure.getMaskedImage().getVariance().getArray()[:, :] = noise**2
        exposure.getMaskedImage().getImage().getArray()[:, :] \
            += numpy.random.randn(exposure.getHeight(), exposure.getWidth())*noise
        catalog = lsst.afw.table.SourceCatalog(schema)
        catalog.extend(self.catalog, mapper=mapper)
        # Loop over sources and generate new HeavyFootprints that divide up the noisy pixels, not the
        # ideal no-noise pixels.
        for record in catalog:
            # parent objects have non-Heavy Footprints, which don't need to be updated after adding noise.
            if record.getParent() == 0:
                continue
            # get flattened arrays that correspond to the no-noise and noisy parent images
            parent = catalog.find(record.getParent())
            footprint = parent.getFootprint()
            parentFluxArrayNoNoise = numpy.zeros(footprint.getArea(), dtype=numpy.float32)
            lsst.afw.detection.flattenArray(
                footprint,
                self.exposure.getMaskedImage().getImage().getArray(),
                parentFluxArrayNoNoise,
                self.exposure.getXY0()
            )
            parentFluxArrayNoisy = numpy.zeros(footprint.getArea(), dtype=numpy.float32)
            lsst.afw.detection.flattenArray(
                footprint,
                exposure.getMaskedImage().getImage().getArray(),
                parentFluxArrayNoisy,
                exposure.getXY0()
            )
            oldHeavy = lsst.afw.detection.HeavyFootprintF.cast(record.getFootprint())
            fraction = (oldHeavy.getImageArray() / parentFluxArrayNoNoise)
            # n.b. this isn't a copy ctor - it's a copy from a vanilla Footprint, so it doesn't copy
            # the arrays we don't want to change, and hence we have to do that ourselves below.
            newHeavy = lsst.afw.detection.HeavyFootprintF(oldHeavy)
            newHeavy.getImageArray()[:] = parentFluxArrayNoisy*fraction
            newHeavy.getMaskArray()[:] = oldHeavy.getMaskArray()
            newHeavy.getVarianceArray()[:] = oldHeavy.getVarianceArray()
            record.setFootprint(newHeavy)
        return exposure, catalog


class AlgorithmTestCase(object):
    # Some tests depend on the noise realization in the test data or from the
    # numpy random number generator. In most cases, they are testing that the
    # measured flux lies within 2 sigma of the correct value, which we should
    # expect to fail sometimes. Some -- but sadly not all -- of these cases
    # have been marked with an "rng dependent" comment.
    #
    # We ensure these tests are provided with data which causes them to pass
    # by seeding the numpy RNG with this value. It can be over-ridden as
    # necessary in subclasses.
    randomSeed = 1234

    @classmethod
    def setUpClass(cls):
        numpy.random.seed(cls.randomSeed)

    def makeSingleFrameMeasurementConfig(self, plugin=None, dependencies=()):
        """Convenience function to create a Config instance for SingleFrameMeasurementTask

        The plugin and its dependencies will be the only plugins run, while the Centroid, Shape,
        and ModelFlux slots will be set to the truth fields generated by the TestDataset class.
        """
        config = SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = "truth"
        config.slots.shape = "truth"
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.instFlux = None
        config.slots.calibFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies)
        return config

    def makeSingleFrameMeasurementTask(self, plugin=None, dependencies=(), config=None, schema=None,
                                       algMetadata=None):
        """Convenience function to create a SingleFrameMeasurementTask with a simple configuration.
        """
        if config is None:
            if plugin is None:
                raise ValueError("Either plugin or config argument must not be None")
            config = self.makeSingleFrameMeasurementConfig(plugin=plugin, dependencies=dependencies)
        if schema is None:
            schema = TestDataset.makeMinimalSchema()
            # Clear all aliases so only those defined by config are set.
            schema.setAliasMap(None)
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        return SingleFrameMeasurementTask(schema=schema, algMetadata=algMetadata, config=config)

    def makeForcedMeasurementConfig(self, plugin=None, dependencies=()):
        """Convenience function to create a Config instance for ForcedMeasurementTask

        In addition to the plugins specified in the plugin and dependencies arguments,
        the TransformedCentroid and TransformedShape plugins will be run and used as the
        Centroid and Shape slots; these simply transform the reference catalog centroid
        and shape to the measurement coordinate system.
        """
        config = ForcedMeasurementTask.ConfigClass()
        config.slots.centroid = "base_TransformedCentroid"
        config.slots.shape = "base_TransformedShape"
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.instFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies) + ("base_TransformedCentroid",
                                                                  "base_TransformedShape")
        return config

    def makeForcedMeasurementTask(self, plugin=None, dependencies=(), config=None, refSchema=None,
                                  algMetadata=None):
        """Convenience function to create a ForcedMeasurementTask with a simple configuration.
        """
        if config is None:
            if plugin is None:
                raise ValueError("Either plugin or config argument must not be None")
            config = self.makeForcedMeasurementConfig(plugin=plugin, dependencies=dependencies)
        if refSchema is None:
            refSchema = TestDataset.makeMinimalSchema()
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        return ForcedMeasurementTask(refSchema=refSchema, algMetadata=algMetadata, config=config)


class TransformTestCase(object):
    """!
    Base class for testing measurement transformations.

    We test both that the transform itself operates successfully (fluxes are
    converted to magnitudes, flags are propagated properly) and that the
    transform is registered as the default for the appropriate measurement
    algorithms.

    In the simple case of one-measurement-per-transformation, the developer
    need not directly write any tests themselves: simply customizing the class
    variables is all that is required. More complex measurements (e.g.
    multiple aperture fluxes) require extra effort.
    """
    # The name used for the measurement algorithm; determines the names of the
    # fields in the resulting catalog. This default should generally be fine,
    # but subclasses can override if required.
    name = "MeasurementTransformTest"

    # These should be customized by subclassing.
    controlClass = None
    algorithmClass = None
    transformClass = None

    # Flags which may be set by the algorithm being tested. Can be customized
    # in subclasses.
    flagNames = ("flag",)

    # The plugin being tested should be registered under these names for
    # single frame and forced measurement. Should be customized by
    # subclassing.
    singleFramePlugins = ()
    forcedPlugins = ()

    def setUp(self):
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0), lsst.afw.geom.Point2I(200, 200))
        self.calexp = TestDataset.makeEmptyExposure(bbox)
        self._setupTransform()

    def tearDown(self):
        del self.calexp
        del self.inputCat
        del self.mapper
        del self.transform
        del self.outputCat

    def _populateCatalog(self, baseNames):
        records = []
        for flagValue in (True, False):
            records.append(self.inputCat.addNew())
            for baseName in baseNames:
                for flagName in self.flagNames:
                    if records[-1].schema.join(baseName, flagName) in records[-1].schema:
                        records[-1].set(records[-1].schema.join(baseName, flagName), flagValue)
        self._setFieldsInRecords(records, baseName)

    def _checkOutput(self, baseNames):
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            for baseName in baseNames:
                self._compareFieldsInRecords(inSrc, outSrc, baseName)
                for flagName in self.flagNames:
                    keyName = outSrc.schema.join(baseName, flagName)
                    if keyName in inSrc.schema:
                        self.assertEqual(outSrc.get(keyName), inSrc.get(keyName))
                    else:
                        self.assertFalse(keyName in outSrc.schema)

    def _runTransform(self, doExtend=True):
        if doExtend:
            self.outputCat.extend(self.inputCat, mapper=self.mapper)
        self.transform(self.inputCat, self.outputCat, self.calexp.getWcs(), self.calexp.getCalib())

    def testTransform(self, baseNames=None):
        """
        Test the operation of the transformation on a catalog containing random data.

        We check that:

        * An appropriate exception is raised on an attempt to transform between catalogs with different
          numbers of rows;
        * Otherwise, all appropriate conversions are properly appled and that flags have been propagated.

        The `baseNames` argument requires some explanation. This should be an iterable of the leading parts of
        the field names for each measurement; that is, everything that appears before `_flux`, `_flag`, etc.
        In the simple case of a single measurement per plugin, this is simply equal to `self.name` (thus
        measurements are stored as `self.name + "_flux"`, etc). More generally, the developer may specify
        whatever iterable they require. For example, to handle multiple apertures, we could have
        `(self.name + "_0", self.name + "_1", ...)`.

        @param[in]  baseNames  Iterable of the initial parts of measurement field names.
        """
        baseNames = baseNames or [self.name]
        self._populateCatalog(baseNames)
        self.assertRaises(lsst.pex.exceptions.LengthError, self._runTransform, False)
        self._runTransform()
        self._checkOutput(baseNames)

    def _checkRegisteredTransform(self, registry, name):
        # If this is a Python-based transform, we can compare directly; if
        # it's wrapped C++, we need to compare the wrapped class.
        self.assertEqual(registry[name].PluginClass.getTransformClass(), self.transformClass)

    def testRegistration(self):
        """
        Test that the transformation is appropriately registered with the relevant measurement algorithms.
        """
        for pluginName in self.singleFramePlugins:
            self._checkRegisteredTransform(lsst.meas.base.SingleFramePlugin.registry, pluginName)
        for pluginName in self.forcedPlugins:
            self._checkRegisteredTransform(lsst.meas.base.ForcedPlugin.registry, pluginName)


class SingleFramePluginTransformSetupHelper(object):

    def _setupTransform(self):
        self.control = self.controlClass()
        inputSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        # Trick algorithms that depend on the slot centroid or alias into thinking they've been defined;
        # it doesn't matter for this test since we won't actually use the plugins for anything besides
        # defining the schema.
        inputSchema.getAliasMap().set("slot_Centroid", "dummy")
        inputSchema.getAliasMap().set("slot_Shape", "dummy")
        self.algorithmClass(self.control, self.name, inputSchema)
        inputSchema.getAliasMap().erase("slot_Centroid")
        inputSchema.getAliasMap().erase("slot_Shape")
        self.inputCat = lsst.afw.table.SourceCatalog(inputSchema)
        self.mapper = lsst.afw.table.SchemaMapper(inputSchema)
        self.transform = self.transformClass(self.control, self.name, self.mapper)
        self.outputCat = lsst.afw.table.BaseCatalog(self.mapper.getOutputSchema())


class ForcedPluginTransformSetupHelper(object):

    def _setupTransform(self):
        self.control = self.controlClass()
        inputMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema(),
                                                  lsst.afw.table.SourceTable.makeMinimalSchema())
        # Trick algorithms that depend on the slot centroid or alias into thinking they've been defined;
        # it doesn't matter for this test since we won't actually use the plugins for anything besides
        # defining the schema.
        inputMapper.editOutputSchema().getAliasMap().set("slot_Centroid", "dummy")
        inputMapper.editOutputSchema().getAliasMap().set("slot_Shape", "dummy")
        self.algorithmClass(self.control, self.name, inputMapper, lsst.daf.base.PropertyList())
        inputMapper.editOutputSchema().getAliasMap().erase("slot_Centroid")
        inputMapper.editOutputSchema().getAliasMap().erase("slot_Shape")
        self.inputCat = lsst.afw.table.SourceCatalog(inputMapper.getOutputSchema())
        self.mapper = lsst.afw.table.SchemaMapper(inputMapper.getOutputSchema())
        self.transform = self.transformClass(self.control, self.name, self.mapper)
        self.outputCat = lsst.afw.table.BaseCatalog(self.mapper.getOutputSchema())


class FluxTransformTestCase(TransformTestCase):

    def _setFieldsInRecords(self, records, name):
        for record in records:
            record[record.schema.join(name, 'flux')] = numpy.random.random()
            record[record.schema.join(name, 'fluxSigma')] = numpy.random.random()

        # Negative fluxes should be converted to NaNs.
        assert len(records) > 1
        records[0][record.schema.join(name, 'flux')] = -1

    def _compareFieldsInRecords(self, inSrc, outSrc, name):
        fluxName, fluxSigmaName = inSrc.schema.join(name, 'flux'), inSrc.schema.join(name, 'fluxSigma')
        if inSrc[fluxName] > 0:
            mag, magErr = self.calexp.getCalib().getMagnitude(inSrc[fluxName], inSrc[fluxSigmaName])
            self.assertEqual(outSrc[outSrc.schema.join(name, 'mag')], mag)
            self.assertEqual(outSrc[outSrc.schema.join(name, 'magErr')], magErr)
        else:
            self.assertTrue(numpy.isnan(outSrc[outSrc.schema.join(name, 'mag')]))
            self.assertTrue(numpy.isnan(outSrc[outSrc.schema.join(name, 'magErr')]))


class CentroidTransformTestCase(TransformTestCase):

    def _setFieldsInRecords(self, records, name):
        for record in records:
            record[record.schema.join(name, 'x')] = numpy.random.random()
            record[record.schema.join(name, 'y')] = numpy.random.random()
            # Some algorithms set no errors; some set only sigma on x & y; some provide
            # a full covariance matrix. Set only those which exist in the schema.
            for fieldSuffix in ('xSigma', 'ySigma', 'x_y_Cov'):
                fieldName = record.schema.join(name, fieldSuffix)
                if fieldName in record.schema:
                    record[fieldName] = numpy.random.random()

    def _compareFieldsInRecords(self, inSrc, outSrc, name):
        centroidResultKey = CentroidResultKey(inSrc.schema[self.name])
        centroidResult = centroidResultKey.get(inSrc)

        coord = lsst.afw.table.CoordKey(outSrc.schema[self.name]).get(outSrc)
        coordTruth = self.calexp.getWcs().pixelToSky(centroidResult.getCentroid())
        self.assertEqual(coordTruth, coord)

        # If the centroid has an associated uncertainty matrix, the coordinate
        # must have one too, and vice versa.
        try:
            coordErr = lsst.afw.table.CovarianceMatrix2fKey(outSrc.schema[self.name],
                                                            ["ra", "dec"]).get(outSrc)
        except lsst.pex.exceptions.NotFoundError:
            self.assertFalse(centroidResultKey.getCentroidErr().isValid())
        else:
            transform = self.calexp.getWcs().linearizePixelToSky(coordTruth, lsst.afw.geom.radians)
            coordErrTruth = numpy.dot(numpy.dot(transform.getLinear().getMatrix(),
                                                centroidResult.getCentroidErr()),
                                      transform.getLinear().getMatrix().transpose())
            numpy.testing.assert_array_almost_equal(numpy.array(coordErrTruth), coordErr)
