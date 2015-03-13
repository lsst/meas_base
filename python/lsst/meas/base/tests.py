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
import lsst.afw.geom.ellipses
import lsst.afw.coord
import lsst.utils.tests

from .sfm import SingleFrameMeasurementTask

class MakeTestData(object):
    @staticmethod
    def drawGaussian(image, flux, ellipse):
        bbox = image.getBBox()
        x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX()),
                              numpy.arange(bbox.getBeginY(), bbox.getEndY()))
        t = ellipse.getGridTransform()
        xt = t[t.XX] * x + t[t.XY] * y + t[t.X]
        yt = t[t.YX] * y + t[t.YY] * y + t[t.Y]
        image.getArray()[:,:] = numpy.exp(-0.5*(xt**2 + yt**2))
        image.getArray()[:,:] *= flux / image.getArray().sum()

    @staticmethod
    def makeCatalog():
        """Create and return the truth catalog and bounding box for the simulation.

        The simulated records will be in image coordinates, and no footprints will be attached.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        nChildKey = schema.addField("deblend_nchild", type=int)
        xKey = schema.addField("truth_x", type=float,
                               doc="true simulated centroid x", units="pixels")
        yKey = schema.addField("truth_y", type=float,
                               doc="true simulated centroid y", units="pixels")
        centroidKey = lsst.afw.table.Point2DKey(xKey, yKey)
        fluxKey = schema.addField("truth_flux", type=float, doc="true flux", units="dn")
        xxKey = schema.addField("truth_xx", type=float,
                                doc="true shape after PSF convolution", units="pixels^2")
        yyKey = schema.addField("truth_yy", type=float,
                                doc="true shape after PSF convolution", units="pixels^2")
        xyKey = schema.addField("truth_xy", type=float,
                                doc="true shape after PSF convolution", units="pixels^2")
        shapeKey = lsst.afw.table.QuadrupoleKey(xxKey, yyKey, xyKey)
        starFlagKey = schema.addField("truth_isStar", type="Flag", doc="set if the object is a star")
        schema.setVersion(1)
        catalog = lsst.afw.table.SourceCatalog(schema)
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0), lsst.afw.geom.Point2I(200, 200))
        # a bright, isolated star near (50, 50)
        record = catalog.addNew()
        record.set(nChildKey, 0)
        record.set(centroidKey, lsst.afw.geom.Point2D(50.1, 49.8))
        record.set(fluxKey, 100000.0)
        record.set(shapeKey, lsst.afw.geom.ellipses.Quadrupole(0.0, 0.0, 0.0))
        record.set(starFlagKey, True)
        # a moderately-resolved, isolated galaxy near (150, 50)
        record = catalog.addNew()
        record.set(nChildKey, 0)
        record.set(centroidKey, lsst.afw.geom.Point2D(150.3, 50.2))
        record.set(fluxKey, 75000.0)
        record.set(shapeKey, lsst.afw.geom.ellipses.Quadrupole(8.0, 6.0, 0.5))
        record.set(starFlagKey, False)
        # a blend of a star and galaxy near (95, 150) and (105, 150)
        parent = catalog.addNew()
        parent.set(nChildKey, 2)
        record1 = catalog.addNew()
        record1.set(nChildKey, 0)
        record1.setParent(parent.getId())
        record1.set(centroidKey, lsst.afw.geom.Point2D(95.4, 149.9))
        record1.set(fluxKey, 80000.0)
        record1.set(shapeKey, lsst.afw.geom.ellipses.Quadrupole(0.0, 0.0, 0.0))
        record1.set(starFlagKey, True)
        record2 = catalog.addNew()
        record2.set(nChildKey, 0)
        record2.setParent(parent.getId())
        record2.set(centroidKey, lsst.afw.geom.Point2D(104.8, 150.2))
        record2.set(fluxKey, 120000.0)
        record2.set(shapeKey, lsst.afw.geom.ellipses.Quadrupole(10.0, 7.0, -0.5))
        record2.set(starFlagKey, False)
        parent.set(fluxKey, record1.get(fluxKey) + record2.get(fluxKey))
        parent.set(starFlagKey, False)
        parent.set(
            centroidKey,
            lsst.afw.geom.Point2D(
                (lsst.afw.geom.Extent2D(record1.get(centroidKey)) * record1.get(fluxKey)
                 + lsst.afw.geom.Extent2D(record2.get(centroidKey)) * record2.get(fluxKey))
                / parent.get(fluxKey)
                )
            )
        # we don't bother setting the truth values for parent's shape, since we don't need them
        # for sims and don't expect to be able to measure them well anyway.
        return catalog, bbox

    @staticmethod
    def makeEmptyExposure(bbox, crval=None, psfSigma=2.0, psfDim=17, fluxMag0=1E12):
        """Create an Exposure, with a Calib, Wcs, and Psf, but no pixel values set.

        @param[in]       bbox        Bounding box of the image (image coordinates) as returned by makeCatalog.
        @param[in]       crval       afw.coord.Coord: center of the TAN WCS attached to the image.
        @param[in]       psfSigma    Radius (sigma) of the Gaussian PSF attached to the image
        @param[in]       psfDim      Width and height of the image's Gaussian PSF attached to the image
        @param[in]       fluxMag0    Flux at magnitude zero (in e-) used to set the Calib of the exposure.
        """
        if crval is None:
            crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        exposure = lsst.afw.image.ExposureF(bbox)
        crpix = lsst.afw.geom.Box2D(bbox).getCenter()
        cdelt = (0.2 * lsst.afw.geom.arcseconds).asDegrees()
        wcs = lsst.afw.image.makeWcs(crval, crpix, cdelt, 0.0, 0.0, cdelt)
        psf = lsst.afw.detection.GaussianPsf(psfDim, psfDim, psfSigma)
        calib = lsst.afw.image.Calib()
        calib.setFluxMag0(fluxMag0)
        exposure.setWcs(wcs)
        exposure.setPsf(psf)
        exposure.setCalib(calib)
        return exposure

    @staticmethod
    def fillImages(catalog, exposure, threshold=None, noise=100.0):
        """Fill a simulated Exposure from makeEmptyExposure() with the objects defined by a truth catalog from
        makeTruthCatalog().  Also attaches Footprints to the SourceRecords in the catalog, sets the
        'coord' field, and fills in the detection mask planes.  The Footprints will be regular Footprints
        for non-blended objects, and HeavyFootprints for deblended objects, intended to mimic
        better-than-realistic outputs from detection and deblending.

        @param[in,out]   catalog     SourceCatalog created by makeCatalog.  HeavyFootprints and 'coord' will
                                     be filled on return, and shapes will be convolved by the PSF shape.
        @param[in,out]   exposure    ExposureF to fill, as created by makeEmptyExposure()
        @param[in]       threshold   afw.detection.Threshold object used to determine the size of the
                                     HeavyFootprints attached to the SourceCatalog (after thresholding, we'll
                                     also grow the Footprints by the PSF sigma).
        @param[in]       noise       Standard deviation of Gaussian noise added to image (constant across the
                                     image; appropriate for sky-dominated limit).
        """
        if threshold is None:
            threshold = lsst.afw.detection.Threshold(10.0, lsst.afw.detection.Threshold.VALUE)

        psf = lsst.afw.detection.GaussianPsf.cast(exposure.getPsf())
        schema = catalog.schema
        nChildKey = schema.find("deblend_nchild").key
        centroidKey = lsst.afw.table.Point2DKey(schema.find("truth_x").key, schema.find("truth_y").key)
        shapeKey = lsst.afw.table.QuadrupoleKey(schema.find("truth_xx").key, schema.find("truth_yy").key,
                                                schema.find("truth_xy").key)
        fluxKey = schema.find("truth_flux").key

        # First pass: generate full-size images, each containing a single object, and add them into
        # the Exposure
        images = {record.getId(): lsst.afw.image.ImageF(exposure.getBBox())
                  for record in catalog}
        for record in catalog:
            if record.get(nChildKey) != 0: continue
            ellipse = lsst.afw.geom.ellipses.Ellipse(record.get(shapeKey).convolve(psf.computeShape()),
                                                     record.get(centroidKey))
            record.set(shapeKey, ellipse.getCore())
            MakeTestData.drawGaussian(images[record.getId()], record.get(fluxKey), ellipse)
            exposure.getMaskedImage().getImage().getArray()[:,:] += images[record.getId()].getArray()[:,:]
            if record.getParent() != 0:  # parent images come from combining child images
                images[record.getParent()] += images[record.getId()]

        exposure.getMaskedImage().getVariance().getArray()[:,:] = noise**2
        exposure.getMaskedImage().getImage().getArray()[:,:] \
            += numpy.random.randn(exposure.getHeight(), exposure.getWidth())*noise

        # Second pass: detect on single-object images to generate better-than-reality Footprints
        for record in catalog:
            # this detection doesn't really match what we'd do on real data (it can't, because we don't
            # have noise yet), but it seems to yield reasonably-sized footprints.
            fpSet = lsst.afw.detection.FootprintSet(images[record.getId()], threshold)
            fpSet = lsst.afw.detection.FootprintSet(fpSet, int(psf.getSigma()+1.0), True)
            if len(fpSet.getFootprints()) != 1:
                raise ValueError("Threshold value results in multiple Footprints for a single object")
            fpSet.setMask(exposure.getMaskedImage().getMask(), "DETECTED")
            record.setFootprint(fpSet.getFootprints()[0])

        # Third pass: for all parent objects, make "perfectly deblended" child HeavyFootprints using the
        # true, noise-free single-object images to divide the flux.
        for parent in catalog.getChildren(0):
            print "Processing parent object", parent.getId()
            for child in catalog.getChildren(parent.getId()):
                print "Processing child object", child.getId()
                fraction = images[child.getId()].getArray() / images[parent.getId()].getArray()
                deblend = lsst.afw.image.MaskedImageF(exposure.getMaskedImage(), True)
                deblend.getImage().getArray()[:,:] *= fraction
                child.setFootprint(lsst.afw.detection.HeavyFootprintF(child.getFootprint(), deblend))
        return images

class AlgorithmTestCase(lsst.utils.tests.TestCase):
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

    def setUp(self):
        catalog, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(catalog, exposure)
        schema = catalog.schema
        self.truth = catalog
        self.calexp = exposure
        self.centroidKey = lsst.afw.table.Point2DKey(schema.find("truth_x").key, schema.find("truth_y").key)
        self.shapeKey = lsst.afw.table.QuadrupoleKey(schema.find("truth_xx").key, schema.find("truth_yy").key,
                                                     schema.find("truth_xy").key)
        self.fluxKey = schema.find("truth_flux").key

    def tearDown(self):
        del self.centroidKey
        del self.shapeKey
        del self.fluxKey
        del self.truth
        del self.calexp

    def runSingleFrameMeasurementTask(self, plugin, dependencies=(), config=None):
        if config is None:
            config = SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = None
        config.slots.shape = None
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.instFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies)
        schemaMapper = lsst.afw.table.SchemaMapper(self.truth.schema)
        schemaMapper.addMinimalSchema(self.truth.schema)
        algMetadata = lsst.daf.base.PropertyList()
        task = SingleFrameMeasurementTask(schema=schemaMapper.editOutputSchema(), algMetadata=algMetadata,
                                          config=config)
        measCat = lsst.afw.table.SourceCatalog(task.schema)
        measCat.table.setMetadata(algMetadata)
        measCat.extend(self.truth, schemaMapper)
        measCat.getTable().defineModelFlux("truth")
        measCat.getTable().defineCentroid("truth")
        measCat.getTable().defineShape("truth")
        task.run(measCat, self.calexp)
        return measCat

class TransformTestCase(lsst.utils.tests.TestCase):
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
        self.calexp = MakeTestData.makeEmptyExposure(MakeTestData.makeCatalog()[1])
        self._setupTransform()

    def tearDown(self):
        del self.calexp
        del self.inputCat
        del self.mapper
        del self.transform
        del self.outputCat

    def _setupTransform(self):
        self.control = self.controlClass()
        inputSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        algo = self.algorithmClass(self.control, self.name, inputSchema)
        self.inputCat = lsst.afw.table.SourceCatalog(inputSchema)
        self.mapper = lsst.afw.table.SchemaMapper(inputSchema)
        self.transform = self.transformClass(self.control, self.name, self.mapper)
        self.outputCat = lsst.afw.table.BaseCatalog(self.mapper.getOutputSchema())

    def _populateCatalog(self, baseNames):
        for flagValue in (True, False):
            record = self.inputCat.addNew()
            for baseName in baseNames:
                record[baseName + '_flux'], record[baseName + '_fluxSigma'] = numpy.random.random(2)
                for flagName in self.flagNames:
                    if baseName + '_' + flagName in record.schema:
                        record.set(baseName + '_' + flagName, flagValue)

    def _checkOutput(self, baseNames):
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            for baseName in baseNames:
                fluxName = baseName + '_flux'
                fluxSigmaName = baseName + '_fluxSigma'
                mag, magErr = self.calexp.getCalib().getMagnitude(inSrc[fluxName], inSrc[fluxSigmaName])
                self.assertEqual(outSrc[baseName + '_mag'], mag)
                self.assertEqual(outSrc[baseName + '_magErr'], magErr)
                for flagName in self.flagNames:
                    keyName = baseName + '_' + flagName
                    if keyName in inSrc.schema:
                        self.assertEqual(outSrc.get(keyName), inSrc.get(keyName))
                    else:
                        self.assertFalse(keyName in outSrc.schema)

    def _runTransform(self):
        self.outputCat.extend(self.inputCat, mapper=self.mapper)
        self.transform(self.inputCat, self.outputCat, self.calexp.getWcs(), self.calexp.getCalib())

    def testTransform(self, baseNames=None):
        """
        Test the operation of the transformation on a catalog containing random data.

        We check that appropriate conversions have been properly applied and that flags have been propagated.

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
        self._runTransform()
        self._checkOutput(baseNames)

    def _checkRegisteredTransform(self, registry, name):
        self.assertEqual(registry[name].PluginClass.getTransformClass(), self.transformClass)

    def testRegistration(self):
        """
        Test that the transformation is appropriately registered with the relevant measurement algorithms.
        """
        for pluginName in self.singleFramePlugins:
            self._checkRegisteredTransform(lsst.meas.base.SingleFramePlugin.registry, pluginName)
        for pluginName in self.forcedPlugins:
            self._checkRegisteredTransform(lsst.meas.base.ForcedPlugin.registry, pluginName)
