from builtins import str
from builtins import object
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
import math

import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pex.config

__all__ = ("NoiseReplacerConfig", "NoiseReplacer", "DummyNoiseReplacer")


class NoiseReplacerConfig(lsst.pex.config.Config):
    noiseSource = lsst.pex.config.ChoiceField(
        doc='How to choose mean and variance of the Gaussian noise we generate?',
        dtype=str,
        allowed={
            'measure': 'Measure clipped mean and variance from the whole image',
            'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
            'variance': "Mean = 0, variance = the image's variance",
        },
        default='measure', optional=False
    )
    noiseOffset = lsst.pex.config.Field(
        dtype=float, optional=False, default=0.,
        doc='Add ann offset to the generated noise.'
    )
    noiseSeedMultiplier = lsst.pex.config.Field(
        dtype=int, default=1,
        doc='The seed multiplier value to use for random number generation.  0 will not set seed.'
    )


class NoiseReplacer(object):
    """!
    Class that handles replacing sources with noise during measurement.

    When measuring a source (or the children associated with a parent source), this class is used
    to replace its neighbors with noise, using the deblender's definition of the sources as stored
    in HeavyFootprints attached to the SourceRecords.  The algorithm works as follows:
     - We start by replacing all pixels that are in source Footprints with artificially
       generated noise (__init__).
     - When we are about to measure a particular source, we add it back in, by inserting that source's
       HeavyFootprint (from the deblender) into the image.
     - When we are done measuring that source, we again replace the HeavyFootprint with (the same)
       artificial noise.
     - After measuring all sources, we return the image to its original state.

    This is a functional copy of the code in the older ReplaceWithNoiseTask, but with a slightly different
    API needed for the new measurement framework; note that it is not a Task, as the lifetime of a
    NoiseReplacer now corresponds to a single exposure, not an entire processing run.
    """

    ConfigClass = NoiseReplacerConfig

    def __init__(self, config, exposure, footprints, noiseImage=None, exposureId=None, log=None):
        """!
        Initialize the NoiseReplacer.

        @param[in]      config       instance of NoiseReplacerConfig
        @param[in,out]  exposure     Exposure to be noise replaced. (All sources replaced on return)
        @param[in]      footprints   dict of {id: (parent, footprint)};
        @param[in]      noiseImage   an afw.image.ImageF used as a predictable noise replacement source
                                     (for tests only)
        @param[in]      log          Log object to use for status messages; no status messages
                                     will be printed if None

        'footprints' is a dict of {id: (parent, footprint)}; when used in SFM, the ID will be the
        source ID, but in forced photometry, this will be the reference ID, as that's what we used to
        determine the deblend families.  This routine should create HeavyFootprints for any non-Heavy
        Footprints, and replace them in the dict.  It should then create a dict of HeavyFootprints
        containing noise, but only for parent objects, then replace all sources with noise.
        This should ignore any footprints that lay outside the bounding box of the exposure,
        and clip those that lie on the border.

        NOTE: as the code currently stands, the heavy footprint for a deblended object must be available
        from the input catalog.  If it is not, it cannot be reproduced here.  In that case, the
        topmost parent in the objects parent chain must be used.  The heavy footprint for that source
        is created in this class from the masked image.
        """
        noiseMeanVar = None
        self.noiseSource = config.noiseSource
        self.noiseOffset = config.noiseOffset
        self.noiseSeedMultiplier = config.noiseSeedMultiplier
        self.noiseGenMean = None
        self.noiseGenStd = None
        self.log = log

        # creates heavies, replaces all footprints with noise
        # We need the source table to be sorted by ID to do the parent lookups
        self.exposure = exposure
        self.footprints = footprints
        mi = exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        # Add temporary Mask planes for THISDET and OTHERDET
        self.removeplanes = []
        bitmasks = []
        for maskname in ['THISDET', 'OTHERDET']:
            try:
                # does it already exist?
                plane = mask.getMaskPlane(maskname)
                if self.log:
                    self.log.debug('Mask plane "%s" already existed', maskname)
            except:
                # if not, add it; we should delete it when done.
                plane = mask.addMaskPlane(maskname)
                self.removeplanes.append(maskname)
            mask.clearMaskPlane(plane)
            bitmask = mask.getPlaneBitMask(maskname)
            bitmasks.append(bitmask)
            if self.log:
                self.log.debug('Mask plane "%s": plane %i, bitmask %i = 0x%x',
                               maskname, plane, bitmask, bitmask)
        self.thisbitmask, self.otherbitmask = bitmasks
        del bitmasks
        self.heavies = {}
        # Start by creating HeavyFootprints for each source which has no parent
        # and just use them for children which do not already have heavy footprints.
        # If a heavy footprint is available for a child, we will use it. Otherwise,
        # we use the first parent in the parent chain which has a heavy footprint,
        # which with the one level deblender will alway be the topmost parent
        # NOTE: heavy footprints get destroyed by the transform process in forcedPhotImage.py,
        # so they are never available for forced measurements.

        # Create in the dict heavies = {id:heavyfootprint}
        for id, fp in footprints.items():
            if fp[1].isHeavy():
                self.heavies[id] = afwDet.cast_HeavyFootprintF(fp[1])
            elif fp[0] == 0:
                self.heavies[id] = afwDet.makeHeavyFootprint(fp[1], mi)

        ### FIXME: the heavy footprint includes the mask
        ### and variance planes, which we shouldn't need
        ### (I don't think we ever want to modify them in
        ### the input image).  Copying them around is
        ### wasteful.

        # We now create a noise HeavyFootprint for each source with has a heavy footprint.
        # We'll put the noise footprints in a dict heavyNoise = {id:heavyNoiseFootprint}
        self.heavyNoise = {}
        noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar, exposureId=exposureId)
        #  The noiseGenMean and Std are used by the unit tests
        self.noiseGenMean = noisegen.mean
        self.noiseGenStd = noisegen.std
        if self.log:
            self.log.debug('Using noise generator: %s', str(noisegen))
        for id in self.heavies:
            fp = footprints[id][1]
            noiseFp = noisegen.getHeavyFootprint(fp)
            self.heavyNoise[id] = noiseFp
            # Also insert the noisy footprint into the image now.
            # Notice that we're just inserting it into "im", ie,
            # the Image, not the MaskedImage.
            noiseFp.insert(im)
            # Also set the OTHERDET bit
            afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def insertSource(self, id):
        """!
        Insert the heavy footprint of a given source into the exposure

        @param[in]  id   id for current source to insert from original footprint dict

        Also adjusts the mask plane to show the source of this footprint.
        """
        # Copy this source's pixels into the image
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        # usedid can point either to this source, or to the first parent in the
        # parent chain which has a heavy footprint (or to the topmost parent,
        # which always has one)
        usedid = id
        while self.footprints[usedid][0] != 0 and usedid not in self.heavies:
            usedid = self.footprints[usedid][0]
        fp = self.heavies[usedid]
        fp.insert(im)
        afwDet.setMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.clearMaskFromFootprint(mask, fp, self.otherbitmask)

    def removeSource(self, id):
        """!
        Remove the heavy footprint of a given source and replace with previous noise

        @param[in]  id   id for current source to insert from original footprint dict

        Also restore the mask plane.
        """
        # remove a single source
        # (Replace this source's pixels by noise again.)
        # Do this by finding the source's top-level ancestor
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()

        # use the same algorithm as in remove Source to find the heavy noise footprint
        # which will undo what insertSource(id) does
        usedid = id
        while self.footprints[usedid][0] != 0 and usedid not in self.heavies:
            usedid = self.footprints[usedid][0]
        # Re-insert the noise pixels
        fp = self.heavyNoise[usedid]
        fp.insert(im)
        # Clear the THISDET mask plane.
        afwDet.clearMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def end(self):
        """!
        End the NoiseReplacer.

        Restore original data to the exposure from the heavies dictionary
        Restore the mask planes to their original state
        """
        # restores original image, cleans up temporaries
        # (ie, replace all the top-level pixels)
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        for id in self.footprints.keys():
            if self.footprints[id][0] != 0:
                continue
            self.heavies[id].insert(im)
        for maskname in self.removeplanes:
            mask.removeAndClearMaskPlane(maskname, True)

        del self.removeplanes
        del self.thisbitmask
        del self.otherbitmask
        del self.heavies
        del self.heavyNoise

    def getNoiseGenerator(self, exposure, noiseImage, noiseMeanVar, exposureId=None):
        """!
        Generate noise image using parameters given
        """
        if noiseImage is not None:
            return ImageNoiseGenerator(noiseImage)
        rand = None
        if self.noiseSeedMultiplier:
            # default plugin, our seed
            if exposureId is not None and exposureId != 0:
                seed = exposureId * self.noiseSeedMultiplier
            else:
                seed = self.noiseSeedMultiplier
            rand = afwMath.Random(afwMath.Random.MT19937, seed)
        if noiseMeanVar is not None:
            try:
                # Assume noiseMeanVar is an iterable of floats
                noiseMean, noiseVar = noiseMeanVar
                noiseMean = float(noiseMean)
                noiseVar = float(noiseVar)
                noiseStd = math.sqrt(noiseVar)
                if self.log:
                    self.log.debug('Using passed-in noise mean = %g, variance = %g -> stdev %g',
                                   noiseMean, noiseVar, noiseStd)
                return FixedGaussianNoiseGenerator(noiseMean, noiseStd, rand=rand)
            except:
                if self.log:
                    self.log.debug('Failed to cast passed-in noiseMeanVar to floats: %s',
                                   str(noiseMeanVar))
        offset = self.noiseOffset
        noiseSource = self.noiseSource

        if noiseSource == 'meta':
            # check the exposure metadata
            meta = exposure.getMetadata()
            # this key name correspond to SubtractBackgroundTask() in meas_algorithms
            try:
                bgMean = meta.getAsDouble('BGMEAN')
                # We would have to adjust for GAIN if ip_isr didn't make it 1.0
                noiseStd = math.sqrt(bgMean)
                if self.log:
                    self.log.debug('Using noise variance = (BGMEAN = %g) from exposure metadata',
                                   bgMean)
                return FixedGaussianNoiseGenerator(offset, noiseStd, rand=rand)
            except:
                if self.log:
                    self.log.debug('Failed to get BGMEAN from exposure metadata')

        if noiseSource == 'variance':
            if self.log:
                self.log.debug('Will draw noise according to the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset, rand=rand)

        # Compute an image-wide clipped variance.
        im = exposure.getMaskedImage().getImage()
        s = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        noiseMean = s.getValue(afwMath.MEANCLIP)
        noiseStd = s.getValue(afwMath.STDEVCLIP)
        if self.log:
            self.log.debug("Measured from image: clipped mean = %g, stdev = %g",
                           noiseMean, noiseStd)
        return FixedGaussianNoiseGenerator(noiseMean + offset, noiseStd, rand=rand)


class NoiseReplacerList(list):
    """Syntactic sugar that makes a list of NoiseReplacers (for multiple exposures)
    behave like a single one.

    This is only used in the multifit driver, but the logic there is already pretty
    complex, so it's nice to have this to simplify it.
    """

    def __init__(self, exposuresById, footprintsByExp):
        # exposuresById --- dict of {exposureId: exposure} (possibly subimages)
        # footprintsByExp --- nested dict of {exposureId: {objId: (parent, footprint)}}
        list.__init__(self)
        for expId, exposure in exposuresById.items():
            self.append(NoiseReplacer(exposure, footprintsByExp[expId]), expId)

    def insertSource(self, id):
        """Insert the original pixels for a given source (by id) into the original exposure.
        """
        for item in self:
            self.insertSource(id)

    def removeSource(self, id):
        """Insert the noise pixels for a given source (by id) into the original exposure.
        """
        for item in self:
            self.removeSource(id)

    def end(self):
        """Cleanup when the use of the Noise replacer is done.
        """
        for item in self:
            self.end()


class NoiseGenerator(object):
    """!
    Base class for noise generators used by the "doReplaceWithNoise" routine:
    these produce HeavyFootprints filled with noise generated in various ways.

    This is an abstract base class.
    """

    def getHeavyFootprint(self, fp):
        bb = fp.getBBox()
        mim = self.getMaskedImage(bb)
        return afwDet.makeHeavyFootprint(fp, mim)

    def getMaskedImage(self, bb):
        im = self.getImage(bb)
        return afwImage.MaskedImageF(im)

    def getImage(self, bb):
        return None


class ImageNoiseGenerator(NoiseGenerator):
    """
    Generates noise by cutting out a subimage from a user-supplied noise Image.
    """

    def __init__(self, img):
        """
        img: an afwImage.ImageF
        """
        self.mim = afwImage.MaskedImageF(img)
        self.mean = afwMath.makeStatistics(img, afwMath.MEAN)
        self.std = afwMath.makeStatistics(img, afwMath.STDEV)

    def getMaskedImage(self, bb):
        return self.mim


class GaussianNoiseGenerator(NoiseGenerator):
    """!
    Generates noise using the afwMath.Random() and afwMath.randomGaussianImage() routines.

    This is an abstract base class.
    """

    def __init__(self, rand=None):
        if rand is None:
            rand = afwMath.Random()
        self.rand = rand

    def getRandomImage(self, bb):
        # Create an Image and fill it with Gaussian noise.
        rim = afwImage.ImageF(bb.getWidth(), bb.getHeight())
        rim.setXY0(bb.getMinX(), bb.getMinY())
        afwMath.randomGaussianImage(rim, self.rand)
        return rim


class FixedGaussianNoiseGenerator(GaussianNoiseGenerator):
    """!
    Generates Gaussian noise with a fixed mean and standard deviation.
    """

    def __init__(self, mean, std, rand=None):
        super(FixedGaussianNoiseGenerator, self).__init__(rand=rand)
        self.mean = mean
        self.std = std

    def __str__(self):
        return 'FixedGaussianNoiseGenerator: mean=%g, std=%g' % (self.mean, self.std)

    def getImage(self, bb):
        rim = self.getRandomImage(bb)
        rim *= self.std
        rim += self.mean
        return rim


class VariancePlaneNoiseGenerator(GaussianNoiseGenerator):
    """!
    Generates Gaussian noise whose variance matches that of the variance plane of the image.
    """

    def __init__(self, var, mean=None, rand=None):
        """
        var: an afwImage.ImageF; the variance plane.
        mean: floating-point or afwImage.Image
        """
        super(VariancePlaneNoiseGenerator, self).__init__(rand=rand)
        self.var = var
        if mean is not None and mean == 0.:
            mean = None
        self.mean = mean

    def __str__(self):
        return 'VariancePlaneNoiseGenerator: mean=' + str(self.mean)

    def getImage(self, bb):
        rim = self.getRandomImage(bb)
        # Use the image's variance plane to scale the noise.
        stdev = afwImage.ImageF(self.var, bb, afwImage.LOCAL, True)
        stdev.sqrt()
        rim *= stdev
        if self.mean is not None:
            rim += self.mean
        return rim


class DummyNoiseReplacer(object):
    """!
    A do-nothing standin for NoiseReplacer, used when we want to disable NoiseReplacer

    DummyNoiseReplacer has all the public methods of NoiseReplacer, but none of them do anything.
    """

    def insertSource(self, id):
        pass

    def removeSource(self, id):
        pass

    def end(self):
        pass
