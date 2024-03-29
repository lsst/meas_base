# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math

import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config

__all__ = ("NoiseReplacerConfig", "NoiseReplacer", "DummyNoiseReplacer")


class NoiseReplacerConfig(lsst.pex.config.Config):
    """Noise replacement configuration."""

    noiseSource = lsst.pex.config.ChoiceField(
        doc='How to choose mean and variance of the Gaussian noise we generate?',
        dtype=str,
        allowed={
            'measure': 'Measure clipped mean and variance from the whole image',
            'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
            'variance': "Mean = 0, variance = the image's variance",
            'variance_median': "Mean = 0, variance = median(variance plane)"
        },
        default='measure', optional=False
    )
    noiseOffset = lsst.pex.config.Field(
        doc='Add ann offset to the generated noise.',
        dtype=float, optional=False, default=0.0
    )
    noiseSeedMultiplier = lsst.pex.config.Field(
        dtype=int, default=1,
        doc="The seed multiplier value to use for random number generation:\n"
            ">= 1: set the seed deterministically based on exposureId\n"
            "0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)"
    )


class NoiseReplacer:
    r"""Replace sources with noise during measurement.

    Parameters
    ----------
    config : `NoiseReplacerConfig`
        Configuration.
    exposure : `lsst.afw.image.Exposure`
        Image in which sources will be replaced by noise. During operation,
        the image will be modified in-place to replace all sources. At the end
        of the measurment procedure, the original sources will be replaced.
    footprints : `dict`
        Mapping of ``id`` to a tuple of ``(parent, Footprint)``. When used in
        single-frame measurement, ``id`` is the source ID, but in forced
        photometry this is the reference ID (as that is used to determine
        deblend families).
    noiseImage : `lsst.afw.image.ImageF`
        An image used as a predictable noise replacement source. Used during
        testing only.
    log : `lsst.log.log.log.Log` or `logging.Logger`, optional
        Logger to use for status messages; no status messages will be recorded
        if `None`.

    Notes
    -----
    When measuring a source (or the children associated with a parent source),
    this class is used to replace its neighbors with noise, using the
    deblender's definition of the sources as stored in
    `~lsst.afw.detection.heavyFootprint.HeavyFootprint`\ s attached to the
    `~lsst.afw.table.SourceRecord`\ s.  The algorithm works as follows:

    #. All pixels in the source `~lsst.afw.detection.Footprint`\ s are replaced
       with artificially generated noise (in `NoiseReplacer.__init__`).
    #. Before each source is measured, we restore the original pixel data by
       inserting that source's
       `~lsst.afw.detection.heavyFootprint.HeavyFootprint` (from the deblender)
       into the image.
    #. After measurement, we again replace the source pixels with (the same)
       artificial noise.
    #. After measuring all sources, the image is returned to its original
       state.

    This is a functional copy of the code in the older
    ``ReplaceWithNoiseTask``, but with a slightly different API needed for the
    new measurement framework; note that it is not an `~lsst.pipe.base.Task`,
    as the lifetime of a ``NoiseReplacer`` now corresponds to a single
    exposure, not an entire processing run.

    When processing the ``footprints`` parameter, this routine should create
    `~lsst.afw.detection.heavyFootprint.HeavyFootprint`\ s for any non-Heavy
    `~lsst.afw.detection.Footprint`\ s, and replace them in the dictionary. It
    should then create a dict of
    `~lsst.afw.detection.heavyFootprint.HeavyFootprint`\ s containing noise,
    but only for parent objects, then replace all sources with noise. This
    should ignore any footprints that lay outside the bounding box of the
    exposure, and clip those that lie on the border.

    As the code currently stands, the heavy footprint for a deblended object
    must be available from the input catalog.  If it is not, it cannot be
    reproduced here. In that case, the topmost parent in the objects parent
    chain must be used. The heavy footprint for that source is created in
    this class from the masked image.
    """

    ConfigClass = NoiseReplacerConfig

    exposure = None
    """Image on which the NoiseReplacer is operating (`lsst.afw.image.Exposure`).
    """

    footprints = None
    """Mapping of ``id`` to a tuple of ``(parent, Footprint)`` (`dict`).
    """

    log = None
    """Logger used for status messages.
    """

    def __init__(self, config, exposure, footprints, noiseImage=None, exposureId=None, log=None):
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
            if maskname in mask.getMaskPlaneDict():
                # does it already exist?
                plane = mask.getMaskPlane(maskname)
                if self.log:
                    self.log.debug('Mask plane "%s" already existed', maskname)
            else:
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
        # NOTE: heavy footprints get destroyed by the transform process in forcedPhotCcd.py
        # or forcedPhotCoadd.py so they are never available for forced measurements.

        # Create in the dict heavies = {id:heavyfootprint}
        for id, fp in footprints.items():
            if fp[1].isHeavy():
                self.heavies[id] = fp[1]
            elif fp[0] == 0:
                self.heavies[id] = afwDet.makeHeavyFootprint(fp[1], mi)

        # ## FIXME: the heavy footprint includes the mask
        # ## and variance planes, which we shouldn't need
        # ## (I don't think we ever want to modify them in
        # ## the input image).  Copying them around is
        # ## wasteful.

        # We now create a noise HeavyFootprint for each source with has a heavy footprint.
        # We'll put the noise footprints in a dict heavyNoise = {id:heavyNoiseFootprint}
        self.heavyNoise = {}
        noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar, exposureId=exposureId)
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
            fp.spans.setMask(mask, self.otherbitmask)

    def insertSource(self, id):
        """Insert the heavy footprint of a given source into the exposure.

        Parameters
        ----------
        id : `int`
            ID of the source to insert from original dictionary of footprints.

        Notes
        -----
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
        fp.spans.setMask(mask, self.thisbitmask)
        fp.spans.clearMask(mask, self.otherbitmask)

    def removeSource(self, id):
        """Replace the heavy footprint of a given source with noise.

        The same artificial noise is used as in the original replacement.

        Parameters
        ----------
        id : `int`
            ID of the source to replace from original dictionary of footprints.

        Notes
        -----
        Also restores the mask plane.
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
        fp.spans.clearMask(mask, self.thisbitmask)
        fp.spans.setMask(mask, self.otherbitmask)

    def end(self):
        """End the NoiseReplacer.

        Restores original data to the exposure from the heavies dictionary and
        the mask planes to their original state.
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
        """Return a generator of artificial noise.

        Returns
        -------
        noiseGenerator : `lsst.afw.image.noiseReplacer.NoiseGenerator`
        """
        if noiseImage is not None:
            return ImageNoiseGenerator(noiseImage)
        rand = None
        if self.noiseSeedMultiplier:
            # default plugin, our seed
            if exposureId is not None and exposureId != 0:
                seed = exposureId*self.noiseSeedMultiplier
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
            except Exception:
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
            except Exception:
                if self.log:
                    self.log.debug('Failed to get BGMEAN from exposure metadata')

        if noiseSource == 'variance':
            if self.log:
                self.log.debug('Will draw noise according to the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset, rand=rand)

        if noiseSource == 'variance_median':
            if self.log:
                self.log.debug('Will draw noise using the median of the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            s = afwMath.makeStatistics(var, afwMath.MEDIAN)
            varMedian = s.getValue(afwMath.MEDIAN)
            if self.log:
                self.log.debug("Measured from variance: median variance = %g",
                               varMedian)
            return FixedGaussianNoiseGenerator(offset, math.sqrt(varMedian), rand=rand)

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
    """Make a list of NoiseReplacers behave like a single one.

    This class provides conenient syntactic sugar for noise replacement across
    multple exposures.

    Notes
    -----
    This is only used in the MultiFit driver, but the logic there is already
    pretty complex, so it's nice to have this to simplify it.
    """

    def __init__(self, exposuresById, footprintsByExp):
        # exposuresById --- dict of {exposureId: exposure} (possibly subimages)
        # footprintsByExp --- nested dict of {exposureId: {objId: (parent, footprint)}}
        list.__init__(self)
        for expId, exposure in exposuresById.items():
            self.append(NoiseReplacer(exposure, footprintsByExp[expId]), expId)

    def insertSource(self, id):
        """Insert original pixels of the given source (by id) into the exposure.
        """
        for item in self:
            self.insertSource(id)

    def removeSource(self, id):
        """Insert noise pixels of the given source (by id) into the exposure.
        """
        for item in self:
            self.removeSource(id)

    def end(self):
        """Clean-up when the use of the noise replacer is done.
        """
        for item in self:
            self.end()


class NoiseGenerator:
    r"""Base class for noise generators.

    Derived classes produce
    `~lsst.afw.detection.heavyFootprint.HeavyFootprint`\ s filled with noise
    generated in various ways.

    Notes
    -----
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
    """Generate noise by extracting a sub-image from a user-supplied image.

    Parameters
    ----------
    img : `lsst.afw.image.ImageF`
        An image to use as the basis of noise generation.
    """

    def __init__(self, img):
        self.mim = afwImage.MaskedImageF(img)
        self.mean = afwMath.makeStatistics(img, afwMath.MEAN)
        self.std = afwMath.makeStatistics(img, afwMath.STDEV)

    def getMaskedImage(self, bb):
        return self.mim


class GaussianNoiseGenerator(NoiseGenerator):
    """Abstract base for Gaussian noise generators.
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
    """Generates Gaussian noise with a fixed mean and standard deviation.
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
    """Generates Gaussian noise with variance matching an image variance plane.

    Parameters
    ----------
    var : `lsst.afw.image.ImageF`
        The input variance image.
    mean : `float` or `lsst.afw.image.Image`, optional.
        Mean value for the generated noise.
    """

    def __init__(self, var, mean=None, rand=None):
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


class DummyNoiseReplacer:
    """A noise replacer which does nothing.

    This is used when we need to disable noise replacement.

    Notes
    -----
    This has all the public methods of `NoiseReplacer`, but none of them do
    anything.
    """

    def insertSource(self, id):
        pass

    def removeSource(self, id):
        pass

    def end(self):
        pass
