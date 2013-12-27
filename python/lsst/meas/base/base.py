import collections

import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base
import lsst.pex.config

import lsst.afw.table

""" Base class for algorithm registries.
    The Algorithm class allowed in the registry is defined on the ctor of the registry
    The intention is that single frame and multi frame algorithms will have different
    registries.
"""
class AlgorithmRegistry(lsst.pex.config.Registry):

    class Configurable(object):
        """ Class used as the actual element in the registry; rather
        than constructing an Algorithm instance, it returns a tuple
        of (runlevel, name, config, AlgorithmClass), which can then
        be sorted before the algorithms are instantiated.
        """

        __slots__ = "AlgorithmClass", "name"
        def __init__(self, name, AlgorithmClass):
            """ Initialize registry with Algorithm Class
           """
            self.name = name
            self.AlgorithmClass = AlgorithmClass

        @property
        def ConfigClass(self): return self.AlgorithmClass.ConfigClass

        def __call__(self, config):
            return (config.executionOrder, self.name, config, self.AlgorithmClass)

    def register(self, name, AlgorithmClass):
        lsst.pex.config.Registry.register(self, name, self.Configurable(name, AlgorithmClass))

    def makeField(self, doc, default=None, optional=False, multi=False):
        return lsst.pex.config.RegistryField(doc, self, default, optional, multi)

class AlgorithmMap(collections.OrderedDict):
    """ Map of algorithms to be run for a task 
        Should later be implemented as Swigged C++ class based on std::map
    """

    def iterSingle(self):
        """ Call each algorithm in the map which has a measureSingle """
        for algorithm in self.itervalues():
            if algorithm.config.doMeasureSingle:
                yield algorithm

    def iterMulti(self):
        """ Call each algorithm in the map which has a measureMulti """
        for algorithm in self.itervalues():
            if algorithm.config.doMeasureMulti:
                yield algorithm

class BaseAlgorithmConfig(lsst.pex.config.Config):
    """Base class for config which should be defined for each measurement algorithm."""

    executionOrder = lsst.pex.config.Field(dtype=float, default=1.0, doc="sets relative order of algorithms")
    doMeasureSingle = lsst.pex.config.Field(dtype=bool, default=True,
                      doc="whether to run this algorithm in single-object mode")
    doMeasureMulti = False  # replace this class attribute with a Field if measureMulti-capable

class BaseAlgorithm(object):
    """Base class for measurement algorithms."""
    pass

class MeasurementDataFlags(object):
    """Flags that describe data to be measured, allowing algorithms with the same signature but
    different requirements to assert their appropriateness for the data they will be run on.
    This should later be implemented as a Swigged C++ enum
    """

    PRECONVOLVED = 0x01  # the image has already been convolved with its PSF;
                         # set for preconvolved difference images and Kaiser-style coadds

    DIFFERENCE = 0x02    # the image is a difference, and hence may have negative surface brightnesses

    COADD = 0x04  # the image is a coadd

    NO_PSF = 0x08 # the image has no Psf object attached

    NO_WCS = 0x10 # the image has no Wcs object attached

class NoiseReplacer(object):
    """ Class that handles replacing sources with noise during measurement.
        This is a functional copy of the code in  ReplaceWithNoiseTask, but with a slightly different API
        needed for the new measurement framework.
    """

    def __init__(self, exposure, footprints, noiseSource, noiseOffset, noiseSeed, log=None):
        """ Initialize the Noise Replacer.
            @param[in,out]  exposure     Exposure to be noise replaced. (All sources replaced on output)
            @param[in]      footprints   dict of {id: (parent, footprint)};
            @param[in]      noiseSource  Source of Gaussian noise fill
            @param[in]      noiseOffset  Offset of the mean of Gaussian noise fill
            @param[in]      noiseSeed    Random number generator seed.

            'footprints' is a dict of {id: (parent, footprint)}; when used in SFM, the ID will be the
            source ID, but in forced photometry, this will be the reference ID, as that's what we used to
            determine the deblend families.  This routine should create HeavyFootprints for any non-Heavy
            Footprints, and replace them in the dict.  It should then create a dict of HeavyFootprints
            containing noise, but only for parent objects, then replace all sources with noise.
            This should ignore any footprints that lay outside the bounding box of the exposure,
            and clip those that lie on the border.
    """
        noiseImage=None
        noiseMeanVar=None
        self.noiseSource = noiseSource
        self.noiseOffset = noiseOffset
        self.noiseSeed = noiseSeed
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
                if self.log: self.log.logdebug('Mask plane "%s" already existed' % maskname)
            except:
                # if not, add it; we should delete it when done.
                plane = mask.addMaskPlane(maskname)
                self.removeplanes.append(maskname)
            mask.clearMaskPlane(plane)
            bitmask = mask.getPlaneBitMask(maskname)
            bitmasks.append(bitmask)
            if self.log: self.log.logdebug('Mask plane "%s": plane %i, bitmask %i = 0x%x'
                % (maskname, plane, bitmask, bitmask))
        self.thisbitmask,self.otherbitmask = bitmasks
        del bitmasks
        self.heavies = {}
        # Start by creating HeavyFootprints for each source which has no parent
        # and just use them for children which already have heavy footprints.
        # Store in the dict heavies = {id:heavyfootprint}
        for id in footprints.keys():
            fp = footprints[id]
            if fp[1].isHeavy():
                self.heavies[id] = afwDet.castHeavyFootprintF(fp[1])
            elif fp[0] == 0:
                self.heavies[id] = afwDet.makeHeavyFootprint(fp[1], mi)
        # We now create a noise HeavyFootprint for each top-level Source.
        # We'll put the noisy footprints in a dict heavyNoise = {id:heavyNoiseFootprint}

        # top-level source: copy pixels from the input
        # image.
        ### FIXME: the heavy footprint includes the mask
        ### and variance planes, which we shouldn't need
        ### (I don't think we ever want to modify them in
        ### the input image).  Copying them around is
        ### wasteful.

        self.heavyNoise = {}
        noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar)
        #  The noiseGenMean and Std are used by the unit tests
        self.noiseGenMean = noisegen.mean
        self.noiseGenStd = noisegen.std
        if self.log: self.log.logdebug('Using noise generator: %s' % (str(noisegen)))
        for id in self.heavies.keys():
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
        """ Insert the heavy footprint of a given source into the exposure
            @param[in]  id   id for current source to insert
                             from original footprint dict
            Fix up the mask plane to show the source of this footprint
        """
        # Copy this source's pixels into the image
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        # usedid can point either to this source, or to the topmost parent 
        usedid = id
        while self.footprints[usedid][0] != 0 and not usedid in self.heavies.keys():
            usedid = self.footprints[usedid][0]
        # Perhaps we should assert this case?
        if not usedid in self.heavies.keys():
            print "footprint not found for ", id
            return
        fp = self.heavies[usedid]
        try:
            fp.insert(im)
        except Exception, e:
            print "BAD"
            import pdb
            pdb.set_trace()
        afwDet.setMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.clearMaskFromFootprint(mask, fp, self.otherbitmask)

    def removeSource(self, id):
        """ Remove the heavy footprint of a given source and replace with previous noise
            @param[in]  id   id for current source to insert
                             from original footprint dict
            Fix up the mask plane to show the source of this footprint
        """
        # remove a single source
        # (Replace this source's pixels by noise again.)
        # Do this by finding the source's top-level ancestor
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        # usedid can point either to this source, or to the topmost parent 
        usedid = id
        while self.footprints[usedid][0] != 0 and not usedid in self.heavies.keys():
            usedid = self.footprints[usedid][0]
        # Perhaps we should assert this case?
        if not usedid in self.heavies.keys():
            print "footprint not found for ", id
            return
        # Re-insert the noise pixels
        fp = self.heavyNoise[usedid]
        fp.insert(im)
        # Clear the THISDET mask plane.
        afwDet.clearMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def end(self):
        """" End the NoiseReplacer.
             Restore original data to the exposure from the heavies dictionary
             Restore the mask planes to their original state
        """
        # restores original image, cleans up temporaries
        # (ie, replace all the top-level pixels)
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        for id in self.footprints.keys():
            fp = self.footprints[id][1]
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

    def getNoiseGenerator(self, exposure, noiseImage, noiseMeanVar):
        """" Generate noise image using parameters given
        """
        if noiseImage is not None:
            return ImageNoiseGenerator(noiseImage)
        rand = None
        if self.noiseSeed:
            # default algorithm, our seed
            rand = afwMath.Random(afwMath.Random.MT19937, self.noiseSeed)
        if noiseMeanVar is not None:
            try:
                # Assume noiseMeanVar is an iterable of floats
                noiseMean,noiseVar = noiseMeanVar
                noiseMean = float(noiseMean)
                noiseVar = float(noiseVar)
                noiseStd = math.sqrt(noiseVar)
                if self.log: self.log.logdebug('Using passed-in noise mean = %g, variance = %g -> stdev %g'
                     % (noiseMean, noiseVar, noiseStd))
                return FixedGaussianNoiseGenerator(noiseMean, noiseStd, rand=rand)
            except:
                if self.log: self.log.logdebug('Failed to cast passed-in noiseMeanVar to floats: %s' 
                    % (str(noiseMeanVar)))
        offset = self.noiseOffset
        noiseSource = self.noiseSource

        if noiseSource == 'meta':
            # check the exposure metadata
            meta = exposure.getMetadata()
            # this key name correspond to estimateBackground() in detection.py
            try:
                bgMean = meta.getAsDouble('BGMEAN')
                # We would have to adjust for GAIN if ip_isr didn't make it 1.0
                noiseStd = math.sqrt(bgMean)
                if self.log: self.log.logdebug('Using noise variance = (BGMEAN = %g) from exposure metadata'
                    % (bgMean))
                return FixedGaussianNoiseGenerator(offset, noiseStd, rand=rand)
            except:
                if self.log: self.log.logdebug('Failed to get BGMEAN from exposure metadata')

        if noiseSource == 'variance':
            if self.log: self.log.logdebug('Will draw noise according to the variance plane.')
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset, rand=rand)

        # Compute an image-wide clipped variance.
        im = exposure.getMaskedImage().getImage()
        s = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        noiseMean = s.getValue(afwMath.MEANCLIP)
        noiseStd = s.getValue(afwMath.STDEVCLIP)
        if self.log: self.log.logdebug("Measured from image: clipped mean = %g, stdev = %g"
            % (noiseMean,noiseStd))
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
        for expId, exposure in exposuresById.iteritems():
            self.append(NoiseReplacer(exposure, footprintsByExp[expId]))

    def insertSource(self, id):
        for item in self: self.insertSource(id)

    def removeSource(self, id):
        for item in self: self.removeSource(id)

    def end(self):
        for item in self: self.end()

class NoiseGenerator(object):
    '''
    Base class for noise generators used by the "doReplaceWithNoise" routine:
    these produce HeavyFootprints filled with noise generated in various ways.

    This is an abstract base class.
    '''
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
    '''
        "Generates" noise by cutting out a subimage from a user-supplied noise Image.
    '''
    def __init__(self, img):
        '''
        img: an afwImage.ImageF
        '''
        self.mim = afwImage.MaskedImageF(img)
    def getMaskedImage(self, bb):
        return self.mim

class GaussianNoiseGenerator(NoiseGenerator):
    '''
        Generates noise using the afwMath.Random() and afwMath.randomGaussianImage() routines.

        This is an abstract base class.
    '''
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
    '''
        Generates Gaussian noise with a fixed mean and standard deviation.
    '''
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
    '''
        Generates Gaussian noise whose variance matches that of the variance plane of the image.
    '''
    def __init__(self, var, mean=None, rand=None):
        '''
        var: an afwImage.ImageF; the variance plane.
        mean: floating-point or afwImage.Image
        '''
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

