import collections

import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base
import lsst.pex.config
import pdb

# We'll need to add a few methods to SourceCatalog objects to support iteration over deblend families.
#  - Catalog.getChildren(parent, *args)
#    Return the subset of the catalog where the parent field equals the given value; additional args
#    are additional catalogs that sould be subset using the same indices.  Because SourceCatalogs
#    should generally be sorted by parent automatically, this should be able to operate via slice
#    indexing.
import lsst.afw.table

class AlgorithmRegistry(lsst.pex.config.Registry):

    class Configurable(object):
        """Class used as the actual element in the registry; rather
        than constructing an Algorithm instance, it returns a tuple
        of (runlevel, name, config, AlgorithmClass), which can then
        be sorted before the algorithms are instantiated.
        """

        __slots__ = "AlgorithmClass", "name"

        def __init__(self, name, AlgorithmClass):
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

# We'd probably actually implement this as Swigged C++ class based on std::map; this is just
# a Python placeholder to demo the API.
class AlgorithmMap(collections.OrderedDict):

    def iterSingle(self):
        for algorithm in self.itervalues():
            if algorithm.doMeasureSingle:
                yield algorithm

    def iterMulti(self):
        for algorithm in self.itervalues():
            if algorithm.doMeasureMulti:
                yield algorithm

class BaseAlgorithmConfig(lsst.pex.config.Config):
    executionOrder = lsst.pex.config.Field(dtype=float, default=1.0, doc="sets relative order of algorithms")
    doMeasureSingle = lsst.pex.config.Field(dtype=bool, default=True,
                                            doc="whether to run this algorithm in single-object mode")
    doMeasureMulti = False  # replace this class attribute with a Field if measureMulti-capable

class BaseAlgorithm(object):
    pass

# This should be a Swigged C++ enum; this is just a placeholder.
class MeasurementDataFlags(object):
    """Flags that describe data to be measured, allowing algorithms with the same signature but
    different requirements to assert their appropriateness for the data they will be run on.
    """

    PRECONVOLVED = 0x01  # the image has already been convolved with its PSF;
                         # set for preconvolved difference images and Kaiser-style coadds

    DIFFERENCE = 0x02    # the image is a difference, and hence may have negative surface brightnesses

    COADD = 0x04  # the image is a coadd

    NO_PSF = 0x08 # the image has no Psf object attached

    NO_WCS = 0x10 # the image has no Wcs object attached

class NoiseReplacer(object):
    """Class that handles replacing sources with noise during measurement.

    This will contain the code currently in ReplaceWithNoiseTask, but with a slightly different API
    needed for the new measurement tasks.
    """

    def __init__(self, exposure, footprints, noiseSource, noiseOffset, noiseSeed):
        # 'footprints' is a dict of {id: (parent, footprint)}; when used in SFM, the ID will be the
        # source ID, but in forced photometry, this will be the reference ID, as that's what we used to
        # determine the deblend families.  This routine should create HeavyFootprints for any non-Heavy
        # Footprints, and replace them in the dict.  It should then create a dict of HeavyFootprints
        # containing noise, but only for parent objects, then replace all sources with noise.  This is
        # mostly like the current ReplaceWithNoiseTask.begin(), but with a different API.
        # This should ignore any footprints that lay outside the bounding box of the exposure, and clip
        # those that lie on the border.
        noiseImage=None
        noiseMeanVar=None
        self.noiseSource = noiseSource
        self.noiseOffset = noiseOffset
        self.noiseSeed = noiseSeed
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
                print 'Mask plane "%s" already existed' % maskname
            except:
                # if not, add it; we should delete it when done.
                plane = mask.addMaskPlane(maskname)
                self.removeplanes.append(maskname)
            mask.clearMaskPlane(plane)
            bitmask = mask.getPlaneBitMask(maskname)
            bitmasks.append(bitmask)
            print 'Mask plane "%s": plane %i, bitmask %i = 0x%x' % (maskname, plane, bitmask, bitmask)
        self.thisbitmask,self.otherbitmask = bitmasks
        del bitmasks

        # Start by creating HeavyFootprints for each source.
        #
        # The "getParent()" checks are here because top-level
        # sources (ie, those with no parents) are not supposed to
        # have HeavyFootprints, but child sources (ie, those that
        # have been deblended) should have HeavyFootprints
        # already.
        self.heavies = {}
        for id in footprints.keys():
            fp = footprints[id][1]
            ### FIXME: the heavy footprint includes the mask
            ### and variance planes, which we shouldn't need
            ### (I don't think we ever want to modify them in
            ### the input image).  Copying them around is
            ### wasteful.
            heavy = afwDet.makeHeavyFootprint(fp, mi)
            self.heavies[id] = heavy

        # We now create a noise HeavyFootprint for each top-level Source.
        # We'll put the noisy footprints in a map from id -> HeavyFootprint:
        self.heavyNoise = {}
        noisegen = self.getNoiseGenerator(exposure, noiseImage, noiseMeanVar)
        print 'Using noise generator: %s' % (str(noisegen))
        for id in footprints.keys():
            fp = footprints[id][1]
            parent = footprints[id][0]
            if parent:
                continue
            heavy = noisegen.getHeavyFootprint(fp)
            self.heavyNoise[id] = heavy
            # Also insert the noisy footprint into the image now.
            # Notice that we're just inserting it into "im", ie,
            # the Image, not the MaskedImage.
            heavy.insert(im)
            # Also set the OTHERDET bit
            afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def insertSource(self, id):
        # Copy this source's pixels into the image
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        fp = self.heavies[id]
        fp.insert(im)
        afwDet.setMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.clearMaskFromFootprint(mask, fp, self.otherbitmask)

    def removeSource(self, id):
        # remove a single source
        # (Replace this source's pixels by noise again.)
        # Do this by finding the source's top-level ancestor
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        # uses the footprints to find the topmost parent
        topId = id
        while topId:
            parentId = self.footprints[topId][0]
            if parentId == 0:
                break
            topId = parentId

        # Re-insert the noise pixels
        fp = self.heavyNoise[topId]
        fp.insert(im)
        # Clear the THISDET mask plane.
        afwDet.clearMaskFromFootprint(mask, fp, self.thisbitmask)
        afwDet.setMaskFromFootprint(mask, fp, self.otherbitmask)

    def end(self):
        # restores original image, cleans up temporaries
        # (ie, replace all the top-level pixels)
        mi = self.exposure.getMaskedImage()
        im = mi.getImage()
        mask = mi.getMask()
        for id in self.footprints.keys():
            fp = self.footprints[id][1]
            parentId = self.footprints[id][0]
            if parentId:
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
                print 'Using passed-in noise mean = %g, variance = %g -> stdev %g' % (noiseMean, noiseVar, noiseStd)
                return FixedGaussianNoiseGenerator(noiseMean, noiseStd, rand=rand)
            except:
                print 'Failed to cast passed-in noiseMeanVar to floats: %s' % (str(noiseMeanVar))
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
                print 'Using noise variance = (BGMEAN = %g) from exposure metadata' % (bgMean)
                return FixedGaussianNoiseGenerator(offset, noiseStd, rand=rand)
            except:
                print 'Failed to get BGMEAN from exposure metadata'

        if noiseSource == 'variance':
            print 'Will draw noise according to the variance plane.'
            var = exposure.getMaskedImage().getVariance()
            return VariancePlaneNoiseGenerator(var, mean=offset, rand=rand)

        # Compute an image-wide clipped variance.
        im = exposure.getMaskedImage().getImage()
        s = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        noiseMean = s.getValue(afwMath.MEANCLIP)
        noiseStd = s.getValue(afwMath.STDEVCLIP)
        print "Measured from image: clipped mean = %g, stdev = %g" % (noiseMean,noiseStd)
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



class ReplaceWithNoiseConfig(lsst.pex.config.Config):
    """
    Configuration for ReplaceWithNoiseTask.
    """
    noiseSource = lsst.pex.config.ChoiceField(doc='How do we choose the mean and variance of the Gaussian noise we generate?',
                                      dtype=str, allowed={
                                          'measure': 'Measure clipped mean and variance from the whole image',
                                          'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
                                          'variance': "Mean = 0, variance = the image's variance",
                                          },
                                      default='measure',
                                      optional=False)

    noiseOffset = lsst.pex.config.Field(dtype=float, optional=False, default=0.,
                                  doc='Add ann offset to the generated noise.')

    noiseSeed = lsst.pex.config.Field(dtype=int, default=0, doc='The seed value to use for random number generation.')


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

