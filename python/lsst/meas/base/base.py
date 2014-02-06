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
""" Base class for plugin registries.
    The Plugin class allowed in the registry is defined on the ctor of the registry
    The intention is that single frame and multi frame algorithms will have different
    registries.
"""
import collections
import math
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pipe.base
import lsst.pex.config
import lsst.pex.exceptions
import lsst.meas.algorithms
import lsst.afw.table

from .baseLib import *

FATAL_EXCEPTIONS = (MemoryError,)  # Exceptions that the framework should always propagate up

def generateAlgorithmName(AlgClass):
    """Generate a string name for an algorithm class that strips away terms that are generally redundant
    while (hopefully) remaining easy to trace to the code.
    """
    name = AlgClass.__name__
    pkg = AlgClass.__module__
    name.strip("Algorithm")
    terms = pkg.split(".")
    if terms[-1].endswith("Lib"):
        terms = terms[:-1]
    if terms[0] == "lsst":
        terms = terms[1:]
    if terms[0] == "meas":
        terms = terms[1:]
    if name.lower().startswith(terms[-1].lower()):
        terms = terms[:-1]
    return "%s_%s" % ("_".join(terms), name)

def callMeasure(task, measRecord, *args, **kwds):
    """Call the measure() method on all plugins in the given task, handling exceptions in a consistent way.

    If all measurement tasks had a common base class, this would probably go there.
    """
    for plugin in task.plugins.iter():
        try:
            plugin.measure(measRecord, *args, **kwds)
        except FATAL_EXCEPTIONS:
            raise
        except lsst.pex.exceptions.LsstCppException as cppError:
            if isinstance(cppError.args[0], MeasurementError):
                error = cppError.args[0]
            else:
                task.log.warn("Error in %s.measure on record %s: %s"
                              % (plugin.name, measRecord.getId(), cppErr))
                error = None
            plugin.fail(measRecord, error)
        except Exception as error:
            task.log.warn("Error in %s.measure on record %s: %s"
                          % (plugin.name, measRecord.getId(), error))
            plugin.fail(measRecord)

def callMeasureN(task, measCat, *args, **kwds):
    """Call the measureN() method on all plugins in the given task, handling exceptions in a consistent way.

    If all measurement tasks had a common base class, this would probably go there.
    """
    for plugin in task.plugins.iterN():
        try:
            plugin.measureN(measCat, *args, **kwds)
        except FATAL_EXCEPTIONS:
            raise
        except lsst.pex.exceptions.LsstCppException as cppError:
            if isinstance(cppError.args[0], MeasurementError):
                error = cppError.args[0]
            else:
                task.log.warn("Error in %s.measureN on records %s-%s: %s"
                              % (plugin.name, measCat[0].getId(), measCat[-1].getId(), cppError))
                error = None
            for measRecord in measCat:
                plugin.fail(measRecord, error)
        except Exception as error:
            for measRecord in measCat:
                plugin.fail(measRecord)
            task.log.warn("Error in %s.measureN on records %s-%s: %s"
                          % (plugin.name, measCat[0].getId(), measCat[-1].getId(), error))

class PluginRegistry(lsst.pex.config.Registry):
    """ Base class for plugin registries

    The Plugin class allowed in the registry is defined on the ctor of the registry
    The intention is that single-frame and multi-frame plugins will have different
    registries.
    """

    class Configurable(object):
        """Class used as the actual element in the registry

        Rather than constructing a Plugin instance, it returns a tuple
        of (runlevel, name, config, PluginClass), which can then
        be sorted before the plugins are instantiated.
        """

        __slots__ = "PluginClass", "name"

        def __init__(self, name, PluginClass):
            """ Initialize registry with Plugin Class
            """
            self.name = name
            self.PluginClass = PluginClass

        @property
        def ConfigClass(self): return self.PluginClass.ConfigClass

        def __call__(self, config):
            return (config.executionOrder, self.name, config, self.PluginClass)

    def register(self, name, PluginClass):
        """Register a Plugin class with the given name.

        The same Plugin may be registered multiple times with different names; this can
        be useful if we often want to run it multiple times with different configuration.

        The name will be used as a prefix for all fields produced by the Plugin, and it
        should generally contain the name of the Plugin or Algorithm class itself
        as well as enough of the namespace to make it clear where to find the code.
        """
        lsst.pex.config.Registry.register(self, name, self.Configurable(name, PluginClass))

    def makeField(self, doc, default=None, optional=False, multi=False):
        return lsst.pex.config.RegistryField(doc, self, default, optional, multi)

class PluginMap(collections.OrderedDict):
    """ Map of plugins to be run for a task

    We assume Plugins are added to the PluginMap according to their executionOrder, so this
    class doesn't actually do any of the sorting (though it does have to maintain that order).

    Should later be implemented as Swigged C++ class so we can pass it to C++-implemented Plugins
    in the future.
    """

    def iter(self):
        """Call each plugin in the map which has a measure() method
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasure:
                yield plugin

    def iterN(self):
        """Call each plugin in the map which has a measureN() method
        """
        for plugin in self.itervalues():
            if plugin.config.doMeasureN:
                yield plugin

class BasePluginConfig(lsst.pex.config.Config):
    """Base class for config which should be defined for each measurement plugin.

    Most derived classes will want to override setDefaults() in order to customize
    the default exceutionOrder.

    A derived class whose correspoding Plugin class implements measureN() should
    additionally add a bool doMeasureN field to replace the bool class attribute
    defined here.
    """

    executionOrder = lsst.pex.config.Field(
        dtype=float, default=2.0,
        doc="""Sets the relative order of plugins (smaller numbers run first).

In general, the following values should be used (intermediate values
are also allowed, but should be avoided unless they are needed):
   0.0 ------ centroids and other algorithms that require only a Footprint and
              its Peaks as input
   1.0 ------ shape measurements and other algorithms that require
              getCentroid() to return a good centroid in addition to a
              Footprint and its Peaks.
   2.0 ------ flux algorithms that require both getShape() and getCentroid()
              in addition to the Footprint and its Peaks
   3.0 ------ Corrections applied to fluxes (i.e. aperture corrections, tying
              model to PSF fluxes). All flux measurements should have an
              executionOrder < 3.0, while all algorithms that rely on corrected
              fluxes (i.e. classification) should have executionOrder > 3.0.
"""
        )


    doMeasure = lsst.pex.config.Field(dtype=bool, default=True,
                                      doc="whether to run this plugin in single-object mode")

    doMeasureN = False  # replace this class attribute with a Field if measureN-capable

class BasePlugin(object):
    """Base class for measurement plugins."""

    def fail(self, measRecord, error=None):
        """Record a failure of the measure or measureN() method.

        When measure() raises an exception, the measurement framework
        will call fail() to allow the plugin to set its failure flag
        field(s).  When measureN() raises an exception, fail() will be
        called repeatedly with all the records that were being
        measured.

        If the exception is a MeasurementError, it will be passed as
        the error argument; in all other cases the error argument will
        be None, and the failure will be logged by the measurement
        framework unless it has been explicitly squashed in config.
        """
        raise NotImplementedError("This algorithm thinks it cannot fail; please report this as a bug.")

class MeasurementDataFlags(object):
    """Flags that describe data to be measured, allowing plugins with the same signature but
    different requirements to assert their appropriateness for the data they will be run on.

    This should later be implemented as a Swigged C++ enum
    """

    PRECONVOLVED = 0x01  # the image has already been convolved with its PSF;
                         # set for preconvolved difference images and Kaiser-style coadds

    DIFFERENCE = 0x02    # the image is a difference, and hence may have negative surface brightnesses

    COADD = 0x04  # the image is a coadd

    NO_PSF = 0x08 # the image has no Psf object attached

    NO_WCS = 0x10 # the image has no Wcs object attached

class SourceSlotConfig(lsst.pex.config.Config):
    """Slot configuration which assigns a particular named plugin to each of a set of
    slots.  Each slot allows a type of measurement to be fetched from the SourceTable
    without knowing which algorithm was used to produced the data.  For example, getCentroid()
    and setCentroid() can be used if the centroid slot is setup, without know that the correct
    key is schema.find("centroid.sdss").key

    NOTE: default for each slot must be registered, even if the default is not used.
    """

    centroid = lsst.pex.config.Field(dtype=str, default="centroid.sdss", optional=True,
                                     doc="the name of the centroiding algorithm used to set source x,y")
    shape = lsst.pex.config.Field(dtype=str, default="shape.sdss", optional=True,
                                  doc="the name of the algorithm used to set source moments parameters")
    apFlux = lsst.pex.config.Field(dtype=str, default="flux.sinc", optional=True,
                                   doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = lsst.pex.config.Field(dtype=str, default="flux.gaussian", optional=True,
                                      doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = lsst.pex.config.Field(dtype=str, default="flux.naive", optional=True,
                                    doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = lsst.pex.config.Field(dtype=str, default="flux.gaussian", optional=True,
                                     doc="the name of the algorithm used to set the source inst flux slot")

    def setupTable(self, table, prefix=None):
        """Convenience method to setup a table's slots according to the config definition.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        if prefix is None: prefix = ""
        if self.centroid is not None: table.defineCentroid(prefix + self.centroid)
        if self.shape is not None: table.defineShape(prefix + self.shape)
        if self.apFlux is not None: table.defineApFlux(prefix + self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(prefix + self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(prefix + self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(prefix + self.instFlux)


class BaseMeasurementConfig(lsst.pex.config.Config):
    """Baseconfig class for all measurement driver tasks."""

    prefix = None

    slots = lsst.pex.config.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source."
        )

    doReplaceWithNoise = lsst.pex.config.Field(dtype=bool, default=True, optional=False,
        doc='When measuring, replace other detected footprints with noise?')

    #  These parameters are to allow the internal NoiseReplacer to be parameterized.
    noiseSource = lsst.pex.config.ChoiceField(
        doc='How to choose mean and variance of the Gaussian noise we generate?',
        dtype=str, allowed={'measure': 'Measure clipped mean and variance from the whole image',
        'meta': 'Mean = 0, variance = the "BGMEAN" metadata entry',
        'variance': "Mean = 0, variance = the image's variance",},
        default='measure', optional=False)

    noiseOffset = lsst.pex.config.Field(dtype=float, optional=False, default=0.,
                                  doc='Add ann offset to the generated noise.')

    noiseSeed = lsst.pex.config.Field(dtype=int, default=0,
        doc='The seed value to use for random number generation.')

    def validate(self):
        lsst.pex.config.Config.validate(self)
        if self.slots.centroid is not None and self.slots.centroid not in self.plugins.names:
            raise ValueError("source centroid slot algorithm is not being run.")
        if self.slots.shape is not None and self.slots.shape not in self.plugins.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None:
                for name in self.plugins.names:
                    if len(name) <= len(slot) and name == slot[:len(name)]:
                        break
                else:
                    raise ValueError("source flux slot algorithm '%s' is not being run." % slot)

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

            NOTE: as the code currently stands, the heavy footprint for a deblended object must be available
            from the input catalog.  If it is not, it cannot be reproduced here.  In that case, the
            topmost parent in the objects parent chain must be used.  The heavy footprint for that source
            is created in this class from the masked image.
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
        # and just use them for children which do not already have heavy footprints.
        # If a heavy footprint is available for a child, we will use it. Otherwise,
        # we use the first parent in the parent chain which has a heavy footprint,
        # which with the one level deblender will alway be the topmost parent
        # NOTE: heavy footprints get destroyed by the transform process in forcedImage.py,
        # so they are never available for forced measurements.

        # Create in the dict heavies = {id:heavyfootprint}
        for id in footprints.keys():
            fp = footprints[id]
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
        # usedid can point either to this source, or to the first parent in the
        # parent chain which has a heavy footprint (or to the topmost parent,
        # which always has one)
        usedid = id
        while self.footprints[usedid][0] != 0 and not usedid in self.heavies.keys():
            usedid = self.footprints[usedid][0]
        fp = self.heavies[usedid]
        fp.insert(im)
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

        # use the same algorithm as in remove Source to find the heavy noise footprint
        # which will undo what insertSource(id) does
        usedid = id
        while self.footprints[usedid][0] != 0 and not usedid in self.heavies.keys():
            usedid = self.footprints[usedid][0]
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
            # default plugin, our seed
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
        """Insert the original pixels for a given source (by id) into the original exposure.
        """
        for item in self: self.insertSource(id)

    def removeSource(self, id):
        """Insert the noise pixels for a given source (by id) into the original exposure.
        """
        for item in self: self.removeSource(id)

    def end(self):
        """Cleanup when the use of the Noise replacer is done.
        """
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

