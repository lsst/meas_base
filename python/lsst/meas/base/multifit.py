"""Base classes for simultaneous multi-epoch measurement ("MultiFit") plugin algorithms and the driver
task for these.

In MultiFit, we assume a reference catalog containing initial measurements is provided (usually from fitting
done on a coadd), and that our job is merely to "tweak up" the result in some way.  As a result, we assume
that all the deblend information is also derived from the reference catalog, and transform reference catalog
Footprints to all measurement frames.  As with forced measurement, if the Footprint transform operation only
supports transforming regular Footprints, and discards the deblend information in HeavyFootprints (as is
currently the case), measurements run in single-object mode may be contaminated by neighbors when run
on objects with parent != 0.

Measurements are recorded in a TAN coordinate system centered very near the object being measured, with
pixel scale set by config.  This is constructed using afw.image.makeLocalWcs(), so it can be reconstructed
using just a position and is aligned with the equatoral celestial coordinate grid (while avoiding problems
with poles and angle-wrapping).  This position is set to the reference catalog position of the parent object
(so all deblend siblings share the same measWcs), and is saved in the output measurement catalog in the
'measOrigin' field.  Note that this position is used not just for the CRVAL of the WCS, but also the
CRPIX - so all position measurements produced by MultiFit algorithms are effectively saved as offsets
from the 'measOrigin'

Similarly, a "measCalib" objects is used to set the flux units for measurements, and
is defined completely by config.
"""

#
# NOTE: like forced measurement, we'll need a command-line driver task in addition to this measurement
# task.  Even more so than with forced measurement, I think it's appropriate for the measurement task
# to actually be the command-line task, and I've made MultiFitTask inherit from CmdLineTask below.
# It will need a few more methods to actually make all that work, but I think they should be pretty
# straightforward (at least once one is familiar with CmdLineTasks in general).
#

import lsst.pipe.base
import lsst.pex.config

from .base import *

class ProjectedSourceData(object):
    # This is the per-epoch data that's passed to each algorithm.
    # It will probably actually have to be a Swigged C++ struct; this is just a placeholder.

    def __init__(self, exposure, footprint, position, fullTransform, localTransform):
        self.exposure = exposure     # Exposure subimage
        self.footprint = footprint   # Transformed footprint
        self.position = position     # Position in Exposure pixel coordinates
        self.fullTransform = fullTransform   # XYTransform from measurement frame to pixel coordinates
        self.localTransform = localTransform # AffineTransform that approximates fullTransform at position

class MultiFitAlgorithmConfig(BaseAlgorithmConfig)

class MultiFitAlgorithm(BaseAlgorithm):

    # All subclasses of MultiFitAlgorithm should be registered here
    registry = AlgorithmRegistry(ForcedAlgorithm)

    ConfigClass = MultiFitAlgorithmConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        """Initialize the measurement object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the algorithm was registered with.
        @param[in,out]  schemaMapper  A SchemaMapper that maps reference catalog fields to output
                                      catalog fields.  Output fields should be added to the
                                      output schema.  While most algorithms will not need to map
                                      fields from the reference schema, if they do so, those fields
                                      will be transferred before any algorithms are run.
        @param[in]  flags        A set of bitflags describing the data that the algorithm
                                 should check to see if it supports.  See MeasuremntDataFlags.
        @param[in]  others       An AlgorithmMap of previously-initialized algorithms
        @param[in]  metadata     Algorithm metadata that will be attached to the output catalog
        """
        self.config = config

    def measureSingle(self, data, measRecord, refRecord, measWcs, refWcs, measCalib, refCalib):
        """Measure the properties of an object using data from multiple epochs.

        @param[in] data            An array of ProjectedSourceData objects, containing the
                                   pixel data to be fit and associated metadata, one for
                                   each epoch.
        @param[in,out] measRecord  lsst.afw.table.SourceRecord to be filled with outputs,
                                   and from which previously-measured quantities can be
                                   retreived.
        @param[in] refRecord       lsst.afw.table.SourceRecord containing previous measurements
                                   to be used to initialize the fit.
        @param[in] measWcs         The coordinate system for output values.
        @param[in] refWcs          The coordinate system for the reference catalog values.
                                   An afw.geom.Angle may be passed, indicating that a
                                   local tangent Wcs should be created for each object
                                   using afw.image.makeLocalWcs and the given angle as
                                   a pixel scale.
        @param[in] measCalib       Photometric system for output values.
        @param[in] refCalib        Photometric system for reference values.
        """
        raise NotImplementedError()

    def measureMulti(self, data, measCat, refCat, measWcs, refWcs, measCalib, refCalib):
        """Measure the properties of a group of blended objects using data from multiple
        epochs.

        @param[in] data            A nested sequence of ProjectedSourceData objects, containing
                                   the pixel data to be fit and associated metadata;
                                   iteration is first over objects and then over epochs, but the
                                   .exposure attributes refer to the same postage stamp
                                   for all objects in each inner sequence.
       @param[in,out] measCat      lsst.afw.table.SourceCatalog to be filled with outputs,
                                   and from which previously-measured quantities can be
                                   retreived.
        @param[in] refCat          lsst.afw.table.SourceCatalog containing previous measurements
                                   to be used to initialize the fit.
        @param[in] measWcs         The coordinate system for output values.
        @param[in] refWcs          The coordinate system for the reference catalog values.
                                   An afw.geom.Angle may be passed, indicating that a
                                   local tangent Wcs should be created for each object
                                   using afw.image.makeLocalWcs and the given angle as
                                   a pixel scale.
        @param[in] measCalib       Photometric system for output values.
        @param[in] refCalib        Photometric system for reference values.

        Arguments are ordered such that zip(data, measCat, refCat) may be used to iterate
        over all sources, yielding arguments that could be passed to measureSingle.
        """
        raise NotImplementedError()

class MultiFitConfig(lsst.pex.config.Config):
    """Config class for simultaneous multi-epoch measurement driver task."""

    algorithms = ForcedAlgorithm.registry.makeField(
        multi=True,
        default=[
            # TODO
            ],
        doc="Plugin algorithms to be run and their configuration"
        )

    measPixelScale = lsst.pex.config.Field(dtype=float, doc="Pixel scale for measurement WCS")

    measFluxMag0 = lsst.pex.config.Field(
        dtype=float,
        doc="Flux at magnitude zero that defines measurement photometric system"
        )

class MultiFitTask(lsst.pipe.base.CmdLineTask):
    """Simultaneous multi-epoch measurement driver task"""

    ConfigClass = MultiFitMeasurementConfig

    def __init__(self, refSchema, flags, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schemaMapper = lsst.afw.table.SchemaMapper(refSchema)
        self.measCalib = lsst.afw.image.Calib(self.config.measFluxMag0)
        self.measOriginKey = self.schemaMapper.addOutputField(
            lsst.afw.table.Field[lsst.afw.coord.Coord](
                "measOrigin",
                ("center of TAN projection used for measurement outputs; "
                 "reconstruct with afw.image.makeLocalWcs()")
                )
            )

    def transformFootprints(self, expCat, refCat, refWcs):
        # Loop over all ExposureRecords in expCat, and transform all Footprints in
        # refCat to each coordinate system.  Return two nested dicts:
        # {refId: {expId: (parent, footprint)}}
        # {expId: {refId: (parent, footprint)}}
        # When a refRecord's transformed Footprint doesn't overlap with an exposure, that
        # (refId, expId) pair should not be included in either list.  If a transformed
        # Footprint goes off the edge of an exposure, it should be clipped.
        # We assume all this fits in memory, because Footprints are highly compressed
        # compared to the pixels, and we can do something fancier someday in the future
        # if necessary.
        raise NotImplementedError()

    def loadExposures(self, butler, footprints):
        # Given a dict of {expId: (parent, footprint)}, load Exposure data from disk and
        # return a dict of {expId: exposure}, where each Exposure contains only the region
        # inside the Footprint bbox (plus some configurable margin)
        raise NotImplementedError()

    def buildAlgorithmInputs(self, measWcs, exposures, footprints, coord):
        # 'data' will actually have to be Swigged C++ container so we can pass it to C++ functions
        data = []
        for expId, (parent, footprint) in footprints.iteritems():
            exposure = exposures[expId],
            position = exposure.getWcs().skyToPixel(coord)
            fullTransform = lsst.afw.image.XYTransformFromWcsPair(exposure.getWcs(), measWcs)
            data.append(
                ProjectedSourceData(
                    exposure=exposure, footprint=footprint, position=position,fullTransform=fullTransform,
                    localTransform=fullTransform.linearizeForwardTransform(measPosition)
                    )
                )
        return measWcs, data

    def run(self, butler, expCat, measCat, refCat, refWcs, refCalib):
        footprintsByRef, footprintsByExp = self.transformFootprints(expCat, refCat, refWcs)

        # Loop over all deblend parents.  We load Exposure subimages from disk at each iteration,
        # and use them for the parent object and its children.
        refParentCat, measParentCat = refCat.getChildren(0, measCat)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(refParentCat, measParentCat):
            exposures = self.loadExposures(butler, footprintsByRef[refParentRecord.getId()])
            noiseReplacers = NoiserReplacerList(exposures, footprintsByExp)

            # As noted in the module docstring, all measurements should be output in this coordinate
            # system, defined by the reference position of the parent.
            measWcs = lsst.afw.image.makeLocalWcs(refParentRecord.getCoord(), self.config.measPixelScale)
            measParentRecord.setCoord(self.measOriginKey, refParentRecord.getCoord())

            # Package up some args to pass to the Algorithm measure* methods; none of these depend
            # on the object itself (they depend at most on its parent)
            algArgs = (measWcs, refWcs, self.measCalib, refCalib)

            # Loop over children and measure them in single-object mode
            # TODO: skip this loop if there are no algorithms configured for single-object mode
            refChildCat, measChildCat = refCat.getChildren(refParentRecord.getId(), measCat)
            childrenData = []  # TODO: probably has to be a swigged C++ container
            for refChildRecord, measChildRecord in zip(refChildRecord, measChildRecord):
                measChildRecord.setCoord(self.measOriginKey, refParentRecord.getCoord())
                childData = self.buildAlgorithmInputs(
                    measWcs,
                    exposures,
                    footprintsByRef[refChildRecord.getId()],
                    refChildRecord.getCoord()
                    )
                childrenData.append(childData)
                noiseReplacers.insertSource(refChildRecord.getId())
                for algorithm in self.algorithms.iterSingle():
                    algorithm.measureSingle(childData, measChildRecord, refChildRecord, *algArgs)
                noiseReplacers.removeSource(refChildRecord.getId())

            # Measure children in multi-object mode, and parent in both single-object and multi-object modes
            noiseReplacers.insertSource(refParentRecord.getId())
            parentData = self.buildAlgorithmInputs(
                measWcs,
                exposures,
                footprintsByRef[refChildRecord.getId()],
                refChildRecord.getCoord()
                )
            for algorithm in self.algorithms.iterSingle():
                algorithm.measureSingle(parentData, measParentRecord, refParentRecord, *algArgs)
            for algorithm in self.algorithms.iterMulti():
                algorithm.measureMulti(childrenData, measChildCat, refChildCat, *algArgs)
                algorithm.measureMulti([parentData], measParentCat[parentIdx:parentIdx+1],
                                       refParentCat[parentIdx:parentIdx+1])
            noiseReplacers.removeSource(refChildRecord.getId())

