"""Base classes for forced measurement plugin algorithms and the driver task for these.

In forced measurement, a reference catalog is used to define restricted measurements (usually just fluxes)
on an image.  As the reference catalog may be deeper than the detection limit of the measurement image, we
do not assume that we can use detection and deblend information from the measurement image.  Instead, we
assume this information is present in the reference catalog and has been "transformed" in some sense to
the measurement frame.  At the very least, this means that Footprints from the reference catalog should
be transformed and installed as Footprints in the output measurement catalog.  If we have a procedure that
can transform HeavyFootprints, we can then proceed with measurement as usual, but using the reference
catalog's id and parent fields to define deblend families.  If this transformation does not preserve
HeavyFootprints (this is currently the case), then we will only be able to replace objects with noise
one deblend family at a time, and hence measurements run in single-object mode may be contaminated by
neighbors when run on objects with parent != 0.

Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).  Note that
the reference catalog may be in a different coordinate system; it is the responsibility of algorithms
to transform the data they need themselves, using the reference WCS provided.  However, for algorithms
that only require a position, they may simply use output SourceCatalog's centroid slot, which will generally
be set to the transformed position of the reference object before any algorithms are run, and hence avoid
using the reference catalog at all.
"""

#
# NOTE: to fully implement forced measurement, we'll need not just this measurement task, but also a
# command-line driver task that uses it as a subtask.  A good prototype for this already exists on
# the HSC branch of pipe_tasks, in the forcedPhot*.py and references.py files.  These should be
# transferred to the LSST side of pipe_tasks, and modified to use the new measurement API described
# here.  Note that it will also be responsible for transforming Footprints (in the generateSources())
# method, and attaching them to sources, as this is no longer something measurement tasks do.  There
# are also some hackish workarounds in that code that could be removed with properly vetted
# modifications to pipe_base (i.e. we need to get a Butler in Task ctors).  We should do that now
# as well.
#
# It's also worth considering merging that command-line-driver task with this measurement task, and
# flattening the hierarchy.  In that case, though, we'd probably want to put the base command-line
# task in meas_base and let the CCD- and coadd-specialized subclasses continue to be in pipe_tasks
# to avoid dependency issues.
#

import lsst.pex.config
import lsst.pipe.base
import lsst.daf.base

from .base import *

class ForcedAlgorithmConfig(BaseAlgorithmConfig):
    """Base class for configs of forced measurement plugin algorithms."""
    pass

class ForcedAlgorithm(BaseAlgorithm):

    # All subclasses of ForcedAlgorithm should be registered here
    registry = AlgorithmRegistry(ForcedAlgorithm)

    ConfigClass = ForcedAlgorithmConfig

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

    def measureSingle(self, exposure, measRecord, refRecord, referenceWcs):
        """Measure the properties of a source on a single image, given data from a
        reference record.

        @param[in] exposure       lsst.afw.image.ExposureF, containing the pixel data to
                                  be measured and the associated Psf, Wcs, etc.  All
                                  other sources in the image will have been replaced by
                                  noise according to deblender outputs.
        @param[in,out] measRecord lsst.afw.table.SourceRecord to be filled with outputs,
                                  and from which previously-measured quantities can be
                                  retreived.
        @param[in] refRecord      lsst.afw.table.SimpleRecord that contains additional
                                  parameters to define the fit, as measured elsewhere.
        @param[in] refWcs         The coordinate system for the reference catalog values.
                                  An afw.geom.Angle may be passed, indicating that a
                                  local tangent Wcs should be created for each object
                                  using afw.image.makeLocalWcs and the given angle as
                                  a pixel scale.

        In the normal mode of operation, the source centroid will be set to the
        WCS-transformed position of the reference object, so algorithms that only
        require a reference position should not have to access the reference object
        at all.
        """
        raise NotImplementedError()

    def measureMulti(self, exposure, measCat, refCat, refWcs):
        """Measure the properties of a group of blended sources on a single image,
        given data from a reference record.

        @param[in] exposure      lsst.afw.image.ExposureF, containing the pixel data to
                                 be measured and the associated Psf, Wcs, etc.  Sources
                                 not in the blended hierarchy to be measured will have
                                 been replaced with noise using deblender outputs.
        @param[in,out] measCat   lsst.afw.table.SourceCatalog to be filled with outputs,
                                 and from which previously-measured quantities can be
                                 retrieved, containing only the sources that should be
                                 measured together in this call.
        @param[in] refCat        lsst.afw.table.SimpleCatalog that contains additional
                                 parameters to define the fit, as measured elsewhere.
                                 Ordered such that zip(sources, references) may be used.
        @param[in] refWcs        The coordinate system for the reference catalog values.
                                 An afw.geom.Angle may be passed, indicating that a
                                 local tangent Wcs should be created for each object
                                 using afw.image.makeLocalWcs and the given Angle as
                                 a pixel scale.

        In the normal mode of operation, the source centroids will be set to the
        WCS-transformed position of the reference object, so algorithms that only
        require a reference position should not have to access the reference object
        at all.
        """
        raise NotImplementedError()

class ForcedMeasurementConfig(lsst.pex.config.Config):
    """Config class for forced measurement driver task."""

    algorithms = ForcedAlgorithm.registry.makeField(
        multi=True,
        default=[
            # TODO
            ],
        doc="Plugin algorithms to be run and their configuration"
        )

    references = ConfigurableField(target=ReferencesTask, doc="Retrieve reference source catalog")

class ForcedMeasurementTask(lsst.pipe.base.Task):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclasses for running on CCDs and Coadds).
    """

    ConfigClass = ForcedMeasurementConfig

    def __init__(self, refSchema, flags, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.makeSubTask("references")
        self.schemaMapper = lsst.afw.table.SchemaMapper(refSchema)
        self.algMetadata = lsst.daf.base.PropertyList()
        self.algorithms = AlgorithmMap()
        for executionOrder, name, config, AlgorithmClass in sorted(self.config.algorithms.apply()):
            self.algorithms[name] = AlgorithmClass(config, name, schemaMapper=schemaMapper, flags=flags,
                                                   others=self.algorithms, metadata=self.algMetadata)

    def run(self, exposure, measCat, refCat, refWcs):
        footprints = {refRecord.getId(): (refRecord.getParent(), measRecord.getFootprint())
                      for refRecord, measRecord in zip(refCat, measCat)}
        noiseReplacer = NoiseReplacer(exposure, footprints)
        refParentCat, measParentCat = refCat.getChildren(0, measCat)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(refParentCat, measParentCat):
            refChildCat, measChildCat = refCat.getChildren(refParentRecord.getId(), measCat)
            # TODO: skip this loop if there are no algorithms configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildRecord, measChildRecord):
                noiseReplacer.insertSource(refChildRecord.getId())
                for algorithm in self.algorithms.iterSingle():
                    algorithm.measureSingle(exposure, measChildRecord, refChildRecord)
                noiseReplacer.removeSource(refChildRecord.getId())
            noiseReplacer.insertSource(refParentRecord.getId())
            for algorithm in self.algorithms.iterSingle():
                algorithm.measureSingle(exposure, measParentRecord, refParentRecord)
            for algorithm in self.algorithms.iterMulti():
                algorithm.measureMulti(exposure, measChildCat, refChildCat)
                algorithm.measureMulti(exposure, measParentCat[parentIdx:parentIdx+1],
                                       refParentCat[parentIdx:parentIdx+1])
            noiseReplacer.removeSource(refChildRecord.getId())
        noiseReplacer.end()

