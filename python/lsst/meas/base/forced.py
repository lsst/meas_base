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
from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
import lsst.daf.base
from lsst.pex.config import DictField,ConfigurableField
from lsst.pipe.tasks.references import CoaddSrcReferencesTask
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer
from .base import *
__all__ = ("ForcedAlgorithmConfig", "ForcedAlgorithm", "ForcedMeasurementConfig", "ForcedMeasurementTask")
class ForcedAlgorithmConfig(BaseAlgorithmConfig):
    """Base class for configs of forced measurement plugin algorithms."""
    pass

class ForcedAlgorithm(BaseAlgorithm):

    # All subclasses of ForcedAlgorithm should be registered here
    registry = AlgorithmRegistry(ForcedAlgorithmConfig)

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
        default=['forced.flux'
            ],
        doc="Plugin algorithms to be run and their configuration"
        )
    references = ConfigurableField(target=CoaddSrcReferencesTask, doc="Retrieve reference source catalog")

    copyColumns = DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent":"parentObjectId"}
        )


    """Config for ProcessCoadd"""
    coaddName = lsst.pex.config.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "goodSeeing",
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



class ForcedMeasurementTask(CmdLineTask):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclasses for running on CCDs and Coadds).
    """

    ConfigClass = ForcedMeasurementConfig
    RunnerClass = ButlerInitializedTaskRunner
    _DefaultName = "forcedMeasurementTask"
    dataPrefix = "deepCoadd_"
    #  The primary input to this init is the butler, which is a keyword argument, 
    #  
    def __init__(self, butler=None, refSchema=None, **kwds):
        CmdLineTask.__init__(self, **kwds)
        self.makeSubtask("references")
        self.algMetadata = lsst.daf.base.PropertyList()
        self.algorithms = AlgorithmMap()
        if not refSchema:
            refSchema = self.references.getSchema(butler)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        minimalSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.mapper.addMinimalSchema(minimalSchema)
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)

        flags = MeasurementDataFlags() 
        for executionOrder, name, config, AlgorithmClass in sorted(self.config.algorithms.apply()):
            self.algorithms[name] = AlgorithmClass(config, name, self.mapper, flags=flags,
                                                   others=self.algorithms, metadata=self.algMetadata)

    def run(self, dataRef):  #exposure, measCat, refCat, refWcs):
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)
        references = list(self.fetchReferences(dataRef, exposure))
        mapper = self.mapper
        retStruct = self.forcedMeasure(exposure, references, refWcs, dataRef)
        self.writeOutput(dataRef, retStruct.sources)
        return retStruct

    def forcedMeasure(self, exposure, references, refWcs, dataRef):

        targetWcs = exposure.getWcs()
        expregion = lsst.afw.geom.Box2I(exposure.getXY0(),
            lsst.afw.geom.Extent2I(exposure.getWidth(), exposure.getHeight()))
        footprints = {refRecord.getId(): (refRecord.getParent(), 
            refRecord.getFootprint().transform(refWcs, targetWcs, expregion, True)) 
            for refRecord in references} 
        for ref in references:
            if ref.getParent():
                if not ref.getParent() in footprints.keys():
                    footprints.pop(ref.getId())
        newref = list()
        
        for ref in references:
            if ref.getId() in footprints.keys():
                newref.append(ref)
        references = newref

        self.log.info("Performing forced measurement on %d sources" % len(references))

        # this is a catalog of just the references in out measurement area
        referenceCat = lsst.afw.table.SourceCatalog(self.mapper.getInputSchema())
        referenceCat.extend(references)

        sources = self.generateSources(dataRef, references)
        outSchema = self.mapper.getOutputSchema()
        # convert the footprints to the coordinate system of the exposure 


        noiseReplacer = NoiseReplacer(exposure, footprints, self.config.noiseSource,
           self.config.noiseOffset, self.config.noiseSeed, log=self.log)
        exposure.getMaskedImage().getImage().writeFits("replaced.fits")

        refParentCat, measParentCat = referenceCat.getChildren(0, sources)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat,measParentCat)):
            refChildCat, measChildCat = referenceCat.getChildren(refParentRecord.getId(), sources)
            # TODO: skip this loop if there are no algorithms configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                #noiseReplacer.insertSource(refChildRecord.getId())
                for algorithm in self.algorithms.iterSingle():
                    algorithm.measureSingle(exposure, measChildRecord, refChildRecord, refWcs)
                noiseReplacer.removeSource(refChildRecord.getId())
            noiseReplacer.insertSource(refParentRecord.getId())
            for algorithm in self.algorithms.iterSingle():
                algorithm.measureSingle(exposure, measParentRecord, refParentRecord, refWcs)
            for algorithm in self.algorithms.iterMulti():
                algorithm.measureMulti(exposure, measChildCat, refChildCat)
                algorithm.measureMulti(exposure, measParentCat[parentIdx:parentIdx+1],
                                       refParentCat[parentIdx:parentIdx+1])
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()

        #self.measurement.run(exposure, sources, references=references, refWcs=refWcs)
        return Struct(sources=sources)

    def makeIdFactory(self, dataRef):
        """Hook for derived classes to define how to make an IdFactory for forced sources.

        Note that this is for forced source IDs, not object IDs, which are usually handled by
        the copyColumns config option.
        """
        raise NotImplementedError()

    def fetchReferences(self, dataRef, exposure):
        """Hook for derived classes to define how to get references objects.

        Derived classes should call one of the fetch* methods on the references subtask,
        but which one they call depends on whether the region to get references for is a
        easy to describe in patches (as it would be when doing forced measurements on a
        coadd), or is just an arbitrary box (as it would be for CCD forced measurements).
        """
        raise NotImplementedError()


    def getExposure(self, dataRef):
        """Read input exposure on which to perform the measurements

        @param dataRef       Data reference from butler
        """
        return dataRef.get(self.dataPrefix + "calexp", immediate=True)

    def writeOutput(self, dataRef, sources):
        """Write forced source table

        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "forced_src")

    def generateSources(self, dataRef, references):
        """Generate sources to be measured, copying any fields in self.config.copyColumns

        @param dataRef     Data reference from butler
        @param references  Sequence (not necessarily a SourceCatalog) of reference sources
        @param idFactory   Factory to generate unique ids for forced sources
        @return Source catalog ready for measurement
        """
        idFactory = self.makeIdFactory(dataRef)
        table = lsst.afw.table.SourceTable.make(self.mapper.getOutputSchema(), idFactory)
        sources = lsst.afw.table.SourceCatalog(table)
        table = sources.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(references))
        for ref in references:
            sources.addNew().assign(ref, self.mapper)
        return sources

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def getPreviousTaskClass(self):
        return MeasureCoaddTask

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return None  #"%s_measureCoadd_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return None   #"%s_measureCoadd_metadata" % (self.config.coaddName,)
