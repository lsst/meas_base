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

r"""Base classes for forced measurement plugins and the driver task for these.

In forced measurement, a reference catalog is used to define restricted
measurements (usually just fluxes) on an image.  As the reference catalog may
be deeper than the detection limit of the measurement image, we do not assume
that we can use detection and deblend information from the measurement image.
Instead, we assume this information is present in the reference catalog and
can be "transformed" in some sense to the measurement frame.  At the very
least, this means that `~lsst.afw.detection.Footprint`\ s from the reference
catalog should be transformed and installed as Footprints in the output
measurement catalog. If we have a procedure that can transform "heavy"
Footprints (ie, including pixel data), we can then proceed with measurement as
usual, but using the reference catalog's ``id`` and ``parent`` fields to
define deblend families.  If this transformation does not preserve
heavy Footprints (this is currently the case, at least for CCD forced
photometry), then we will only be able to replace objects with noise one
deblend family at a time, and hence measurements run in single-object mode may
be contaminated by neighbors when run on objects with ``parent != 0``.

Measurements are generally recorded in the coordinate system of the image
being measured (and all slot-eligible fields must be), but non-slot fields may
be recorded in other coordinate systems if necessary to avoid information loss
(this should, of course, be indicated in the field documentation).  Note that
the reference catalog may be in a different coordinate system; it is the
responsibility of plugins to transform the data they need themselves, using
the reference WCS provided.  However, for plugins that only require a position
or shape, they may simply use output `~lsst.afw.table.SourceCatalog`\'s
centroid or shape slots, which will generally be set to the transformed
position of the reference object before any other plugins are run, and hence
avoid using the reference catalog at all.

Command-line driver tasks for forced measurement include
`ForcedPhotImageTask`, `ForcedPhotCcdTask`, and `ForcedPhotCoaddTask`.
"""

import lsst.pex.config
import lsst.pipe.base

from .pluginRegistry import PluginRegistry
from .baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
                              BaseMeasurementConfig, BaseMeasurementTask)
from .noiseReplacer import NoiseReplacer, DummyNoiseReplacer

__all__ = ("ForcedPluginConfig", "ForcedPlugin",
           "ForcedMeasurementConfig", "ForcedMeasurementTask")


class ForcedPluginConfig(BaseMeasurementPluginConfig):
    """Base class for configs of forced measurement plugins."""

    pass


class ForcedPlugin(BaseMeasurementPlugin):
    """Base class for forced measurement plugins.

    Parameters
    ----------
    config : `ForcedPlugin.ConfigClass`
        Configuration for this plugin.
    name : `str`
        The string with which the plugin was registered.
    schemaMapper : `lsst.afw.table.SchemaMapper`
        A mapping from reference catalog fields to output catalog fields.
        Output fields should be added to the output schema. While most plugins
        will not need to map fields from the reference schema, if they do so,
        those fields will be transferred before any plugins are run.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog.
    logName : `str`, optional
        Name to use when logging errors.
    """

    registry = PluginRegistry(ForcedPluginConfig)
    """Subclasses of `ForcedPlugin` must be registered here (`PluginRegistry`).
    """

    ConfigClass = ForcedPluginConfig

    def __init__(self, config, name, schemaMapper, metadata, logName=None):
        BaseMeasurementPlugin.__init__(self, config, name, logName=logName)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        """Measure the properties of a source given an image and a reference.

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            The pixel data to be measured, together with the associated PSF,
            WCS, etc. All other sources in the image should have been replaced
            by noise according to deblender outputs.
        measRecord : `lsst.afw.table.SourceRecord`
            Record describing the object being measured. Previously-measured
            quantities will be retrieved from here, and it will be updated
            in-place with the outputs of this plugin.
        refRecord : `lsst.afw.table.SimpleRecord`
            Additional parameters to define the fit, as measured elsewhere.
        refWcs : `lsst.afw.geom.SkyWcs` or `lsst.afw.geom.Angle`
            The coordinate system for the reference catalog values. An
            `~lsst.geom.Angle` may be passed, indicating that a local tangent
            WCS should be created for each object using the given angle as a
            pixel scale.

        Notes
        -----
        In the normal mode of operation, the source centroid will be set to
        the WCS-transformed position of the reference object, so plugins that
        only require a reference position should not have to access the
        reference object at all.
        """
        raise NotImplementedError()

    def measureN(self, measCat, exposure, refCat, refWcs):
        """Measure the properties of blended sources from image & reference.

        This operates on all members of a blend family at once.

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            The pixel data to be measured, together with the associated PSF,
            WCS, etc. Sources not in the blended hierarchy to be measured
            should have been replaced with noise using deblender outputs.
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog describing the objects (and only those objects) being
            measured. Previously-measured quantities will be retrieved from
            here, and it will be updated in-place with the outputs of this
            plugin.
        refCat : `lsst.afw.table.SimpleCatalog`
            Additional parameters to define the fit, as measured elsewhere.
            Ordered such that ``zip(measCat, refcat)`` may be used.
        refWcs : `lsst.afw.geom.SkyWcs` or `lsst.afw.geom.Angle`
            The coordinate system for the reference catalog values. An
            `~lsst.geom.Angle` may be passed, indicating that a local tangent
            WCS should be created for each object using the given angle as a
            pixel scale.

        Notes
        -----
        In the normal mode of operation, the source centroids will be set to
        the WCS-transformed position of the reference object, so plugins that
        only require a reference position should not have to access the
        reference object at all.
        """
        raise NotImplementedError()


class ForcedMeasurementConfig(BaseMeasurementConfig):
    """Config class for forced measurement driver task.
    """

    plugins = ForcedPlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_TransformedCentroid",
                 "base_SdssCentroid",
                 "base_TransformedShape",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_CircularApertureFlux",
                 "base_PsfFlux",
                 "base_LocalBackground",
                 ],
        doc="Plugins to be run and their configuration"
    )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")
    undeblended = ForcedPlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )

    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent": "parentObjectId", "deblend_nChild": "deblend_nChild",
                 "coord_ra": "coord_ra", "coord_dec": "coord_dec",
                 "base_ClassificationExtendedness_value": "ref_base_ClassificationExtendedness_value",
                 "base_ClassificationExtendedness_flag": "ref_base_ClassificationExtendedness_flag",
                 "calib_psf_used": "ref_calib_psf_used", "calib_photometry_used": "ref_calib_photometry_used",
                 "calib_astrometry_used": "ref_calib_astrometry_used",
                 "base_PsfFlux_instFlux": "ref_base_PsfFlux_instFlux",
                 "base_PsfFlux_instFluxErr": "ref_base_PsfFlux_instFluxErr",
                 "base_PsfFlux_flag": "ref_base_PsfFlux_flag",
                 "modelfit_CModel_instFlux": "ref_modelfit_CModel_instFlux",
                 "modelfit_CModel_instFluxErr": "ref_modelfit_CModel_instFluxErr",
                 "modelfit_CModel_flag": "ref_modelfit_CModel_flag",
                 "base_CircularApertureFlux_12_0_instFlux": "ref_base_CircularApertureFlux_12_0_instFlux",
                 "base_CircularApertureFlux_12_0_instFluxErr":
                     "ref_base_CircularApertureFlux_12_0_instFluxErr",
                 "base_CircularApertureFlux_12_0_flag":
                     "ref_base_CircularApertureFlux_12_0_flag",
                 "base_Blendedness_raw": "ref_base_Blendedness_raw",
                 "detect_isPrimary": "ref_detect_isPrimary"}
    )

    checkUnitsParseStrict = lsst.pex.config.Field(
        doc="Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'",
        dtype=str,
        default="raise",
    )

    def setDefaults(self):
        self.slots.centroid = "base_TransformedCentroid"
        self.slots.shape = "base_TransformedShape"
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = None
        self.slots.gaussianFlux = None
        self.slots.calibFlux = None


class ForcedMeasurementTask(BaseMeasurementTask):
    """Measure sources on an image, constrained by a reference catalog.

    A subtask for measuring the properties of sources on a single image,
    using an existing "reference" catalog to constrain some aspects of the
    measurement.

    Parameters
    ----------
    refSchema : `lsst.afw.table.Schema`
        Schema of the reference catalog.  Must match the catalog later passed
        to 'ForcedMeasurementTask.generateMeasCat` and/or
        `ForcedMeasurementTask.run`.
    algMetadata : `lsst.daf.base.PropertyList` or `None`
        Will be updated in place to to record information about each
        algorithm. An empty `~lsst.daf.base.PropertyList` will be created if
        `None`.
    **kwds
        Keyword arguments are passed to the supertask constructor.

    Notes
    -----
    Note that while `SingleFrameMeasurementTask` is passed an initial
    `~lsst.afw.table.Schema` that is appended to in order to create the output
    `~lsst.afw.table.Schema`, `ForcedMeasurementTask` is initialized with the
    `~lsst.afw.table.Schema` of the reference catalog, from which a new
    `~lsst.afw.table.Schema` for the output catalog is created.  Fields to be
    copied directly from the reference `~lsst.afw.table.Schema` are added
    before ``Plugin`` fields are added.
    """

    ConfigClass = ForcedMeasurementConfig

    def __init__(self, refSchema, algMetadata=None, **kwds):
        super(ForcedMeasurementTask, self).__init__(algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(), False)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)
        # Copy over all merge_measurement_FILTER and merge_peak_sky columns
        mergeKeys = refSchema.extract("merge_measurement*")
        mergeKeys.update(refSchema.extract("merge_peak_sky"))
        for mergeName in mergeKeys.keys():
            refItem = refSchema.find(mergeName)
            self.mapper.addMapping(refItem.key, mergeName)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        self.initializePlugins(schemaMapper=self.mapper)
        self.schema = self.mapper.getOutputSchema()
        self.schema.checkUnits(parse_strict=self.config.checkUnitsParseStrict)

    def run(self, measCat, exposure, refCat, refWcs, exposureId=None, beginOrder=None, endOrder=None):
        r"""Perform forced measurement.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposureF`
            Image to be measured. Must have at least a `lsst.afw.geom.SkyWcs`
            attached.
        measCat : `lsst.afw.table.SourceCatalog`
            Source catalog for measurement results; must be initialized with
            empty records already corresponding to those in ``refCat`` (via
            e.g. `generateMeasCat`).
        refCat : `lsst.afw.table.SourceCatalog`
            A sequence of `lsst.afw.table.SourceRecord` objects that provide
            reference information for the measurement.  These will be passed
            to each plugin in addition to the output
            `~lsst.afw.table.SourceRecord`.
        refWcs : `lsst.afw.geom.SkyWcs`
            Defines the X,Y coordinate system of ``refCat``.
        exposureId : `int`, optional
            Optional unique exposureId used to calculate random number
            generator seed in the NoiseReplacer.
        beginOrder : `int`, optional
            Beginning execution order (inclusive). Algorithms with
            ``executionOrder`` < ``beginOrder`` are not executed. `None` for no limit.
        endOrder : `int`, optional
            Ending execution order (exclusive). Algorithms with
            ``executionOrder`` >= ``endOrder`` are not executed. `None` for no limit.

        Notes
        -----
        Fills the initial empty `~lsst.afw.table.SourceCatalog` with forced
        measurement results.  Two steps must occur before `run` can be called:

        - `generateMeasCat` must be called to create the output ``measCat``
          argument.
        - `~lsst.afw.detection.Footprint`\ s appropriate for the forced sources
          must be attached to the ``measCat`` records. The
          `attachTransformedFootprints` method can be used to do this, but
          this degrades "heavy" (i.e., including pixel values)
          `~lsst.afw.detection.Footprint`\s to regular
          `~lsst.afw.detection.Footprint`\s, leading to non-deblended
          measurement, so most callers should provide
          `~lsst.afw.detection.Footprint`\s some other way. Typically, calling
          code will have access to information that will allow them to provide
          heavy footprints - for instance, `ForcedPhotCoaddTask` uses the
          heavy footprints from deblending run in the same band just before
          non-forced is run measurement in that band.
        """
        # First check that the reference catalog does not contain any children
        # for which any member of their parent chain is not within the list.
        # This can occur at boundaries when the parent is outside and one of
        # the children is within.  Currently, the parent chain is always only
        # one deep, but just in case, this code checks for any case where when
        # the parent chain to a child's topmost parent is broken and raises an
        # exception if it occurs.
        #
        # I.e. this code checks that this precondition is satisfied by
        # whatever reference catalog provider is being paired with it.
        refCatIdDict = {ref.getId(): ref.getParent() for ref in refCat}
        for ref in refCat:
            refId = ref.getId()
            topId = refId
            while(topId > 0):
                if topId not in refCatIdDict:
                    raise RuntimeError("Reference catalog contains a child for which at least "
                                       "one parent in its parent chain is not in the catalog.")
                topId = refCatIdDict[topId]

        # Construct a footprints dict which looks like
        # {ref.getId(): (ref.getParent(), source.getFootprint())}
        # (i.e. getting the footprint from the transformed source footprint)
        footprints = {ref.getId(): (ref.getParent(), measRecord.getFootprint())
                      for (ref, measRecord) in zip(refCat, measCat)}

        self.log.info("Performing forced measurement on %d source%s", len(refCat),
                      "" if len(refCat) == 1 else "s")

        if self.config.doReplaceWithNoise:
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, exposure,
                                          footprints, log=self.log, exposureId=exposureId)
            algMetadata = measCat.getTable().getMetadata()
            if algMetadata is not None:
                algMetadata.addInt("NOISE_SEED_MULTIPLIER", self.config.noiseReplacer.noiseSeedMultiplier)
                algMetadata.addString("NOISE_SOURCE", self.config.noiseReplacer.noiseSource)
                algMetadata.addDouble("NOISE_OFFSET", self.config.noiseReplacer.noiseOffset)
                if exposureId is not None:
                    algMetadata.addLong("NOISE_EXPOSURE_ID", exposureId)
        else:
            noiseReplacer = DummyNoiseReplacer()

        # Create parent cat which slices both the refCat and measCat (sources)
        # first, get the reference and source records which have no parent
        refParentCat, measParentCat = refCat.getChildren(0, measCat)
        for parentIdx, (refParentRecord, measParentRecord) in enumerate(zip(refParentCat, measParentCat)):

            # first process the records which have the current parent as children
            refChildCat, measChildCat = refCat.getChildren(refParentRecord.getId(), measCat)
            # TODO: skip this loop if there are no plugins configured for single-object mode
            for refChildRecord, measChildRecord in zip(refChildCat, measChildCat):
                noiseReplacer.insertSource(refChildRecord.getId())
                self.callMeasure(measChildRecord, exposure, refChildRecord, refWcs,
                                 beginOrder=beginOrder, endOrder=endOrder)
                noiseReplacer.removeSource(refChildRecord.getId())

            # then process the parent record
            noiseReplacer.insertSource(refParentRecord.getId())
            self.callMeasure(measParentRecord, exposure, refParentRecord, refWcs,
                             beginOrder=beginOrder, endOrder=endOrder)
            self.callMeasureN(measParentCat[parentIdx:parentIdx+1], exposure,
                              refParentCat[parentIdx:parentIdx+1],
                              beginOrder=beginOrder, endOrder=endOrder)
            # measure all the children simultaneously
            self.callMeasureN(measChildCat, exposure, refChildCat,
                              beginOrder=beginOrder, endOrder=endOrder)
            noiseReplacer.removeSource(refParentRecord.getId())
        noiseReplacer.end()

        # Undeblended plugins only fire if we're running everything
        if endOrder is None:
            for measRecord, refRecord in zip(measCat, refCat):
                for plugin in self.undeblendedPlugins.iter():
                    self.doMeasurement(plugin, measRecord, exposure, refRecord, refWcs)

    def generateMeasCat(self, exposure, refCat, refWcs, idFactory=None):
        r"""Initialize an output catalog from the reference catalog.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposureF`
            Image to be measured.
        refCat : iterable of `lsst.afw.table.SourceRecord`
            Catalog of reference sources.
        refWcs : `lsst.afw.geom.SkyWcs`
            Defines the X,Y coordinate system of ``refCat``.
        idFactory : `lsst.afw.table.IdFactory`, optional
            Factory for creating IDs for sources.

        Returns
        -------
        meascat : `lsst.afw.table.SourceCatalog`
            Source catalog ready for measurement.

        Notes
        -----
        This generates a new blank `~lsst.afw.table.SourceRecord` for each
        record in ``refCat``. Note that this method does not attach any
        `~lsst.afw.detection.Footprint`\ s.  Doing so is up to the caller (who
        may call `attachedTransformedFootprints` or define their own method -
        see `run` for more information).
        """
        if idFactory is None:
            idFactory = lsst.afw.table.IdFactory.makeSimple()
        table = lsst.afw.table.SourceTable.make(self.schema, idFactory)
        measCat = lsst.afw.table.SourceCatalog(table)
        table = measCat.table
        table.setMetadata(self.algMetadata)
        table.preallocate(len(refCat))
        for ref in refCat:
            newSource = measCat.addNew()
            newSource.assign(ref, self.mapper)
        return measCat

    def attachTransformedFootprints(self, sources, refCat, exposure, refWcs):
        r"""Attach Footprints to blank sources prior to measurement.

        Notes
        -----
        `~lsst.afw.detection.Footprint`\s for forced photometry must be in the
        pixel coordinate system of the image being measured, while the actual
        detections may start out in a different coordinate system. This
        default implementation transforms the Footprints from the reference
        catalog from the WCS to the exposure's WCS, which downgrades
        ``HeavyFootprint``\s into regular `~lsst.afw.detection.Footprint`\s,
        destroying deblend information.

        Note that `ForcedPhotImageTask` delegates to this method in its own
        `~ForcedPhotImageTask.attachFootprints` method.  This method can be
        overridden by its subclasses to define how their
        `~lsst.afw.detection.Footprint`\s should be generated.

        See the documentation for `run` for information about the
        relationships between `run`, `generateMeasCat`, and
        `attachTransformedFootprints`.
        """
        exposureWcs = exposure.getWcs()
        region = exposure.getBBox(lsst.afw.image.PARENT)
        for srcRecord, refRecord in zip(sources, refCat):
            srcRecord.setFootprint(refRecord.getFootprint().transform(refWcs, exposureWcs, region))
