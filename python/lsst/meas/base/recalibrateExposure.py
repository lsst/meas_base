from lsst.pex.config import Config, Field, ChoiceField
from lsst.pipe.base import Task, Struct
import lsst.pipe.base.connectionTypes as cT

from lsst.afw.image import PhotoCalib

__all__ = ("RecalibrateExposureConfig", "RecalibrateExposureTask")


class MissingExposureError(Exception):
    """Raised when data cannot be retrieved for an exposure.

    When processing patches, sometimes one exposure is missing; this lets us
    distinguish bewteen that case, and other errors.
    """
    pass


class RecalibrateExposureConfig(Config):
    doApplyExternalPhotoCalib = Field(
        dtype=bool,
        default=False,
        doc=("Whether to apply external photometric calibration via an "
             "`lsst.afw.image.PhotoCalib` object. Uses the "
             "``externalPhotoCalibName`` field to determine which calibration "
             "to load."),
    )
    doApplyExternalSkyWcs = Field(
        dtype=bool,
        default=False,
        doc=("Whether to apply external astrometric calibration via an "
             "`lsst.afw.geom.SkyWcs` object. Uses ``externalSkyWcsName`` "
             "field to determine which calibration to load."),
    )
    doApplySkyCorr = Field(
        dtype=bool,
        default=False,
        doc="Apply sky correction?",
    )
    includePhotoCalibVar = Field(
        dtype=bool,
        default=False,
        doc="Add photometric calibration variance to warp variance plane?",
    )
    externalPhotoCalibName = ChoiceField(
        dtype=str,
        doc=("Type of external PhotoCalib if ``doApplyExternalPhotoCalib`` is True. "
             "Unused for Gen3 middleware."),
        default="jointcal",
        allowed={
            "jointcal": "Use jointcal_photoCalib",
            "fgcm": "Use fgcm_photoCalib",
            "fgcm_tract": "Use fgcm_tract_photoCalib"
        },
    )
    externalSkyWcsName = ChoiceField(
        dtype=str,
        doc="Type of external SkyWcs if ``doApplyExternalSkyWcs`` is True. Unused for Gen3 middleware.",
        default="jointcal",
        allowed={
            "jointcal": "Use jointcal_wcs"
        },
    )


class RecalibrateExposureTask(Task):
    """Read and apply calibrations to a exposure

    This is used to apply calibrations from jointcal or fgcm to a exposure.

    If ``config.doApplyExternalPhotoCalib`` is ``True``, the photometric
    calibration (``photoCalib``) is taken from ``config.externalPhotoCalibName``
    via the ``<name>_photoCalib`` dataset.  Otherwise, the photometric
    calibration is retrieved from the processed exposure.

    When ``config.doApplyExternalSkyWcs` is ``True``, the astrometric
    calibration is taken from ``config.externalSkyWcsName`` with the
    ``<name>_wcs`` dataset. Otherwise, the astrometric calibration is taken from
    the processed exposure.

    If ``config.doApplySkyCorr`` is ``True``, the sky background is read from
    dataset ``skyCorr``.

    Use in Gen3
    -----------

    Wrap your connections class in the decorator provided by
    ``RecalibrateExposureTask.addConnections``, and add a template for
    ``calibSource``. Then when in ``runQuantum`` you call
    ``butlerQC.get(inputRefs)``, you will get the calibrations, and you can
    use the ``extractCalibs`` and ``run`` methods to apply the calibrations.

        class FooConfig(lsst.pex.config.Config):
            recalibrate = lsst.pex.config.ConfigurableField(
                target=RecalibrateExposureTask,
                doc="Read calibrated exposure"
            )

        @RecalibrateExposureTask.addConnections()
        class FooTaskConnections(lsst.pipe.base.PipelineTaskConnections,
                                dimensions=("instrument", "visit", "detector",
                                            "skymap", "tract"),
                                defaultTemplates={"calibSource": "jointcal"}):
            ...

        class FooTask(lsst.pipe.base.PipelineTask):
            def runQuantum(self, butlerQC, inputRefs, outputRefs):
                inputs = butlerQC.get(inputRefs)
                calibs = self.recalibrate.extractCalibs(**inputs)
                inputs["exposure"] = self.recalibrate.run(inputs["exposure"],
                                                          **calibs)
                outputs = self.run(**inputs)
                butlerQC.put(outputs, outputRefs)
    """
    ConfigClass = RecalibrateExposureConfig

    @classmethod
    def addConnections(cls, recalibrate="recalibrate", **kwargs):
        """Return a decorator for subclasses of `PipelineTaskConnections` to
        add connections for Gen3 middleware

        This method is intended for use within the Gen3 middleware only.

        Parameters
        ----------
        recalibrate : `str`
            Name of the ``RecalibrateExposureTask`` in the task
            configuration.
        **kwargs : `dict`
            Additional arguments to use for all the connections, e.g.,
            ``multiple=True``.

        Returns
        -------
        decorator : callable
            Class decorator for `lsst.pipe.base.PipelineTaskConnections`
            subclasses.
        """
        def decorator(ConnectionsClass):
            """Class decorator to add connections required for
            `RecalibrateExposureTask`

            The ``__init__`` method is also wrapped to filter the connections
            using the run-time task configuration.

            Parameters
            ----------
            ConnectionsClass : subclass of `lsst.pipe.base.PipelineTaskConnections`
                Class to which to add connections.

            Returns
            -------
            ConnectionsClass : subclass of `lsst.pipe.base.PipelineTaskConnections`
                Class with added connections.
            """
            connections = dict(
                photoCalib=cT.Input(
                    doc="Photometric calibration",
                    name="{calibSource}_photoCalib",
                    storageClass="PhotoCalib",
                    dimensions=["instrument", "visit", "detector", "tract"],
                    **kwargs
                ),
                skyWcs=cT.Input(
                    doc="Astrometric calibration",
                    name="{calibSource}_skyWcs",
                    storageClass="SkyWcs",
                    dimensions=["instrument", "visit", "detector", "tract"],
                    **kwargs
                ),
                skyCorr=cT.Input(
                    doc="Sky correction",
                    name="skyCorr",
                    storageClass="Background",
                    dimensions=["instrument", "visit", "detector"],
                    **kwargs
                ),
            )
            for conn in connections:
                setattr(ConnectionsClass, conn, connections[conn])

            # Wrap the __init__ so it will run our 'filterConnnections' method
            originalInit = ConnectionsClass.__init__

            def init(self, *args, **kwargs):
                originalInit(self, *args, **kwargs)
                cls.filterConnections(self, getattr(self.config, recalibrate))

            ConnectionsClass.__init__ = init

            return ConnectionsClass
        return decorator

    @classmethod
    def filterConnections(cls, connections, config):
        """Filter connections based on run-time task configuration

        This method is intended for use within the Gen3 middleware only.

        Parameters
        ----------
        connections : subclass of `lsst.pipe.base.PipelineTaskConnections`
            Declared datasets of interest.
        config : subclass of `lsst.pex.config.Config`
            Run-time task configuration.
        """
        if not config.doApplyExternalPhotoCalib:
            connections.inputs.discard("photoCalib")
        if not config.doApplyExternalSkyWcs:
            connections.inputs.discard("skyWcs")
        if not config.doApplySkyCorr:
            connections.inputs.discard("skyCorr")

    def extractCalibs(self, inputs):
        """Extract only the datasets of interest

        The result can be provided to the ``run`` method.

        Parameters
        ----------
        inputs : `dict` mapping `str` to some data
            Datasets that have been read from the butler. Will be modified,
            removing the calibs.

        Returns
        -------
        calibs : `dict`
            Calibration inputs for the ``run`` method:

            ``photoCalib``
                Photometric calibration (`lsst.afw.image.PhotoCalib`), or
                ``None`` if not ``doApplyExternalPhotoCalib``.
            ``skyWcs``
                Astrometric calibration (`lsst.afw.geom.SkyWcs`), or ``None``
                if not ``doApplyExternalSkyWcs``.
            ``skyCorr``
                Sky correction to apply, or ``None`` if not ``doApplySkyCorr``.
        """
        return {name: inputs.pop(name, None) for name in ("photoCalib", "skyWcs", "skyCorr")}

    def readCalibs(self, dataRef):
        """Read all inputs

        This includes the exposure, photometric and astrometric calibrations,
        and sky correction.

        This method is intended for use with the Gen2 middleware only.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Inputs that were read, including the following elements:

            ``photoCalib``
                Photometric calibration (`lsst.afw.image.PhotoCalib`), or
                ``None`` if not ``doApplyExternalPhotoCalib``.
            ``skyWcs``
                Astrometric calibration (`lsst.afw.geom.SkyWcs`), or ``None``
                if not ``doApplyExternalSkyWcs``.
            ``skyCorr``
                Sky correction to apply, or ``None`` if not ``doApplySkyCorr``.
        """
        photoCalib = None
        if self.config.doApplyExternalPhotoCalib:
            photoCalib = dataRef.get(f"{self.config.externalPhotoCalibName}_photoCalib")

        skyWcs = None
        if self.config.doApplyExternalSkyWcs:
            skyWcs = dataRef.get(f"{self.config.externalSkyWcsName}_wcs")

        skyCorr = None
        if self.config.doApplySkyCorr:
            skyCorr = dataRef.get("skyCorr")

        return Struct(photoCalib=photoCalib, skyWcs=skyWcs, skyCorr=skyCorr)

    def run(self, exposure, photoCalib=None, skyWcs=None, skyCorr=None, rescale=True):
        """Apply calibrations to an exposure

        The result is an exposure with units of nJy, unless ``rescale``, in
        which case the mean calibration is used.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to recalibrate. Will be modified.
        photoCalib : `lsst.afw.image.PhotoCalib`, optional
            Photometric calibration, or ``None`` if no photometric calibration
            is to be applied.
        skyWcs : `lsst.afw.geom.SkyWcs`, optional
            Astrometric calibration, or ``None`` if no astrometric calibration
            is to be applied.
        skyCorr : `lsst.afw.math.BackgroundList`, optional
            Sky correction, or ``None`` if no sky correction is to be applied.
        rescale : `bool`
            Rescale to mean zero-point?

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            Recalibrated exposure.
        """
        if photoCalib is None:
            photoCalib = exposure.getPhotoCalib()
        exposure.maskedImage = photoCalib.calibrateImage(exposure.maskedImage,
                                                         self.config.includePhotoCalibVar)
        if rescale:
            scaling = photoCalib.getCalibrationMean()
            exposure.maskedImage /= scaling
            exposure.setPhotoCalib(PhotoCalib(1.0/scaling))
        else:
            exposure.setPhotoCalib(PhotoCalib(1.0))  # calibrateImage makes 1 ADU = 1 nJy everywhere
        if skyWcs is not None:
            exposure.setWcs(skyWcs)
        if skyCorr is not None:
            exposure.maskedImage -= skyCorr.getImage()
        return exposure

    def runDataRef(self, dataRef, exposure, rescale=True):
        """Return a re-calibrated exposure

        This method is intended for use with the Gen2 middleware only.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        exposure : `lsst.afw.image.Exposure`
            Exposure to recalibrate.
        rescale : `bool`
            Rescale to an exposure with units of nJy?

        Returns
        -------
        exposure : `lsst.afw.image.Exposure`
            Re-calibrated exposure.

        Raises
        ------
        MissingExposureError
            If data for the exposure is not available.
        """
        calibs = self.readCalibs(dataRef)
        return self.run(exposure, rescale=rescale, **calibs.getDict())
