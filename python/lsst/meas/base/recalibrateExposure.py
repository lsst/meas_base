from lsst.pex.config import Config, Field, ChoiceField
from lsst.pipe.base import Task, Struct, PipelineTaskConnections
import lsst.pipe.base.connectionTypes as cT

from lsst.afw.image import PhotoCalib

__all__ = ("RecalibrateExposureConfig", "RecalibrateExposureTask", "RecalibrateExposureConnections")


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


class RecalibrateExposureConnections(PipelineTaskConnections,
                                     dimensions=("instrument", "visit", "detector", "skymap", "tract"),
                                     defaultTemplates={"calibSource": "jointcal"}):
    photoCalib = cT.PrerequisiteInput(
        doc="Photometric calibration",
        name="{calibSource}_photoCalib",
        storageClass="PhotoCalib",
        dimensions=["instrument", "visit", "detector", "tract"],
        multiple=True,
    )
    skyWcs = cT.PrerequisiteInput(
        doc="Astrometric calibration",
        name="{calibSource}_skyWcs",
        storageClass="SkyWcs",
        dimensions=["instrument", "visit", "detector", "tract"],
        multiple=True,
    )
    skyCorr = cT.PrerequisiteInput(
        doc="Sky correction",
        name="skyCorr",
        storageClass="Background",
        dimensions=["instrument", "visit", "detector"],
        multiple=True,
    )

    def __init__(self, *args, recalibrate="recalibrate", **kwargs):
        """Initialise

        Filters connections based on run-time task configuration.

        Parameters
        ----------
        recalibrate : `str`
            Name of the `RecalibrateExposureTask` in the ``config``.
        *args, **kwargs
            Usual construction arguments for `PipelineTaskConnections`.
        """
        super().__init__(*args, **kwargs)
        config = getattr(self.config, recalibrate)
        if not config.doApplyExternalPhotoCalib:
            self.prerequisiteInputs.discard("photoCalib")
        if not config.doApplyExternalSkyWcs:
            self.prerequisiteInputs.discard("skyWcs")
        if not config.doApplySkyCorr:
            self.prerequisiteInputs.discard("skyCorr")


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

        class FooConfig(lsst.pex.config.Config):
            recalibrate = lsst.pex.config.ConfigurableField(
                target=RecalibrateExposureTask,
                doc="Read calibrated exposure"
            )

        class FooTaskConnections(RecalibrateExposureConnections):
            ...

        class FooTask(lsst.pipe.base.PipelineTask):
            def runQuantum(self, butlerQC, inputRefs, outputRefs):
                inputs = butlerQC.get(inputRefs)
                calibs = self.recalibrate.extractCalibs(inputs, len(inputs["exposure"]))
                self.recalibrate.runMultiple(inputs["exposure"], **calibs)
                outputs = self.run(**inputs)
                ...
    """
    ConfigClass = RecalibrateExposureConfig

    def extractCalibs(self, inputs, number=None, single=False):
        """Extract only the datasets of interest

        This is a convenience function for the Gen3 middleware. The calibs are
        removed from the ``inputs``, so the ``inputs`` contains only those
        datasets for further operations.

        The result can be provided to the ``run`` method.

        Parameters
        ----------
        inputs : `dict` mapping `str` to some data
            Datasets that have been read from the butler. Will be modified,
            removing the calibs.
        number : `int`, optional
            Expected number of each type of calib.
        single : `bool`, optional
            Is there a single set only? If so, we'll return a single set of
            outputs, rather than a list. This allows passing the results to the
            ``run`` method rather than ``runMultiple``.

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

        Raises
        ------
        RuntimeError
            If the declared ``number`` of calibs doesn't match what's present.
        """
        if single:
            if number is not None and number != 1:
                raise RuntimeError(f"number ({number}) vs single ({single}) mismatch")
            number = 1
        calibs = {}
        for dataset in ("photoCalib", "skyWcs", "skyCorr"):
            data = inputs.pop(dataset, None)
            if data is not None and number is not None and len(data) != number:
                raise RuntimeError(f"Declared number of calibs ({number}) doesn't match actual number of "
                                   f"{dataset} ({len(data)})")
            calibs[dataset] = data
        if single:
            for name in calibs:
                if calibs[name] is None:
                    continue
                calibs[name] = calibs[name].pop()
        return calibs

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
            exposure.setPhotoCalib(PhotoCalib(scaling, photoCalib.getCalibrationErr()))
        else:
            exposure.setPhotoCalib(PhotoCalib(1.0))  # calibrateImage makes 1 ADU = 1 nJy everywhere
        if skyWcs is not None:
            exposure.setWcs(skyWcs)
        if skyCorr is not None:
            exposure.maskedImage -= skyCorr.getImage()
        return exposure

    def runMultiple(self, exposure, photoCalib=None, skyWcs=None, skyCorr=None, rescale=True):
        num = len(exposure)
        if photoCalib is None:
            photoCalib = [None]*num
        if skyWcs is None:
            skyWcs = [None]*num
        if skyCorr is None:
            skyCorr = [None]*num
        lengths = (len(photoCalib), len(skyWcs), len(skyCorr))
        if len(set(lengths)) != 1 or lengths[0] != num:
            raise RuntimeError(f"Length mismatch between exposure ({num}) and "
                               f"photoCalib, skyWcs, skyCorr {lengths}")
        for args in zip(exposure, photoCalib, skyWcs, skyCorr):
            self.run(*args, rescale=rescale)
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
        """
        calibs = self.readCalibs(dataRef)
        return self.run(exposure, rescale=rescale, **calibs.getDict())
