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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


__all__ = ("EvaluateLocalCalibrationConfig",
           "EvaluateLocalCalibrationTask")


class EvaluateLocalCalibrationConfig(pexConfig.Config):
    pass


class EvaluateLocalCalibrationTask(pipeBase.Task):
    """Compute local calibrations from an input calexp and store the values per
    row in a source catalog.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Input SourceTable schema.
    """
    ConfigClass = EvaluateLocalCalibrationConfig
    _DefaultName = "evaluateLocalCalibration"

    def __init__(self, schema, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.schema = schema

        self.photoKey = self.schema.addField(
            "base_localPhotoCalib",
            type="D",
            doc="Local approximation of the PhotoCalib calibration factor at "
                "the location of the src.")
        self.photoErrKey = self.schema.addField(
            "base_localPhotoCalibErr",
            type="D",
            doc="Error on the local approximation of the PhotoCalib "
                "calibration factor at the location of the src.")

        self.cdMatrix11Key = self.schema.addField(
            "base_CDMatrix_1_1",
            type="D",
            doc="(1, 1) element of the CDMatrix for the linear approximation "
                "of the WCS at the src location.")
        self.cdMatrix12Key = self.schema.addField(
            "base_CDMatrix_1_2",
            type="D",
            doc="(1, 2) element of the CDMatrix for the linear approximation "
                "of the WCS at the src location.")
        self.cdMatrix21Key = self.schema.addField(
            "base_CDMatrix_2_1",
            type="D",
            doc="(2, 1) element of the CDMatrix for the linear approximation "
                "of the WCS at the src location.")
        self.cdMatrix22Key = self.schema.addField(
            "base_CDMatrix_2_2",
            type="D",
            doc="(2, 2) element of the CDMatrix for the linear approximation "
                "of the WCS at the src location.")

    def run(self, sources, exposure):
        """Add calibration products to the source catalog given the calibrated
        exposure.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceTable`
            Catalog of sources to add local calibrations to.
        exposure : `lsst.afw.image.Exposure`
            Calibrated exposure to strip local calibrations from.
        """
        wcs = exposure.getWcs()
        photoCalib = exposure.getPhotoCalib()
        for srcRec in sources:
            self._storePhotoCalib(srcRec, photoCalib)
            self._storeWcs(srcRec, wcs)

    def _storePhotoCalib(self, srcRec, photoCalib):
        """Compute local calibration (counts -> Nanojansky) for the input
        source record.

        Parameters
        ----------
        srcRec : `lsst.afw.table.SourceRecord`
            Record with centroid to store local calibrations in.
        photoCalib : `lsst.afw.image.PhotoCalib`
            Photometric calibration object to compute local calibration from.
        """
        calib = photoCalib.getLocalCalibration(srcRec.getCentroid())
        srcRec.set(self.photoKey, calib)

        calibErr = photoCalib.getCalibrationErr()
        srcRec.set(self.photoErrKey, calibErr)

    def _storeWcs(self, srcRec, wcs):
        """Compute and store the CDMatrix at the location of the source record.

        Parameters
        ----------
        srcRec : `lsst.afw.table.SourceRecord`
            Source record with centroid slot set.
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs to compute a linear approximation of.
        """
        localCDMatrix = wcs.getCdMatrix(srcRec.getCentroid())
        srcRec.set(self.cdMatrix11Key, localCDMatrix[0, 0])
        srcRec.set(self.cdMatrix12Key, localCDMatrix[0, 1])
        srcRec.set(self.cdMatrix21Key, localCDMatrix[1, 0])
        srcRec.set(self.cdMatrix22Key, localCDMatrix[1, 1])
