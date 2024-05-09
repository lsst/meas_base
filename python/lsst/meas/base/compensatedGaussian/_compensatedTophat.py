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

from __future__ import annotations

__all__ = (
    "SingleFrameCompensatedTophatFluxConfig",
    "SingleFrameCompensatedTophatFluxPlugin",
)

import numpy as np
import math

from lsst.pex.config import RangeField, ListField
from lsst.geom import Point2I
import lsst.afw.geom

from ..sfm import SingleFramePlugin, SingleFramePluginConfig
from ..pluginRegistry import register

from .._measBaseLib import ApertureFluxAlgorithm, FlagHandler, FlagDefinitionList


class OutOfBoundsError(Exception):
    pass


class SingleFrameCompensatedTophatFluxConfig(SingleFramePluginConfig):
    apertures = ListField(
        doc="The aperture radii (in pixels) to measure the top-hats.",
        dtype=int,
        minLength=1,
        default=[12,],
    )
    inner_scale = RangeField(
        doc="Inner background annulus scale (relative to aperture).",
        dtype=float,
        default=1.5,
        min=1.0,
    )
    outer_scale = RangeField(
        doc="Outer background annulus scale (relative to aperture).",
        dtype=float,
        default=2.0,
        min=1.0,
    )

    def validate(self):
        super().validate()

        if self.outer_scale <= self.inner_scale:
            raise ValueError("The outer_scale must be greater than the inner_scale")


@register("base_CompensatedTophatFlux")
class SingleFrameCompensatedTophatFluxPlugin(SingleFramePlugin):
    ConfigClass = SingleFrameCompensatedTophatFluxConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def __init__(
        self,
        config: SingleFrameCompensatedTophatFluxConfig,
        name: str,
        schema,
        metadata,
        logName=None,
        **kwds,
    ):
        super().__init__(config, name, schema, metadata, logName, **kwds)

        flagDefs = FlagDefinitionList()

        self.aperture_keys = {}
        self._rads = {}
        self._inner_scale = config.inner_scale
        self._outer_scale = config.outer_scale
        for aperture in config.apertures:
            base_key = f"{name}_{aperture}"

            # flux
            flux_str = f"{base_key}_instFlux"
            flux_key = schema.addField(
                flux_str,
                type="D",
                doc="Compensated Tophat flux measurement.",
                units="count",
            )

            # flux error
            err_str = f"{base_key}_instFluxErr"
            err_key = schema.addField(
                err_str,
                type="D",
                doc="Compensated Tophat flux error.",
                units="count",
            )

            # mask bits
            mask_str = f"{base_key}_mask_bits"
            mask_key = schema.addField(mask_str, type=np.int32, doc="Mask bits set within aperture.")

            # failure flags
            failure_flag = flagDefs.add(f"{aperture}_flag", "Compensated Tophat measurement failed")
            oob_flag = flagDefs.add(f"{aperture}_flag_bounds", "Compensated Tophat out-of-bounds")

            self.aperture_keys[aperture] = (flux_key, err_key, mask_key, failure_flag, oob_flag)
            self._rads[aperture] = int(math.ceil(self._outer_scale*aperture))

        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)
        self._max_rad = max(self._rads)

    def fail(self, measRecord, error=None):
        """Record failure

        See also
        --------
        lsst.meas.base.SingleFramePlugin.fail
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()
        bbox = exposure.getBBox()

        y = center.getY() - bbox.beginY
        x = center.getX() - bbox.beginX

        y_floor = math.floor(y)
        x_floor = math.floor(x)

        ctrl = ApertureFluxAlgorithm.Control()

        for aperture, (flux_key, err_key, mask_key, failure_flag, oob_flag) in self.aperture_keys.items():
            rad = self._rads[aperture]

            if Point2I(center) not in exposure.getBBox().erodedBy(rad):
                self.flagHandler.setValue(measRecord, failure_flag.number, True)
                self.flagHandler.setValue(measRecord, oob_flag.number, True)
                continue

            y_slice = slice(y_floor - rad, y_floor + rad + 1, 1)
            x_slice = slice(x_floor - rad, x_floor + rad + 1, 1)

            # We will need the mask below, we can use this to test bounds as
            # well.
            sub_mask = exposure.mask.array[y_slice, x_slice]

            if sub_mask.size == 0 or sub_mask.shape[0] != sub_mask.shape[1] or (sub_mask.shape[0] % 2) == 0:
                self.flagHandler.setValue(measRecord, failure_flag.number, True)
                self.flagHandler.setValue(measRecord, oob_flag.number, True)
                continue

            # Compute three aperture fluxes.
            ellipse = lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(float(aperture),
                                                                        float(aperture), 0.0),
                                            center)
            tophat = ApertureFluxAlgorithm.computeFlux(exposure.maskedImage, ellipse, ctrl)
            ellipse.grow((self._inner_scale - 1.0)*aperture)
            inner = ApertureFluxAlgorithm.computeFlux(exposure.maskedImage, ellipse, ctrl)
            ellipse.grow((self._outer_scale - self._inner_scale)*aperture)
            outer = ApertureFluxAlgorithm.computeFlux(exposure.maskedImage, ellipse, ctrl)

            # We have flux in 3 circular apertures, a_0, a_1, a_2 with
            # associated variances \sigma_{a_0}^2, \sigma_{a_1}^2,
            # \sigma_{a_2)^2.
            # We transform these to annular fluxes:
            # b_0 = a_0
            # \sigma_{b_0}^2 = \sigma_{a_0}^2
            # b_1 = a_1 - a_0
            # \sigma_{b_1}^2 = \sigma_{a_1}^2 - \sigma_{a_0}^2
            # b_2 = a_2 - a_1
            # \sigma_{b_2}^2 = \sigma_{a_2}^2 - \sigma_{a_1}^2
            # Generally, the flux is then a weighted combination:
            # f = s_0*b_0 + s_1*b_1 + s_2*b_2
            # \sigma_f^2 = s_0^2*\sigma_{b_0}^2 + s_1^2*\sigma_{b_1}^2
            #              + s_2^2*\sigma_{b_2}^2
            # The inner aperture we use as-is, so s_0 = 1.0
            # We do not need the middle annulus, so s_1 = 0.0
            # The outer annulus is scaled by s_2 = -area_0 / (area_2 - area_1)

            a_0 = tophat.instFlux
            var_a_0 = tophat.instFluxErr*tophat.instFluxErr
            a_1 = inner.instFlux
            var_a_1 = inner.instFluxErr*inner.instFluxErr
            a_2 = outer.instFlux
            var_a_2 = outer.instFluxErr*outer.instFluxErr

            b_2 = a_2 - a_1
            var_b_2 = var_a_2 - var_a_1
            s_2 = 1.0/(self._outer_scale**2. - self._inner_scale**2.)

            flux = a_0 - s_2*b_2
            err = np.sqrt(var_a_0 + s_2*s_2*var_b_2)

            measRecord.set(flux_key, flux)
            measRecord.set(err_key, err)
            measRecord.set(mask_key, np.bitwise_or.reduce(sub_mask, axis=None))
