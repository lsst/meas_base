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
    "SingleFrameCompensatedGaussianFluxConfig",
    "SingleFrameCompensatedGaussianFluxPlugin",
)

import math
import numpy as np
import scipy.stats as sps

from lsst.pex.config import Field, ListField
from lsst.geom import Point2I

from ..sfm import SingleFramePlugin, SingleFramePluginConfig
from ..pluginRegistry import register

from .._measBaseLib import _compensatedGaussianFiltInnerProduct


class OutOfBoundsError(Exception):
    pass


class SingleFrameCompensatedGaussianFluxConfig(SingleFramePluginConfig):
    kernel_widths = ListField(
        doc="The widths (in pixels) of the kernels for which to measure compensated apertures.",
        dtype=int,
        minLength=1,
        default=[3, 5]
    )

    t = Field(
        doc="Scale parameter of outer Gaussian compared to inner Gaussian.",
        dtype=float,
        default=2.0,
    )


@register("base_CompensatedGaussianFlux")
class SingleFrameCompensatedGaussianFluxPlugin(SingleFramePlugin):
    ConfigClass = SingleFrameCompensatedGaussianFluxConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def __init__(
        self,
        config: SingleFrameCompensatedGaussianFluxConfig,
        name: str,
        schema,
        metadata,
        logName=None,
        **kwds,
    ):
        super().__init__(config, name, schema, metadata, logName, **kwds)

        # create generic failure key
        self.fatalFailKey = schema.addField(
            f"{name}_flag", type="Flag", doc="Set to 1 for any fatal failure."
        )

        # Out of bounds failure key
        self.ooBoundsKey = schema.addField(
            f"{name}_bounds_flag",
            type="Flag",
            doc="Flag set to 1 if not all filters fit within exposure.",
        )

        self.width_keys = {}
        self._rads = {}
        self._flux_corrections = {}
        self._variance_corrections = {}
        self._t = config.t
        for width in config.kernel_widths:
            base_key = f"{name}_{width}"

            # flux
            flux_str = f"{base_key}_instFlux"
            flux_key = schema.addField(
                flux_str,
                type="D",
                doc="Compensated Gaussian flux measurement.",
                units="count",
            )

            # flux error
            err_str = f"{base_key}_instFluxErr"
            err_key = schema.addField(
                err_str,
                type="D",
                doc="Compensated Gaussian flux error.",
                units="count",
            )

            # mask bits
            mask_str = f"{base_key}_mask_bits"
            mask_key = schema.addField(mask_str, type=np.int32, doc="Mask bits set within aperture.")

            self.width_keys[width] = (flux_key, err_key, mask_key)
            self._rads[width] = math.ceil(sps.norm.ppf((0.995,), scale=width * config.t)[0])

        self._max_rad = max(self._rads)

    def fail(self, measRecord, error=None):
        if isinstance(error, OutOfBoundsError):
            measRecord.set(self.ooBoundsKey, True)
        measRecord.set(self.fatalFailKey, True)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()
        bbox = exposure.getBBox()

        if Point2I(center) not in exposure.getBBox().erodedBy(self._max_rad):
            raise OutOfBoundsError("Not all the kernels for this source fit inside the exposure.")

        y = center.getY() - bbox.beginY
        x = center.getX() - bbox.beginX

        y_floor = math.floor(y)
        x_floor = math.floor(x)

        for width, (flux_key, err_key, mask_key) in self.width_keys.items():
            rad = self._rads[width]
            y_slice = slice(y_floor - rad, y_floor + rad + 1, 1)
            x_slice = slice(x_floor - rad, x_floor + rad + 1, 1)
            y_mean = y - y_floor + rad
            x_mean = x - x_floor + rad

            flux, var = _compensatedGaussianFiltInnerProduct(
                exposure.image.array[y_slice, x_slice],
                exposure.variance.array[y_slice, x_slice],
                x_mean,
                y_mean,
                width,
                self._t,
            )
            measRecord.set(flux_key, flux)
            measRecord.set(err_key, np.sqrt(var))
            measRecord.set(mask_key, np.bitwise_or.reduce(exposure.mask.array[y_slice, x_slice], axis=None))
