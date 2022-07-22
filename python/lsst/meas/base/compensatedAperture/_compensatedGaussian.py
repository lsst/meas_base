from __future__ import annotations

__all__ = (
    "CompensatedGaussianAperturePluginConfig",
    "CompensatedGaussianAperturePlugin",
)

import math
import numpy as np
import scipy.stats as sps

from lsst.pex.config import Field
from lsst.geom import Point2I

from ..pluginRegistry import register

from ._compensatedAperture import (
    CompensatedAperturePluginConfig,
    CompensatedAperturePlugin,
    OutOfBoundsException,
)
from ._calcGaussian import _gaussianFiltInnerProduct


class CompensatedGaussianAperturePluginConfig(CompensatedAperturePluginConfig):
    t = Field(
        doc="Scale parameter of outer Gaussian compaired to inner",
        dtype=float,
        default=2,
    )
    alpha = Field(doc="Approximate width of the psf if it were a symmetric gaussian", dtype=float, default=3)


@register("base_CompensatedGaussianAp")
class CompensatedGaussianAperturePlugin(CompensatedAperturePlugin):
    ConfigClass = CompensatedGaussianAperturePluginConfig

    def __init__(
        self,
        config: CompensatedGaussianAperturePluginConfig,
        name: str,
        schema,
        metadata,
        logName=None,
        **kwds,
    ):
        super().__init__(config, name, schema, metadata, logName, **kwds)

        self.width_keys = {}
        self._rads = {}
        self._flux_corrections = {}
        self._variance_corrections = {}
        self._t = config.t
        for width in config.kernel_widths:
            base_key = f"{name}_{width}"

            # flux
            flux_str = f"{base_key}_instFlux"
            flux_key = schema.addField(flux_str, type="D", doc="Compensated flux measurement")

            # uncertainty
            uncert_str = f"{base_key}_instFluxErr"
            uncert_key = schema.addField(uncert_str, type="D", doc="Compensated flux uncertainty")

            # mask bits
            mask_str = f"{base_key}_mask_bits"
            mask_key = schema.addField(mask_str, type=np.int32, doc="Mask bits set within Aperture")

            # normalization, set to 1.0 here and modified later
            norm_str = f"{base_key}_normalization"
            norm_key = schema.addField(norm_str, type="D", doc="Compensated flux normalization; set later")

            self.width_keys[width] = (flux_key, uncert_key, mask_key, norm_key)
            self._rads[width] = math.ceil(sps.norm.ppf((0.995,), scale=width * config.t)[0])

            alpha = self.config.alpha
            self._flux_corrections[width] = (
                2
                * np.pi
                * (alpha**4 + alpha**2 * width**2 * (self._t**2 + 1) + width**4 * self._t**2)
                / (width**2 * (self._t**2 - 1))
            )
            self._variance_corrections[width] = 4 * np.pi * width**2 * (self._t**2 + 1) / self._t**2

        self._max_rad = max(self._rads)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()

        if Point2I(center) not in exposure.getBBox().erodedBy(self._max_rad):
            raise OutOfBoundsException("Not all the kernels for this source fit inside the detector")

        y = center.getY()
        x = center.getX()

        y_floor = math.floor(y)
        x_floor = math.floor(x)

        for width, (flux_key, uncert_key, mask_key, norm_key) in self.width_keys.items():
            rad = self._rads[width]
            y_slice = slice(y_floor - rad, y_floor + rad + 1, 1)
            x_slice = slice(x_floor - rad, x_floor + rad + 1, 1)
            y_mean = y - y_floor + rad
            x_mean = x - x_floor + rad
            flux_uncal, uncert_uncal = _gaussianFiltInnerProduct(
                exposure.image.array[y_slice, x_slice], x_mean, y_mean, width, self._t
            )
            measRecord.set(flux_key, flux_uncal * self._flux_corrections[width])
            measRecord.set(uncert_key, np.sqrt(uncert_uncal * self._variance_corrections[width]))
            measRecord.set(mask_key, np.bitwise_or.reduce(exposure.mask.array, axis=None))
            measRecord.set(norm_key, self._flux_corrections[width])
