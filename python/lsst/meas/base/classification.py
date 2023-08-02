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

"""Definition and registration of classification plugins.
"""

import numpy as np

import lsst.pex.config
from .catalogCalculation import CatalogCalculationPluginConfig, CatalogCalculationPlugin
from .pluginRegistry import register

__all__ = (
    "CatalogCalculationClassificationConfig", "CatalogCalculationClassificationPlugin",
)


class CatalogCalculationClassificationConfig(CatalogCalculationPluginConfig):
    """Configuration for catalog classification plugin.
    """

    fluxRatio = lsst.pex.config.Field(dtype=float, default=.95, optional=True,
                                      doc="critical ratio of model to psf flux")
    modelErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                           doc="correction factor for modelFlux error")
    psfErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                         doc="correction factor for psfFlux error")


@register("base_ClassificationExtendedness")
class CatalogCalculationClassificationPlugin(CatalogCalculationPlugin):
    """Plugin which calculates a binary measure of source extendedness.

    Extendedness is based on a simple cut of the ratio of the PSF flux to the
    model flux.

    Notes
    -----
    Because the fluxes on which this algorithm is based on are slot
    measurements, they can be provided by different algorithms, and the
    `~CatalogCalculationClassificationConfig.fluxRatio` threshold used by this
    algorithm should generally be set differently for different algorithms.
    To do this, plot the difference between the PSF magnitude and the model
    magnitude vs. the PSF magnitude, and look for where the cloud of galaxies
    begins.
    """

    ConfigClass = CatalogCalculationClassificationConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def __init__(self, config, name, schema, metadata):
        CatalogCalculationPlugin.__init__(self, config, name, schema, metadata)
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for extended sources, 0 for point sources.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def calculate(self, measRecord):
        modelFlux = measRecord.getModelInstFlux()
        psfFlux = measRecord.getPsfInstFlux()
        modelFluxFlag = (measRecord.getModelFluxFlag()
                         if measRecord.table.getModelFluxSlot().isValid()
                         else False)
        psfFluxFlag = (measRecord.getPsfFluxFlag()
                       if measRecord.table.getPsfFluxSlot().isValid()
                       else False)
        flux1 = self.config.fluxRatio*modelFlux
        if self.config.modelErrFactor != 0:
            flux1 += self.config.modelErrFactor*measRecord.getModelInstFluxErr()
        flux2 = psfFlux
        if not self.config.psfErrFactor == 0:
            flux2 += self.config.psfErrFactor*measRecord.getPsfInstFluxErr()

        # A generic failure occurs when either FluxFlag is set to True
        # A generic failure also occurs if either calculated flux value is NaN:
        #     this can occur if the Flux field itself is NaN,
        #     or the ErrFactor != 0 and the FluxErr is NaN
        if np.isnan(flux1) or np.isnan(flux2) or modelFluxFlag or psfFluxFlag:
            self.fail(measRecord)
        else:
            measRecord.set(self.keyProbability, 0.0 if flux1 < flux2 else 1.0)

    def fail(self, measRecord, error=None):
        measRecord.set(self.keyFlag, True)
