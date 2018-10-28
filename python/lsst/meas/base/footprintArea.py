#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy as np

from .catalogCalculation import (CatalogCalculationPluginConfig,
                                 CatalogCalculationPlugin)
from .pluginRegistry import register

__all__ = (
    "CatalogCalculationFootprintAreaConfig",
    "CatalogCalculationFootprintAreaPlugin",
)


class CatalogCalculationFootprintAreaConfig(CatalogCalculationPluginConfig):
    """Configuration for footprint area catalog calculation plugin.
    """

    pass


@register("base_FootprintArea")
class CatalogCalculationFootprintAreaPlugin(CatalogCalculationPlugin):
    """Catalog calculation plugin to record the area of a source's footprint.
    """

    ConfigClass = CatalogCalculationFootprintAreaConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def __init__(self, config, name, schema, metadata):
        CatalogCalculationPlugin.__init__(self, config, name, schema, metadata)
        self.key = schema.addField(
            schema.join(name, "value"),
            type=np.int32,
            doc="Number of pixels in the source's detection footprint.",
            units="pixel"
        )

    def calculate(self, measRecord):
        measRecord.set(self.key, measRecord.getFootprint().getArea())

    def fail(self, measRecord, error=None):
        # Should be impossible for this algorithm to fail unless there is no
        # Footprint (and that's a precondition for measurement).
        pass
