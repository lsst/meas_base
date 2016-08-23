#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import traceback

from builtins import object

import lsst.pex.config
from .transforms import PassThroughTransform

__all__ = ("BasePluginConfig", "BasePlugin")


class BasePluginConfig(lsst.pex.config.Config):
    """!
    Base class measurement Plugin config classes.

    Most derived classes will want to set defaults that make sense for the plugin type
    """
    pass


class BasePlugin(object):
    """!
    Base class for measurement plugins.

    This is the base class for SingleFramePlugin and ForcedPlugin; derived classes should inherit
    from one of those.
    """
    # named class constants for execution order
    CENTROID_ORDER = 0.0
    SHAPE_ORDER = 1.0
    FLUX_ORDER = 2.0
    APCORR_ORDER = 3.0
    DEFAULT_AFTERBURNER = 4.0

    @classmethod
    def getExecutionOrder(cls):
        """Sets the relative order of plugins (smaller numbers run first).

        In general, the following class constants should be used (other values
        are also allowed, but should be avoided unless they are needed):
        CENTROID_ORDER      centroids and other algorithms that require only a Footprint and its Peaks as input
        SHAPE_ORDER         shape measurements and other algorithms that require getCentroid() to return
                            a good centroid (in addition to a Footprint and its Peaks).
        FLUX_ORDER          flux algorithms that require both getShape() and getCentroid(),
                            in addition to a Footprint and its Peaks
        DEFAULT_AFTERBURNER plugins that only operate on the catalog

        Must be reimplemented as a class method by concrete derived classes.

        This approach was chosen instead of a full graph-based analysis of dependencies
        because algorithm dependencies are usually both quite simple and entirely substitutable:
        an algorithm that requires a centroid can typically make use of any centroid algorithms
        outputs.  That makes it relatively easy to figure out the correct value to use for any
        particular algorithm.
        """
        raise NotImplementedError("All plugins must implement getExecutionOrder()")

    def __init__(self, config, name):
        """!
        Initialize the plugin object.

        @param[in]  config       An instance of this class's ConfigClass.
        @param[in]  name         The string the plugin was registered with.
        """
        object.__init__(self)
        self.config = config
        self.name = name

    def fail(self, measRecord, error=None):
        """!
        Record a failure of the measure or measureN() method.

        When the plugin raises an exception, framework will call
        fail() to allow the plugin to set its failure flag
        field(s).  When measureN() raises an exception, fail() will be
        called repeatedly with all the records that were being
        measured.

        If the exception is a MeasurementError, it will be passed as
        the error argument; in all other cases the error argument will
        be None, and the failure will be logged by the measurement
        framework as a warning.
        """
        traceback.print_exc()
        message = ("The algorithm '%s' thinks it cannot fail, but it did; "
                   "please report this as a bug (the full traceback is above)."
                   % self.__class__.__name__)
        raise NotImplementedError(message)

    @staticmethod
    def getTransformClass():
        """!
        Get the measurement transformation appropriate to this plugin.

        This returns a subclass of MeasurementTransform, which may be
        instantiated with details of the algorithm configuration and then
        called with information about calibration and WCS to convert from raw
        measurement quantities to calibrated units. Calibrated data is then
        provided in a separate output table.

        By default, we copy everything from the input to the output without
        transformation.
        """
        return PassThroughTransform
