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

import lsst.pex.config
from .transforms import PassThroughTransform

__all__ = ("BasePluginConfig", "BasePlugin")


class BasePluginConfig(lsst.pex.config.Config):
    """
    Base class measurement plugin config classes.

    Notes
    -----
    Most derived classes will want to set defaults that make sense for the
    plugin type.
    """
    pass


class BasePlugin:
    """
    Base class for measurement plugins.

    This is the base class for `SingleFramePlugin` and `ForcedPlugin`; derived
    classes should inherit from one of those.

    Parameters
    ----------
    config : `BasePluginConfig`
        Plugin configuration.
    name : `str`
        Plugin name.
    logName : `str`
        Logger name.

    Notes
    -----
    Relative execution orders are defined by a series of named constants
    defined in this class: plugins with a lower execution number are run
    first.

    This approach was chosen instead of a full graph-based analysis of
    dependencies because algorithm dependencies are usually both quite simple
    and entirely substitutable: an algorithm that requires a centroid can
    typically make use of any centroid algorithms outputs.  That makes it
    relatively easy to figure out the correct value to use for any particular
    algorithm.
    """

    CENTROID_ORDER = 0.0
    """Order for algorithms which require only Footprint and Peaks (`float`).

    Notes
    -----
    Algorithms with this execution order include centroids.
    """

    SHAPE_ORDER = 1.0
    """Order for algorithms which require a centroid (`float`).

    Notes
    -----
    These algorithms may refer assume that `getCentroid` will return a good
    centroid, and that a Footprint and its Peaks are available.
    """

    FLUX_ORDER = 2.0
    """Order for algorithms which require a shape and a centroid (`float`).

    Notes
    -----
    These algorithms may assume that both `getCentroid` and `getShape` will
    return good values, and that a Footprint and its Peaks are available.
    """

    APCORR_ORDER = 3.0
    """Order for algorithms which require shape, centroid and flux (`float`).

    Notes
    -----
    These algorithms may assume that `getCentroid` and `getShape` will return
    good values, that flux has been measured, and that and that a Footprint
    and its Peaks are available.
    """

    DEFAULT_CATALOGCALCULATION = 4.0
    """Order for catalog calculation plugins.

    Notes
    -----
    These plugins only operate on catalogs; they may not access pixel values.
    """

    ConfigClass = BasePluginConfig
    """Plugin configuration information (`lsst.pex.config.Config`).
    """

    @classmethod
    def getExecutionOrder(cls):
        """Get the relative execution order of this plugin.

        Must be reimplemented as a class method by concrete derived classes.
        """
        raise NotImplementedError("All plugins must implement getExecutionOrder()")

    def __init__(self, config, name, logName=None):
        object.__init__(self)
        self.config = config
        self.name = name
        self.logName = logName

    def getLogName(self):
        return self.logName

    def fail(self, measRecord, error=None):
        """Record a failure of the `measure` or `measureN` method.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Table record describing the source being measured.
        error : `MeasurementError`, optional
            Only provided if the measurement failed due to a
            `MeasurementError` being raised; otherwise, will be `None`.

        Notes
        -----
        When the plugin raises an exception, framework will call
        `BasePlugin.fail` to allow the plugin to set its failure flag
        field(s).  When `BasePlugin.measureN` raises an exception,
        `BasePlugin.fail` will be called repeatedly with all the records that
        were being measured.

        If the exception is an `MeasurementError`, it will be passed as the
        error argument; in all other cases the error argument will be `None`,
        and the failure will be logged by the measurement framework as a
        warning.

        """
        traceback.print_exc()
        message = ("The algorithm '%s' thinks it cannot fail, but it did; "
                   "please report this as a bug (the full traceback is above)."
                   % (self.__class__.__name__,))
        raise NotImplementedError(message)

    @staticmethod
    def getTransformClass():
        """Get the measurement transformation appropriate to this plugin.

        This returns a subclass of `transforms.MeasurementTransform`, which
        may be instantiated with details of the algorithm configuration and
        then called with information about calibration and WCS to convert from
        raw measurement quantities to calibrated units. Calibrated data is
        then provided in a separate output table.

        Notes
        -----
        By default, we copy everything from the input to the output without
        transformation.
        """
        return PassThroughTransform
