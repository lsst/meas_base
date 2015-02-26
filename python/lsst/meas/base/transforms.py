#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
"""Measurement transformations.

When a measurement plugin is run, it provides raw, uncalibrated outputs such
as pixel positions. A transformation may be run as a post-processing step to
convert those outputs to calibrated quantities, such as celestial coordinates.

At construction, the transformation is passed the configuration and name of
the plugin whose outputs it will be transformaing (all fields in the input
table produced by that plugin will have their field names prefixed by the
plugin name) and a `SchemaMapper` which holds the schemata for the input and
output catalogs and which may be used to directly map fields between the
catalogs.

When a transformer is called, it is handed a `SourceCatalog` containing the
measurements to be transformed, a `BaseCatalog` in which to store results, and
information about the WCS and calibration of the data. It may be safely
assumed that both are contiguous in memory, thus a ColumnView may be used for
efficient processing.

Transformations can be defined in Python or in C++. Python code should inherit
from `MeasurementTransform`, following its interface. C++ code inherits from
`BaseTransform`. `TransformWrapper` wraps the C++ implementation so that it
provides a uniform interface with Python.
"""

__all__ = ("NullTransform", "PassThroughTransform", "TransformWrapper")

class MeasurementTransform(object):
    """!
    Base class for measurement transformations.

    Create transformations by deriving from this class, implementing
    `__call__()` and (optionally) augmenting `__init__()`.
    """
    def __init__(self, config, name, mapper):
        self.name = name
        self.config = config

    def __call__(self, inputCatalog, outputCatalog, wcs, calib):
        raise NotImplementedError()


class NullTransform(MeasurementTransform):
    """!
    The null transform transfers no data from input to output.

    This is intended as the default for measurements for which no other
    transformation is specified.
    """
    def __call__(self, inputCatalog, outputCatalog, wcs, calib):
        pass


class PassThroughTransform(MeasurementTransform):
    """!
    Copy all fields named after the measurement plugin from input to output, without transformation.
    """
    def __init__(self, config, name, mapper):
        MeasurementTransform.__init__(self, config, name, mapper)
        for key, field in mapper.getInputSchema().extract(name + "*").itervalues():
            mapper.addMapping(key)

    def __call__(self, inputCatalog, outputCatalog, wcs, calib):
        pass


class TransformWrapper(object):
    """!
    Wrap a C++-style transformation.

    C++ code requires a `Control`, rather than a `Config`, object. This class
    converts the latter to the former: when the TransformWrapper is invoked
    Note that this is only possible for measurement plugins which have control
    objects defined in C++: pure Python measurement plugins must use Python
    transformations.
    """
    def __init__(self, transform):
        self.transform = transform

    def __call__(self, config, name, mapper):
        return self.transform(config.makeControl(), name, mapper)
