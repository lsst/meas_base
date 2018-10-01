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
"""Registry for names of instrument flux fields that should be aperture corrected
"""

__all__ = ("addApCorrName", "getApCorrNameSet")

# Set of names of algorithms that measure instrument flux that can be aperture corrected
_ApCorrNameSet = set()


def addApCorrName(name):
    """Add to the set of field name prefixes for instrument flux that should be aperture corrected

    Parameters
    ----------
    name: `str`
        Field name prefix for a flux that should be aperture corrected.
        The corresponding field names are {name}_flux, {name}_fluxErr and {name}_flag.
        For example name "base_PsfFlux" corresponds to fields base_PsfFlux_flux,
        base_PsfFlux_fluxErr and base_PsfFlux_flag.
    """
    global _ApCorrNameSet
    _ApCorrNameSet.add(str(name))


def getApCorrNameSet():
    """Return a copy of the set of field name prefixes for instrument flux that should be aperture corrected.

    Notes
    -----
    For example the returned set will likely include "base_PsfFlux" and "base_GaussianFlux".
    """
    global _ApCorrNameSet
    return _ApCorrNameSet.copy()
