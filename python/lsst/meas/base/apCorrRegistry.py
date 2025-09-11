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

"""Registry of instrument flux fields that should be aperture corrected.
"""

__all__ = ("addApCorrName", "getApCorrNameSet")

# Set of names of algorithms that measure instrument flux that can be aperture corrected
_ApCorrNameSet = set()


def addApCorrName(name):
    """Register an instrumental flux field name prefix for aperture correction.

    Parameters
    ----------
    name : `str`
        Field name prefix for an instrumental flux that should be aperture
        corrected.

    Notes
    -----
    The prefix ``name`` corresponds to the fields ``name_instFlux``,
    ``name_instFluxErr`` and ``name_flag``. For example, specifying
    ``base_PsfFlux`` will select the fields ``base_PsfFlux_instFlux``,
    ``base_PsfFlux_instFluxErr`` and ``base_PsfFlux_flag``.
    """
    _ApCorrNameSet.add(str(name))


def getApCorrNameSet():
    """Get a copy of the field name prefixes which will be aperture corrected.

    Returns
    -------
    apCorrNameSet : `set`
        Field prefixes which will be aperture corrected.

    Notes
    -----
    For example, the returned set may include ``base_PsfFlux`` and
    ``base_GaussianFlux``.
    """
    return _ApCorrNameSet.copy()
