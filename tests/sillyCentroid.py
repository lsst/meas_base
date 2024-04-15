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

from lsst.meas.base import wrapSimpleAlgorithm, wrapTransform, BasePlugin
from _sillyCentroid import SillyCentroidAlgorithm, SillyCentroidControl, SillyTransform

# Do not register SillyCentroid in plugins.py, as it's not part of meas_base
wrapSimpleAlgorithm(SillyCentroidAlgorithm, Control=SillyCentroidControl, name="testLib_SillyCentroid",
                    TransformClass=SillyTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapTransform(SillyTransform)
