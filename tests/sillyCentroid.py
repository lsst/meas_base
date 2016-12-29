from __future__ import absolute_import

from lsst.meas.base import SimpleAlgorithm, wrapSimpleAlgorithm, wrapTransform, BasePlugin
from _sillyCentroid import SillyCentroidAlgorithm, SillyCentroidControl, SillyTransform

# Do not register SillyCentroid in plugins.py, as it's not part of meas_base
wrapSimpleAlgorithm(SillyCentroidAlgorithm, Control=SillyCentroidControl,
                    TransformClass=SillyTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapTransform(SillyTransform)
