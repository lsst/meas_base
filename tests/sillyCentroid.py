from __future__ import absolute_import

from lsst.meas.base import SimpleAlgorithm, wrapSimpleAlgorithm, BasePlugin
from _sillyCentroid import SillyCentroidAlgorithm, SillyCentroidControl, SillyTransform

# Do not register SillyCentroid in plugins.py, as it's not part of meas_base
wrapSimpleAlgorithm(SillyCentroidAlgorithm, Control=SillyCentroidControl,
                    TransformClass=SillyTransform, executionOrder=BasePlugin.CENTROID_ORDER)

# Emulate Swig's %feature("pythonprepend")
SillyTransform.__init = SillyTransform.__init__
def _init(self, ctrl, name, mapper):
    if hasattr(ctrl, "makeControl"):
        ctrl = ctrl.makeControl()
    self.__init(ctrl, name, mapper)
SillyTransform.__init__ = _init
