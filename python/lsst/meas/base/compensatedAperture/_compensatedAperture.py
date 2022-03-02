__all__ = (
    "OutOfBoundsException",
    "CompensatedAperturePluginConfig",
    "CompensatedAperturePlugin",
)

from lsst.pex.config import ListField, FieldValidationError

from ..sfm import SingleFramePlugin, SingleFramePluginConfig


class OutOfBoundsException(Exception):
    pass


class CompensatedAperturePluginConfig(SingleFramePluginConfig):
    kernel_widths = ListField(
        doc="The widths of the kernels for which to " "measure compensated apertures",
        dtype=int,
        minLength=1,
        default=[1, 2, 3, 4, 7, 11, 17, 23, 27, 50]
    )

    # commented out until I can create a ticket, since this actually breaks
    # validation for pex_config
    #def validate(self, instance):
        #if len(self.kernel_widths) != len(set(self)):
            #raise FieldValidationError(
                #type(self).kernel_widths, self, "All values for kernel_width must be unique"
            #)


class CompensatedAperturePlugin(SingleFramePlugin):
    ConfigClass = CompensatedAperturePluginConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata, logName=None, **kwds):
        super().__init__(config, name, schema, metadata)

        # create generic failure key
        self.fatalFailKey = schema.addField(
            f"{name}_flag", type="Flag", doc="Set to 1 for any fatal failure"
        )

        # Out of bounds failure key
        self.ooBoundsKey = schema.addField(
            f"{name}_bounds_flag",
            type="Flag",
            doc="Flag set to 1 if now all filters fit within exposure",
        )

    def fail(self, measRecord, error=None):
        if isinstance(error, OutOfBoundsException):
            measRecord.set(self.ooBoundsKey, True)
        measRecord.set(self.fatalFailKey, True)
