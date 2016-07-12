from lsst.meas.base import FlagDefinition, FlagDefinitionVector, FlagHandler


def addFlagHandler(*args):
    '''
    Class decorator to create a flag handler for a plugin. Adds the class variables FLAGDEFS and
    ErrEnum. An instnace variable flagHandler is added to the __init__ function to be created at
    initialization. The arguments to this function are tuples which have the name of the failure
    as the first element, and the documentation for the failure as the second element.

    Usage:
    @addFlagHandler(("name_of_failure", "Doc for failure), ("name_of_second_faiure", "Doc of"
                    " second failure"), .....)
    '''
    def classFactory(cls):
        cls.FLAGDEFS = [FlagDefinition(name, doc) for name, doc in args]
        # Verify all flag names are unique
        names = [entry[0] for entry in args]
        if len(names) != len(set(names)):
            raise ValueError("All flag names must be unique, given {}".format(names))
        # Better to have a scoped enumeration rather than attach variables strait to the class to
        # prevent shadowing
        cls.ErrEnum = type('ErrEnum', (), {entry.name: i for i, entry in enumerate(cls.FLAGDEFS)})
        oldInit = cls.__init__

        def newInit(self, *args, **kwargs):
            if 'schema' in kwargs:
                schema = kwargs['schema']
            else:
                schema = args[2]
            if 'name' in kwargs:
                name = kwargs['name']
            else:
                name = args[1]
            self.flagHandler = FlagHandler.addFields(schema, name, FlagDefinitionVector(self.FLAGDEFS))
            oldInit(self, *args, **kwargs)
        cls.__init__ = newInit
        return cls
    return classFactory
