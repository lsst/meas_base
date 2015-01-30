from .baseLib import NaiveCentroidTransform

class TransformPlugin(object):
    def __init__(self, name, mapper, cfg):
        self.name = name
        self.cfg = cfg

    def __call__(self, oldCatalog, newCatalog, wcs, calib):
        raise NotImplementedError()

class NullTransform(TransformPlugin):
    def __call__(self, oldCatalog, newCatalog, wcs, calib):
        pass

class PassThrough(TransformPlugin):
    def __init__(self, name, mapper, cfg):
        TransformPlugin.__init__(self, name, mapper, cfg)
        for key, field in mapper.getInputSchema().extract(name + "*").itervalues():
            mapper.addMapping(key)

    def __call__(self, oldCatalog, newCatalog, wcs, calib):
        pass

class ReverseCentroid(PassThrough):
    def __init__(self, name, mapper, cfg):
        PassThrough.__init__(self, name, mapper, cfg)
        newSchema = mapper.editOutputSchema()
        self.keyRevX = newSchema.addField(self.name + "_revX", type="D", doc="reversed centroid", units="pixels")
        self.keyRevY = newSchema.addField(self.name + "_revY", type="D", doc="reversed centroid", units="pixels")

    def __call__(self, oldCatalog, newCatalog, wcs, calib):
        oldColumns = oldCatalog.getColumnView()
        newColumns = newCatalog.getColumnView()
        newColumns[self.keyRevX] = -1.0 * oldColumns['base_PeakCentroid_x']
        newColumns[self.keyRevY] = -1.0 * oldColumns['base_PeakCentroid_y']

class TransformWrapper(object):
    """
    Contains a C++ transform plugin and converts the Python config into a C++
    control object when called.
    """
    def __init__(self, transform):
        self.transform = transform

    def __call__(self, name, mapper, cfg):
        return self.transform(name, mapper, cfg.makeControl())
