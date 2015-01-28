class TransformRegistry(object):
    def __init__(self):
        self._contents = {}

    def register(self, name, obj):
        self._contents[name] = obj

    def get(self, name):
        return self._contents[name] if name in self._contents else None

class TransformPlugin(object):
    def __init__(self, name, mapper, cfg, wcs, calib):
        self.name = name
        self.cfg = cfg
        self.wcs = wcs
        self.calib = calib

    def __call__(self, oldRecord, newRecord):
        raise NotImplementedError()

class PassThrough(TransformPlugin):
    def __init__(self, name, mapper, cfg, wcs, calib):
        TransformPlugin.__init__(self, name, mapper, cfg, wcs, calib)
        for key, field in mapper.getInputSchema().extract(name + "*").itervalues():
            mapper.addMapping(key)

    def __call__(self, oldRecord, newRecord):
        pass

class ReverseCentroid(PassThrough):
    def __init__(self, name, mapper, cfg, wcs, calib):
        PassThrough.__init__(self, name, mapper, cfg, wcs, calib)
        newSchema = mapper.editOutputSchema()
        self.keyRevX = newSchema.addField(self.name + "_revX", type="D", doc="reversed centroid", units="pixels")
        self.keyRevY = newSchema.addField(self.name + "_revY", type="D", doc="reversed centroid", units="pixels")

    def __call__(self, oldRecord, newRecord):
        newRecord.set(self.keyRevX, -1.0 * oldRecord['base_PeakCentroid_x'])
        newRecord.set(self.keyRevY, -1.0 * oldRecord['base_PeakCentroid_y'])

registry = TransformRegistry()
registry.register("base_PeakCentroid", ReverseCentroid)
