class Costam:

    class Config(dict):
        def __init__(self, **kwargs):
            self._dict = {}
            self._dict.update(kwargs)

    def __init__(self, **kwargs):
        self.config = self.Config(**kwargs)

    def __repr__(self):
        return self.config

print Costam()