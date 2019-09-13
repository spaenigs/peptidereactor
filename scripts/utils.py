class Trigger:

    def trigger(self, wildcards):
        return self.swf(self.files)

    def __init__(self, swf, files):
        self.swf = swf
        self.files = files
