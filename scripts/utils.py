class Trigger:

    def trigger(self, wildcards):
        return self.swf(self.output)

    def __init__(self, swf, output):
        self.swf = swf
        self.output = output
