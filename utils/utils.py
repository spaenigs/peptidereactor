class Trigger:

    def set_input(self):
        pass

    def trigger(self, wildcards):
        # print(self.swf)
        # print(self.swf(self.output))
        return self.swf(self.output)

    def __init__(self, swf, output):
        self.swf = swf
        self.output = output
