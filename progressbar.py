class progressbar:
    "A progress bar for very long operations."

    def __init__(self):
        self.start = False
        self.length = 0

    def display(self, message):
        if self.start:
            pad = self.length - len(message)
            if pad < 0:
                pad = 0
            print("\r" + message + " " * pad, end="")
        else:
            print(message, end="")
        self.start = True
        self.length = len(message)

    def done(self):
        print("")

