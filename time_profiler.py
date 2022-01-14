import time

class Timer:
    def __init__(self):
        self.timeli = []
        self.timers = []
        
    def start_timer(self, entry):
        if len(self.timeli) <= entry:
            print("New timer: ", entry)
            self.timeli = self.timeli + [0]*(entry - len(self.timeli) + 1)
            self.timers = self.timers + [0]*(entry - len(self.timers) + 1)
        try:
            self.timers[entry] = time.time()
        except:
            pdb.set_trace()
        
    def end_timer(self, entry):
        self.timeli[entry] += time.time() - self.timers[entry]

    def print_results(self):
        print(self.timeli)

