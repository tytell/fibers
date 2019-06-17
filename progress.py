'''
Created on Apr 18, 2012

@author: eric
'''

import time, sys, math
from IPython.core.display import clear_output


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        choice = input(question + prompt).lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


class ProgressBar:
    '''
    Command line progress display.
    From http://code.activestate.com/recipes/473899-progress-meter/
    
    Here is a silly example of its usage:
    
    import progress
    import time
    import random
    
    total = 1000
    p = progress.ProgressMeter(total=total)
    
    while total > 0:
        cnt = random.randint(1, 25)
        p.update(cnt)
        total -= cnt
        time.sleep(random.random())
    
    
    Here is an example of its output:
    
    [------------------------->                                   ] 41%  821.2/sec    
    '''
    
    def __init__(self, total, refreshrate=1, ticks=60):
        self.total = total
        self.count = 0
        self.remaining = total
        if self.total > 0:
            self.meter_ticks = int(ticks)
            self.meter_division = float(self.total) / self.meter_ticks
            self.meter_value = int(self.count / self.meter_division)
        self.refreshrate = refreshrate
        self.abort = False
        self.last_update = None
        self.last_refresh = 0
        
    def __enter__(self):
        self.last_update = None
        self.count = 0
        self.start_time = time.time()
        self.elapsed = 0
        return self
    
    def __exit__(self, type, value, traceback):
        if (self.count < self.total):
            self.abort = True
        self.refresh()
        
    def update(self, count):
        now = time.time()

        self.elapsed = now - self.start_time
        self.count += count

        self.remaining = self.elapsed / self.count * (self.total - self.count)
        
        # Device Total by meter division
        value = int(self.count / self.meter_division)
        if value > self.meter_value:
            self.meter_value = value
            
        if self.last_refresh:
            if (now - self.last_refresh) > self.refreshrate or \
                (self.count >= self.total):
                    self.refresh()
        else:
            self.refresh()
            
        return self

    def get_meter(self):
        bar = '-' * self.meter_value
        pad = ' ' * (self.meter_ticks - self.meter_value)
        perc = (float(self.count) / self.total) * 100
        return '[%s>%s] %d%%  %s elapsed, %s remaining' \
          % (bar, pad, perc, self.format_sec(self.elapsed), self.format_sec(self.remaining))

    def format_sec(self, sec):
        if sec > 3600:
            minutes = int(sec/60)
            hours = int(minutes/60)
            minutes = minutes - hours*60
            return '{:d}h {:d}m'.format(hours, minutes)
        elif sec > 60:
            minutes = int(sec/60)
            sec = int(sec) - minutes*60
            return '{:d}m {:d}sec'.format(minutes, sec)
        else:
            return '{:.0f}sec'.format(sec)
        
    def refresh(self):
        if self.total == 0:
            return
        
        # Clear line
        try:
            clear_output()
        except Exception:
            pass
        print('\r', self.get_meter(), sys.stdout.flush())

        # Timestamp
        self.last_refresh = time.time()
        
