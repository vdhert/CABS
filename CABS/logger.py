import os, sys
from threading import Thread
from threading import Event
from subprocess import check_output
from sys import stderr
from os.path import exists
from time import time, strftime, gmtime, sleep
import textwrap,trace

_name = "Logger"

colors = {
    "blue": '\033[94m',
    "yellow": '\033[93m',
    'green': '\033[92m',
    'red': '\033[91m',
    'light_blue': '\033[96m',
    'purple': '\033[95m',
    'end': '\033[0m',
}

log_levels = {

   -1: "CRITICAL",
    0: "WARNING",
    1: "INFO",
    2: "OUT FILES",
    3: "DEBUG"
}

color_prefix = {
   -1: colors["red"] + "[CRITICAL]" + colors["end"],
    0: colors["yellow"] + "[WARNING]" + colors["end"],
    1: colors["blue"] + "[INFO]" + colors["end"],
    2: colors["purple"] + "[OUT FILES]" + colors["end"],
    3: colors["green"] + "[DEBUG]" + colors["end"],
}

_init_time = time()
log_level = 1
color = sys.stderr.isatty()
stream = sys.stderr
if color: prefix = color_prefix
else: prefix = log_levels


def setup_log_level(new_level):
    global log_level
    if type(new_level) is int and new_level < 4 and new_level >= -1:
        log_level = new_level
    else:
        warning(module_name=_name,msg="Verbose should be a number between -1 and 3")
    info(module_name=_name ,msg="Verbosity set to: " + str(log_level) + " - " + log_levels[log_level])


def coloring(color_name = "light_blue", msg = ""):
    if color:
        return colors[color_name] + msg + colors["end"]
    return msg


def log(module_name = "MISC", msg = "Processing ", l_level = 1, out = stream):
    if l_level <= log_level:
        t = gmtime(time() - _init_time)
        if len(msg) < 76:
            msg = '%-20s %-19s%-75s %s\n' % (
                prefix[l_level] , coloring(msg=module_name+":",color_name='light_blue') , msg,strftime('(%H:%M:%S)', t))
            out.write(msg)
            out.flush()
        else:
            lines = textwrap.wrap(msg,width=75)
            firstLine = '%-20s %-19s%-75s \n' % (prefix[l_level] , coloring(msg=module_name+":",color_name='light_blue') , lines[0])
            out.write(firstLine)
            for lineNumber in xrange(1,len(lines)-1):
                line = '%-22s%-75s \n' % (" ", lines[lineNumber])
                out.write(line)
            finalLine = '%-22s%-75s %s\n' % (" ", lines[-1], strftime('(%H:%M:%S)', t))
            out.write(finalLine)
            out.flush()


def critical(module_name = "_name", msg = ""):
    log(module_name=module_name, msg=msg, l_level=-1)


def warning(module_name = "_name", msg = ""):
    log(module_name=module_name, msg=msg, l_level=0)


def info(module_name = "_name", msg = ""):
    log(module_name=module_name, msg=msg, l_level=1)


def log_file(module_name = "_name", msg = ""):
    log(module_name=module_name, msg=msg, l_level=2)


def debug(module_name = "_name", msg = ""):
    log(module_name=module_name, msg=msg, l_level=3)


def to_file(filename='',content='',msg='',allowErr=True,traceback=True):
    """

    :param filename: path for the file to be saved
    :param content: a string to be saved (be careful not to pass a string that is too long
    :param msg: optional: a message to be logged
    :param allowErr: if True a warning will be logged on OSError, if False program exit call will be made
    :param traceback: if True an Exception will be raised on exit call
    :return:
    """
    if filename and content:
        try:
            if os.path.isfile(filename): log_file(module_name=_name,msg = "Overwriting %s" % filename)
            with open(filename,'w') as f:
                f.write(content)
        except IOError:
            if allowErr:
                warning(module_name=_name,msg ="IOError while writing to: %s" % filename)
            else:
                exit_program(module_name=_name,msg ="IOError while writing to: %s" % filename,traceback=traceback)
    if msg:
        log_file(module_name=_name,msg = msg)


def exit_program(module_name =_name, msg="Shutting down",traceback=True,exc=None):
    """
    Exits the program, depending on options might do it quietly, raise, or print traceback
    :param module_name: string with  the calling module's name
    :param msg: string, message to be printed when the program exits
    :param traceback: bool, if log level is high traceback will be printed
    :param exc: a specific exception passed by the caller
    :return: None
    """
    critical(module_name=module_name,msg=msg)
    if log_level > 2 and traceback:
        if exc:
            debug(module_name=module_name, msg="Exception caught. Printing traceback:")
            for line in trace.format_stack():
                stream.write(line)
        else:
            raise
    sys.exit(1)


class ProgressBar:
    ''' This class assumes a manual call to done() will be made to exit the bar'''
    WIDTH = 65
    FORMAT = '%-20s %-19s[%s] %.1f%%\r'
    BAR0 = ' '
    BAR1 = '#'

    def __init__(self, total=100, module_name='', job_name='', stream=stderr, delay=0, start_msg=''):
        if log_level >=1 :
            self.stream = stream
        else: self.stream = open('/dev/null', 'w')
        self.total = total
        self.current = 0.
        self.job_name = job_name
        self.is_done = False
        self.module_name = module_name
        self.prefix = prefix[1]
        if start_msg:
            self.stream.write(coloring(msg=start_msg) + '\n')
        if self.job_name:
            log(module_name=self.module_name, msg=self.job_name + " running..." ,out=self.stream)
        self.start_time = time()
        self.update(state=0.)
        sleep(delay)

    def write(self):
        percent = 1.0 * self.current / self.total
        num = int(self.WIDTH * percent)
        percent = round(100. * percent, 1)
        bar = self.BAR1 * num + self.BAR0 * (self.WIDTH - num)
        self.stream.write(self.FORMAT % (self.prefix,coloring(msg=self.module_name+":"),bar, percent))
        self.stream.flush()

    def update(self, state=False):
        if state is False:
            self.current += 1
        else:
            self.current = state
        if self.current >= self.total:
            return True
        self.write()
        return False

    def finish(self):
        if self.current < self.total:
            for i in xrange(int(self.total - self.current)):
                self.update()
                sleep(0.001)

    def done(self, show_time=True):
        if not self.is_done:
            self.finish()
            self.stream.write(" "*80+"\r")
            if show_time:
                t = gmtime(time() - self.start_time)
                log(module_name=self.module_name, msg=self.job_name + ' done in %s' % strftime('%H:%M:%S', t),out=self.stream)
            self.stream.flush()
            self.is_done = True


class cabs_observer(Thread):

    def __init__ (self,interval=0.5,traj='',n_lines = 0, job_name ='CABS simulation',msg =''):
        Thread.__init__(self)
        self.exit_event=Event()
        self.interval=interval
        self.progress_bar = ProgressBar(module_name='CABS',job_name=job_name,start_msg=msg)
        self.traj = traj
        self.n_lines = n_lines
        self.daemon = True  # In case main program ends abruptly
        self.start()

    def exit(self):
        self.progress_bar.done()
        self.exit_event.set()

    def run(self):
        while not self.exit_event.isSet():
            if self.progress_bar.update(self.status()):
                self.exit()
            sleep(self.interval)

    def status(self):
        if not exists(self.traj):
            progress = 0.
        else:
            progress = 100. * int(check_output(['wc', '-l', self.traj]).split()[0]) / self.n_lines
        return progress


