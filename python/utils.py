import sys, os

# Aquí tenemos algunas funciones útiles

def printerr(line):
    sys.stderr.write(u'%s\n' % line)

def suppress_stderr(func):
    def wrapper(*args, **kwargs):
        sys.stderr = open(os.devnull, "w")
        ret = func(*args, **kwargs)
        sys.stderr = sys.__stderr__
        return ret
    return wrapper
