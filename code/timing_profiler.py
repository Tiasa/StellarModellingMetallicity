import time

def timing(f):
    def wrap(*args):
        t_i = time.time()
        ret = f(*args)
        t_f = time.time()
        print '[timing_profiler] {func} function took {t:.3f} s'.format(func = f.func_name, t = (t_f-t_i))
        return ret
    return wrap