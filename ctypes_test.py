import ctypes as ct
import numpy as np

mandelbrot = ct.CDLL("./iteration_subroutines.so")

render = mandelbrot.render_ssaa

render.argtypes = [np.ctypeslib.ndpointer(ct.c_double, flags = "F_CONTIGUOUS"),\
                   np.ctypeslib.ndpointer(ct.c_double, flags = "F_CONTIGUOUS"),\
                   ct.c_int, ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_double,\
                   np.ctypeslib.ndpointer(ct.c_double, flags = "F_CONTIGUOUS")]

m = 1080
n = int(m*3/2)

re = np.linspace(-3,1,n)
deltar = ct.c_double(re[1] - re[0])
re = np.asfortranarray(re)
#re = re.ctypes.data_as(ct.POINTER(ct.c_double))

im = np.linspace(-1.333333333333,1.333333333333,m)
deltai = ct.c_double(im[1] - im[0])
im = np.asfortranarray(im)
#im = im.ctypes.data_as(ct.POINTER(ct.c_double))

maxiters = ct.c_int(255)
depth = ct.c_int(255)
ssaa = ct.c_int(3)
gamma = ct.c_double(1)

result = np.asfortranarray(np.zeros((n,m)))
#result = result.ctypes.data_as(ct.POINTER(ct.c_double))

n = ct.c_int(n)
m = ct.c_int(m)

result = render(re, im, maxiters, depth, ssaa, deltar, deltai, gamma, result)
