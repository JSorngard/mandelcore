**Multi threaded mandelbrot generator**

Generates a plot of the mandelbrot set using all available cpu cores. Does the iterating in Fortran using OpenMP.

Compile mandel_subroutine.f90 as:
f2py -c mandel_subroutine.f90 -m mandelfortran --f90flags="-fopenmp" -lgomp

**Future:**

Would like to make the image colored, like https://preshing.com/20110926/high-resolution-mandelbrot-in-obfuscated-python/.  

Would like to not distort the image just because it is not a square.
