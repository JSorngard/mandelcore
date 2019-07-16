**Multi threaded mandelbrot generator**

Generates a plot of the mandelbrot set using all available cpu cores. Does the iterating and colouring in Fortran using OpenMP.

Compile mandel_subroutine.f90 as:
f2py -c mandel_subroutine.f90 -m mandelfortran --f90flags="-fopenmp" -lgomp

and image_subroutine.f90 as:
f2py -c image_subroutine.f90 -m imagefortran --f90flags="-fopenmp" -lgomp
