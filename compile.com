#!/bin/bash
f2py -c image_subroutine.f90 -m imagefortran --f90flags="-fopenmp" -lgomp
f2py -c mandel_subroutine.f90 -m mandelfortran --f90flags="-fopenmp" -lgomp
