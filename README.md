**Multi threaded mandelbrot generator**

Generates a plot of the mandelbrot set using all available cpu cores. Does the iterating and colouring in Fortran using OpenMP.

Command line arguments:  
-c CENTER, --center CENTER  
    Specify the point in the complex plane to center the image on. Defaults to -75.

-y YRESOLUTION, --yresolution YRESOLUTION  
    Specify the number of pixels along the y-axis. Defaults to 1440.

-r ASPECTRATIO, --aspectratio ASPECTRATIO  
    Specifies the aspect ratio of the image. Defaults to 3/2.

-z ZOOM, --zoom ZOOM  
    The zoom factor to reach. E.g. if -z 2 is used the distance along any axis will be 1/2 times the orignal.  
    Defaults to 1.

-C, --colourshift  
    Use this flag if you want to shift the colouring so that the fastest escaping point is always blue.

-f FRAMES, --frames FRAMES  
    The number of frames until the full zoom is achieved.

-g GAMMA, --gamma GAMMA  
    Raise the output of the mandelbrot iterations to this number. Works as a gamma between .4 and 1. Defaults to 1.

-s SSAAFACTOR, --ssaafactor SSAAFACTOR  
    Supersample each pixel this many times per axis. If set to 1, no supersampling will be used. Defaults to 3.

-e FILEEXTENSION, --fileextension FILEEXTENSION  
    Set the file extension of the generated image to this (include the dot). Defaults to .bmp.  
    Ignored if frames is set to a number larger than 1, in which case it is set to .gif.

--saveresult  
    Use this argument if you want to save the result of the iterations to a .dat.gz file.

--debug  
    Use this argument if you want more detailed information on the computation.  
    Is on by default if making only one image.

---

Compile mandel_subroutine.f90 as:  
f2py -c mandel_subroutine.f90 -m mandelfortran --f90flags="-fopenmp" -lgomp

and image_subroutine.f90 as:  
f2py -c image_subroutine.f90 -m imagefortran --f90flags="-fopenmp" -lgomp
