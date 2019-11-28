**Multi threaded mandelbrot generator**

Generates a plot of the mandelbrot set using all available cpu cores. Does the iterating and colouring in Fortran using OpenMP.
Allows the user to specify resolution, whether the image should be coloured or black and white, where to center the image and how much to zoom among other things. Also allows the generation of animated gifs of zooming in to the fractal.

Compile by running the compile.com script.

Command line arguments:  
  -cr RCENTER, --Rcenter RCENTER  
                        Specify the real part of the point in the complex  
                        plane to center the image on. Defaults to -0.75.
                        
  -ci ICENTER, --Icenter ICENTER  
                        Specify the imaginary part of the point in the complex  
                        plane to center the image on. Defaults to 0.0.
                        
  -y YRESOLUTION, --yresolution YRESOLUTION  
                        Specify the y-axis resolution of the image. Defaults  
                        to 2160.
                        
  -r ASPECTRATIO, --aspectratio ASPECTRATIO  
                        Specifies the aspect ratio of the image. Defaults to 1.5.
                        
  -z ZOOM, --zoom ZOOM  The zoom in factor to reach.  
  
  -b, --blackwhite      Use this command if you want the image to be black and  
                        white.  
                        
  -C, --colourshift     Use this flag if you want to shift the colouring so  
                        that the fastest escaping point is always blue.  
                        
  -f FRAMES, --frames FRAMES  
                        The number of frames per until the full zoom is  
                        achieved.  
                        
  -d DURATION, --duration DURATION  
                        The number of seconds the animation should take  
                        (applicable if frames and zoom > 1).  
                        
  -g GAMMA, --gamma GAMMA  
                        Raises the output of the mandelbrot iterations to this  
                        number. Works as a gamma between 0.4 and 1.  
                        
  -s SSAAFACTOR, --ssaafactor SSAAFACTOR  
                        Supersample each pixel this many times squared.  
                        Defaults to 3. If 1, no SSAA will be computed.  
                        
  -e FILEEXTENSION, --fileextension FILEEXTENSION  
                        Set the file extension of the generated image.  
                        Defaults to bmp.  
                        
  -l COMPRESSLEVEL, --compresslevel COMPRESSLEVEL  
                        Specify the compression level of the image if the file  
                        extension is set to png. Defaults to 6.  
                        
  --saveresult          Use this argument if you want to save the results of  
                        the mandelbrot iterations to a .dat.gz file.  
                        
  --noimage             Use this argument if you do not want the program to  
                        output an image file.  
                        
  -v, --verbose         Use this argument if you want more detailed  
                        information on the computation.  
