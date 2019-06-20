import numpy as np
import multiprocessing as mp
import sys
import imageio
import os
from mandelfortran import *
#from mandelfortran import mandel_calc
#from mandelfortran import mandel_calc_array
#from mandelfortran import colorize as fcolor
import time

#The highest allowed number of iterations.
iters = 100

#Sets the aspect ratio of the image.
aspect_ratio = 3./2.

#Defines the "window" to look at the fractal in. Make sure to match the apect ratio.
start = -2.7-1.333j
end = 1.3+1.333j

#Number of points per axis to compute.
im_eval_points = 15000 #y-axis. Must be even.
re_eval_points = int(aspect_ratio*im_eval_points) #x-axis

#Compute it multithreaded.
multicore = True

#The method of parallelization.
#If true it will be openmp in fortran, otherwise multiprocessing in python.
fortran_omp = True

#Save the result as an image.
saveimage = True

#Make a pass of gaussian blur.
blur = False #Blurring requires allocating an extra image in memory, just as large as the main one. The image can therefore not be as big with this option turned on.
radius = 1 #the radius of the gaussian blur. 1 is probably optimal.

#Make the image in color. Only relevant if saveimage is True.
colorize = False #This option requires the allocation of an extra image in memory, three times as large as the main one. The image must therefore be much smaller with this option turned on.

#Save the resulting iteration grid to file
saveresult = False

#What file type to save the image as.
image_file_ext = ".png"

#What file type to save the raw data as. If it ends in .gz it will be compressed.
data_file_ext = ".dat.gz"

#Color depth.
depth = 255

#Print extra information about memory use.
memory_debug = False

def mandel_func(c,maxiterations=iters):
	"""Takes a complex number and iterates the mandelbrot function on it until either its magnitude is larger than six (2 gives worse colour fade), or it has iterated enough."""
	x = np.real(c)
	y = np.imag(c)
	y2 = y**2.
	#q=(x-.25)**2. + y2
	mag2=x**2 + y2

	#Filter out all points inside the main bulb and the period-2 bulb.
	if((x + 1.)**2 + y2 < 0.0625 or mag2*(8.*mag2-3.) <= .09375 - x):
		return maxiterations

	z = 0+0j
	iterations=0
	while(np.real(z)**2. + np.imag(z)**2. <= 36. and iterations < maxiterations):
		iterations += 1
		z = z**2. + c

	return iterations

if(__name__ == "__main__"):

	if(not saveresult and not saveimage):
		print("Note: program will produce no output.")

	#Determines the working directory of the program.
	path = os.path.dirname(os.path.abspath(__file__))

	#Determines whether to use \ or / for file paths.
	pathdelim = "\\" if sys.platform == "win32" else "/"

	#Multiplatform clock
	get_timer = time.perf_counter if sys.platform == "win32" else time.time
	
	total_time = get_timer()

	colorname = "_color" if colorize else "_bw"

	blurname = "_blur="+str(radius) if blur else ""

	#Only compute half the points along the imaginary axis since there is a reflection symmetry.
	if(np.mod(im_eval_points,2)==0):
		im_eval_points = int(im_eval_points/2)
	else:
		print("Number of imaginary points must be even.")
		exit()

	re_points= np.linspace(np.real(start),np.real(end),re_eval_points)
	im_points= np.linspace(0,np.imag(end),im_eval_points)

	re_grid,im_grid = np.meshgrid(re_points,im_points*1j,sparse=True)

	grid = re_grid + im_grid
	
	if(memory_debug):		
		gridshape = np.shape(grid)
		elements = gridshape[0]*gridshape[1]
		cmplxsize = sys.getsizeof(1+1j)
		cmplxnparraysize = sys.getsizeof(np.array(1+1j))
		cmplxnparray10size = sys.getsizeof((1+1j)*np.ones(1,dtype=complex))

		print("Size of a complex number: "+str(cmplxsize)+" B.")
		print("Size of a numpy array with a complex number: "+str(cmplxnparraysize)+" B.")
		print("Size of a numpy array with 10 complex numbers: "+str(cmplxnparray10size)+" B.")
		print("Elements in grid: "+str(elements)+".")
		print("Grid 'should' take up roughly "+str(elements*cmplxsize/1000000)+" MB.")
		print("Size of grid: "+str(sys.getsizeof(grid)/1e6)+" MB.")

	if(multicore):
		cores = mp.cpu_count()
		print("Attempting to evaluate on "+str(cores)+" cores...")
		eval_type = "_multicore"

		if(fortran_omp):
			print("Computing...")
			time = get_timer()
			grid = mandel_calc_array(grid,iters)
			result = grid
			time = get_timer() - time
			print("Done in "+str(time)[:4]+" seconds.")			

		else:
			#Create a pool with the number of threads equal to the number of processor cores.
			print("Creating thread pool...")
			time = get_timer()
			pool = mp.Pool(processes=cores)
			#Warm up the pool
			pool.map(mandel_helper,np.ones((10,cores)))
			time = get_timer() - time
			print("Done in "+str(time)[:5]+" seconds.")

			print("Computing...")
			time = get_timer()
			result = pool.map(mandel_helper,grid)
			time = get_timer() - time
			print("Done in "+str(time)[:4]+" seconds.")

	else:
		print("Evaluating on a single core...")
		eval_type = "_singlecore"
		mandel_vector = np.vectorize(mandel_func)

		print("Computing...")
		time = get_timer()
		result = mandel_vector(grid,iters)
		time = get_timer() - time
		print("Done in "+str(time)[:4]+" seconds.")

	if(saveresult):
		print("Writing raw data...")
		if(data_file_ext[-3:] == ".gz"):
			print(" compressing...")
		time = get_timer()
		#Write image to file. If the file name ends in .gz numpy automatically compresses it.
		np.savetxt(path+pathdelim+"mandel"+colorname+eval_type+data_file_ext,result,delimiter=' ')
		time = get_timer() - time
		print("Done in "+str(time)[:4]+" seconds.")

	if(saveimage):
		print("Performing image manipulations...")
		time = get_timer()
		
		if(blur):
			print(" blurring...")
			blurred = np.zeros(np.shape(result))
			blurred = fastgauss(np.real(result),radius)
			result = blurred

		#Normalize result to 0-1
		
		print(" normalizing...")
		result = np.array(result)/float(iters)
		
		

		if(colorize):
			print(" colouring...")
			colorized = np.zeros((im_eval_points,re_eval_points,3))
			colorized = fcolor(colorized,result)
			result = colorized

		else:
			print(" fitting to color depth...")
			#Scale up to 0-depth. What should be black is now depth.
			result *= depth
			#Invert so that black is 0 and white is depth.
			result -= depth 
			result = np.abs(result) 
	
		#Convert to uints for imageio.
		result = result.astype(np.uint8)
		print(" mirroring...")
		#Adds a flipped copy of the image to the top.
		result = np.concatenate((np.flip(result,axis=0),result))

		time = get_timer() - time
		print("Done in "+str(time)[:4]+" seconds.")

		print("Writing image...")
		time = get_timer()
		#Write image to file.
		imageio.imwrite(path+pathdelim+"mandelbrot_"+str(iters)+"_iterations"+colorname+eval_type+blurname+image_file_ext,result)
		time = get_timer() - time
		print("Done in "+str(time)[:4]+" seconds.")
		
		print("Total time consumption: "+str(get_timer() - total_time)[:4]+" seconds.")
