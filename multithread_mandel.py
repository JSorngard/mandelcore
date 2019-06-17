import numpy as np
import multiprocessing as mp
import sys
import imageio
import os
from mandelfortran import mandel_calc
from mandelfortran import colorize as fcolor
import time


#The highest allowed number of iterations.
iters = 100

#Defines the "window" to look at the fractal in.
start = -2.5 - 1.6j
end = .8 + 1.6j

#Number of points per axis to compute.
re_eval_points = 5000
im_eval_points = re_eval_points

#Compute it multithreaded.
multicore = True

#Make the image in color.
colorize = True

#What file type to save the image as.
file_ext = ".png"

#Color depth.
depth = 255

def mandel_helper(cs,maxiterations=iters):
	"""Takes in an array of values and calls the Fortran function for each of them."""
	result = np.zeros(np.shape(cs))
	for index,c in enumerate(cs):
		result[index] = mandel_calc(np.real(c),np.imag(c),iters)
	return result

"""def mandel_func(c,maxiterations=iters):
	#Takes a complex number and iterates the mandelbrot function on it until either 
	#its magnitude is larger than two, or it has iterated enough.
	x = np.real(c)
	y = np.imag(c)
	y2 = y**2.
	q=(x-.25)**2. + y2

	#Filter out all points inside the main bulb and the period-2 bulb.
	if((x + 1.)**2 + y2 < 0.0625 or q + (x - .25) < .25*y2):
		return maxiterations

	z = 0+0j
	iterations=0
	while(np.real(z)**2. + np.imag(z)**2. <= 4. and iterations < maxiterations):
		iterations += 1
		z = z**2. + c

	return iterations"""

#def colorize_iters(itermatrix,x,y,maxiters=iters,i=1):
#	"""Replaces the number of iteration with an RGB triplet."""
#	colormatrix = np.zeros((x,y,3))
#	for row_id,row in enumerate(itermatrix):
#		for col_id,T in enumerate(row):
#			
#			#Maps 0 to black and other numbers between 0 and 1 to a range from brown to blue.
#			colormatrix[row_id,col_id] = [T*80 + T**9*i - 950*T**99, T*70 - 880*T**18 + 701*T**9, T*i**(1 - T**45*2)]
#
#	return colormatrix

if(__name__ == "__main__"):

	#Multiplatform clock
	get_timer = time.clock if sys.platform == "win32" else time.time


	if(colorize):
		colorname = "_color"
	else:
		colorname = "_bw"

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

	print("Performing image manipulations...")
	time = get_timer()
	#Normalize result to 0-1
	result = np.array(result)/float(iters)

	if(colorize):
		colorized = np.zeros((im_eval_points,re_eval_points,3))
		colorized = fcolor(colorized,result)
		result = colorized

	else:
		#Scale up to 0-depth. What should be black is now depth.
		result *= depth 
		#Invert so that black is 0 and white is depth.
		result -= depth 
		result = np.abs(result) 
	
	#Convert to uints for imageio.
	result = result.astype(np.uint8)
	
	#Adds a flipped copy of the image to the top.
	result = np.concatenate((np.flip(result,axis=0),result))

	time = get_timer() - time
	print("Done in "+str(time)[:4]+" seconds.")

	print("Writing image...")
	time = get_timer()
	path = os.path.dirname(os.path.abspath(__file__))
	pathdelim = "\\" if sys.platform == "win32" else "/"
	imageio.imwrite(path+pathdelim+"mandel"+colorname+eval_type+file_ext,result)
	time = get_timer() - time
	print("Done in "+str(time)[:4]+" seconds.")
