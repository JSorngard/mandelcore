import numpy as np
import multiprocessing as mp
import sys
import imageio
import os

#The highest allowed number of iterations.
iters = 100

def mandel_helper(cs,maxiterations=iters):
	"""Takes in an array of values and calls mandel_func for each of them."""
	result = np.zeros(np.shape(cs))
	for index,c in enumerate(cs):
		result[index] = mandel_func(c,iters)
	return result

def mandel_func(c,maxiterations=iters):
	"""Takes a complex number and iterates the mandelbrot function on it until either 
	its magnitude is larger than two, or it has iterated enough."""
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

	return iterations

if(__name__=="__main__"):

	#Multiplatform clock
	#get_timer = time.clock if sys.platform == "win32" else time.time

	#Defines the "window" to look at the fractal in.
	start = -2.5 - 1.6j
	end = .8 + 1.6j

	#Number of points per axis to compute.
	re_points = 7000
	im_points = re_points

	multicore = True

	if(np.mod(im_points,2)==0):
		im_points = int(im_points/2)
	else:
		print("Number of imaginatry points must be even.")
		exit()

	re_points= np.linspace(np.real(start),np.real(end),re_points)
	im_points= np.linspace(0,np.imag(end),im_points)

	re_grid,im_grid = np.meshgrid(re_points,im_points*1j,sparse=True)

	grid = re_grid + im_grid
		
	gridshape = np.shape(grid)
	elements = gridshape[0]*gridshape[1]
	cmplxsize = sys.getsizeof(1+1j)
	cmplxnparraysize = sys.getsizeof(np.array(1+1j))
	cmplxnparray10size = sys.getsizeof((1+1j)*np.ones(10,dtype=complex))

	print("Size of a complex number: "+str(cmplxsize)+" B.")
	print("Size of a numpy array with a complex number: "+str(cmplxnparraysize)+" B.")
	print("Size of a numpy array with 10 complex numbers: "+str(cmplxnparray10size)+" B.")
	print("Elements in grid: "+str(elements)+".")
	print("Grid 'should' take up roughly "+str(elements*cmplxsize/1000000)+" MB.")
	print("Size of grid: "+str(sys.getsizeof(grid)/1e6)+" MB.")

	if(multicore):
		print("Attempting to evaluate on many cores...")
		filename_end = "multicore.png"

		#Create a pool with the number of threads equal to the number of processor cores.
		print("Creating thread pool...")
		cores = mp.cpu_count()
		pool = mp.Pool(processes=cores)
		#Warm up the pool
		pool.map(mandel_func,range(cores))
		print("Done.")

		print("Computing...")
		result = pool.map(mandel_helper,grid)
		print("Done.")

	else:
		print("Evaluating on a single core...")
		filename_end = "singlecore.png"
		mandel_vector = np.vectorize(mandel_func)

		print("Computing...")
		result = mandel_vector(grid,iters)
		print("Done.")

	print("Performing image manipulations...")

	#Adds a flipped copy of the image to the bottom.
	result = np.concatenate((np.flip(result,axis=0),result))

	#Normalize to 0-1.
	result = result/float(iters)
	#Scale up to 0-255. What should be black is now 255.
	result *= 255 
	#Invert so that black is 0 and white is 255.
	result -= 255 
	result = np.abs(result) 
	#Convert to uints for iamgeio.
	result = result.astype(np.uint8)


	#def colorize(T):
	#	return (T*80+T**9*255-950*T**99,T*70-880*T**18+701*T**9,T*255**(1-T**45*2))


	print("Done.")
	print("Writing image...")

	path = os.path.dirname(os.path.abspath(__file__))
	imageio.imwrite(path+"\\mandel_"+filename_end,result)

	print("Done.")
