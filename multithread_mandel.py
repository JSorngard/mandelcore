import numpy as np
import multiprocessing as mp
import sys
import os
import mandelfortran
import imagefortran
import time
import argparse
#This is probably poor form. But now it runs on a computer where I don't have the priviliges to use pip!
try:
	import imageio
	has_imageio = True
except ImportError:
	from PIL import Image
	has_imageio = False


#DEFAULT SETTINGS:

#Color depth.
depth = 255 #Integer. Set to 255 if using colorize.

#The highest allowed number of iterations.
#Set to color depth if using the scaled algorithm and colorize.
iters = depth #Integer

#Sets the aspect ratio of the image.
aspect_ratio = 3/2

#Defines parameters needed to determine the "window" to look at the fractal in.
im_dist = 8/3 #The distance along the imaginary axis where the fractal is.
fractal_center = -3/4+0j #The point to center the view on.

#Number of points per axis to compute.
im_eval_points = 1440 #y-axis. Must be an even integer.

#Compute it multithreaded.
multicore = True

#The method of parallelization.
#If true it will be openmp in Fortran, otherwise multiprocessing in python.
#Fortran is ~100 times faster.
fortran_omp = True

#Save the result as an image.
saveimage = True

#Make a pass of gaussian blur.
#I managed to get SSAA working, so there's very little reason to do this.
blur = False #Blurring requires allocating two extra images in memory,
#just as large as the main one.
#The image can therefore not be as big with this option turned on.
radius = 1 #Integer. The radius of the gaussian blur.

#Compute each pixel multiple times in different locations to smooth jaggies.
ssaa = True
#The number of sampled points along each dimension for each pixel.
ssfactor = 5 #Integer. The computation will run slower by a factor of this number squared.

#Make the image in color. Only relevant if saveimage is True.
colorize = True
#This option requires the allocation of an extra image in memory,
#three times as large as the main one.
#The image must therefore be much smaller with this option turned on.

#Zoom in on the fractal in a gif.
zoom = 1 #The factor by which to zoom.
frames = 1 #The number of frames before full zoom is achieved.
duration = 1 #The number of seconds the animation should last for.

#Raises the result of the mandelbrot iterations to this number.
gamma = .75 #Closer to 0 means a darker image.
#I take no artistic responsibillity for numbers larger than 1.
#Due to how coloring works too low gamma will result in weirdness.
#0.5 works, 0.4 does not. 0.75 looks nice.
#This isn't really a normal gamma value, it just pretends that points
#escaped later or earlier than they did and colours the image after that.

#Save the resulting iteration grid to file
saveresult = False

#What file type to save the image as.
image_file_ext = ".bmp"
#Trials at im_eval_points = 10000, aspect_ratio = 3/2:
#bmp: 47.5 seconds to encode a 450 MB image. No compression.
#jpg: 54 seconds to encode a 3.8 MB image. Lossy compression.
#png: 182 seconds to encode an 18.2 MB image. Lossless compression.
#Only relative numbers are probably useful:
#bmp: 1 to encode a size 1 image.
#jpg: 1.14 to encode a size 0.008 lossy image.
#png: 3.83 to encode a size 0.040 lossless image.

#What file type to save the raw data as. If it ends in .gz it will be compressed.
data_file_ext = ".dat.gz"

#Whether to print more detailed information during image generation.
#Will be set to true if frames == 1.
debug = False

#Print extra information about memory use.
memory_debug = True
#Set to False to print one less line. Hooray!

#Define the different possible command line paramters.
parser = argparse.ArgumentParser(description="Computes and saves an image of the mandelbrot set.")
parser.add_argument("-c","--center",required=False,type=complex,default=fractal_center,help="Specify the point in the complex plane to center the image on. Defaults to "+str(fractal_center)+".")
parser.add_argument("-y","--yresolution",required=False,type=int,default=im_eval_points,help="Specify the y-axis resolution of the image. Defaults to "+str(im_eval_points)+".")
parser.add_argument("-r","--aspectratio",required=False,type=float,default=aspect_ratio,help="Specifies the aspect ratio of the image. Defaults to "+str(aspect_ratio)+".")
parser.add_argument("-z","--zoom",required=False,type=float,default=1.,help="The zoom in factor to reach.")
parser.add_argument("-f","--frames",required=False,type=int,default=1,help="The number of frames per until the full zoom is achieved.")
parser.add_argument("-d","--duration",required=False,type=float,default=duration,help="The number of seconds the animation should take (applicable if frames and zoom > 1).")
parser.add_argument("-g","--gamma",required=False,type=float,default=gamma,help="Raises the output of the mandelbrot iterations to this number. Works as a gamma between 0.4 and 1.")
parser.add_argument("-s","--ssaafactor",required=False,type=int,default=ssfactor,help="Supersample each pixel this many times squared. Defaults to "+str(ssfactor)+".")
parser.add_argument("-e","--fileextension",required=False,default=image_file_ext,help="Set the file extension of the generated image. Defaults to "+image_file_ext[1:]+".")
parser.add_argument("--saveresult",required=False,action="store_true",help="Use this argument if you want to save the results of the mandelbrot iterations to a "+data_file_ext+" file.")
parser.add_argument("--noimage",required=False,action="store_false",help="Use this argument if you do not want the program to output an image file.")
parser.add_argument("--debug",required=False,action="store_true",help="Use this argument if you want more detailed information on the computation. Will be on if only generating one image.")

#Extract the given arguments.
args=vars(parser.parse_args())

#Put them in appropriate variables.
im_eval_points = args["yresolution"]
re_eval_points = int(aspect_ratio*im_eval_points)
gamma = args["gamma"]
ssfactor = args["ssaafactor"]
zoom = args["zoom"]
frames = args["frames"]
duration = args["duration"]
image_file_ext = args["fileextension"]
saveresult = args["saveresult"]
saveimage = args["noimage"]
fractal_center = args["center"]
aspect_ratio = args["aspectratio"]
debug = args["debug"]

#start = -2.7-1.333j #Good for 3/2 aspect ratio.
#end = 1.3+1.333j

#start = -3.5 - 4/3j #Good for 16/9.
#end = 1.6 + 4/3j

#Multiplatform clock
get_time = time.perf_counter if sys.platform == "win32" else time.time

#Determines the working directory of the program.
path = os.path.dirname(os.path.abspath(__file__))

def get_window(fractal_center,im_dist,aspect_ratio,zoom):
	im_dist *= 1/zoom
	re_dist = im_dist*aspect_ratio
	start = fractal_center - re_dist/2 - im_dist/2*1j
	end = fractal_center + re_dist/2 + im_dist/2*1j
	return start,end

if(not fortran_omp):
	def mandel_func(c,maxiterations=iters,colordepth=float(depth)):
		"""Takes a complex number and iterates the mandelbrot function on it until either its magnitude is larger than six (2 gives worse colour fade), or it has iterated enough."""
		x = np.real(c)
		y = np.imag(c)
		y2 = y**2.
		#q=(x-.25)**2. + y2
		mag2=x**2. + y2

		#Filter out all points inside the main bulb and the period-2 bulb.
		if((x + 1.)**2 + y2 < 0.0625 or mag2*(8.*mag2-3.) <= .09375 - x):
			return 0.

		z = 0.+0.j
		iterations=0
		while(np.real(z)**2. + np.imag(z)**2. <= 36. and iterations < maxiterations):
			iterations += 1
			z = z**2. + c
		return (2+(maxiterations-iterations-2.)-4.*np.abs(z)**(-.4))/colordepth if iterations != maxiterations else 0.

	def mandel_helper(cs,maxiterations=iters):
		return [mandel_func(c,iters) for c in cs]

def mandelbrot(start,end,re_eval_points,im_eval_points,aspect_ratio,depth=depth,iters=iters,multicore=multicore,saveimage=saveimage,blur=blur,radius=radius,ssaa=ssaa,ssfactor=ssfactor,colorize=colorize,gamma=gamma,path=path,image_file_ext=image_file_ext,data_file_ext=data_file_ext,memory_debug=memory_debug,debug=False):

		#Only compute half the points along the imaginary axis since there is a reflection symmetry.
		if(np.mod(im_eval_points,2)==0):
			im_eval_points = int(im_eval_points/2)
		else:
			print("Number of imaginary points must be even.")
			return 0
		
		if(debug):
			print("Generating "+str(re_eval_points)+" by "+str(im_eval_points)+" grid...")
			time = get_time()

		re_points= np.linspace(np.real(start),np.real(end),re_eval_points)
		deltar = re_points[1] - re_points[0]
		im_points= np.linspace(0,np.imag(end),im_eval_points)
		deltai = im_points[1] - im_points[0]

		re_grid,im_grid = np.meshgrid(re_points,im_points*1j,sparse=True)

		try:
			grid = re_grid + im_grid
		except MemoryError:
			print("Out of memory when allocating grid.")
			grid = None
			re_points = None
			im_points = None
			return 0

		re_points = None
		im_points = None

		
		if(debug):
			time = get_time() - time
			print("Done in "+str(time)[:4]+" seconds.")
		
		if(memory_debug and debug):		
			gridshape = np.shape(grid)
			elements = gridshape[0]*gridshape[1]
			cmplxsize = sys.getsizeof(1+1j)
			#cmplxnparraysize = sys.getsizeof(np.array(1+1j))
			#cmplxnparray10size = sys.getsizeof((1+1j)*np.ones(1,dtype=complex))

			#print("Size of a complex number: "+str(cmplxsize)+" B.")
			#print("Size of a numpy array with a complex number: "+str(cmplxnparraysize)+" B.")
			#print("Size of a numpy array with 10 complex numbers: "+str(cmplxnparray10size)+" B.")
			#print("Elements in grid: "+str(elements)+".")
			print("Grid should take up roughly "+str(elements*cmplxsize/1e6)+" MB in RAM.")
			#print("Size of grid: "+str(sys.getsizeof(grid)/1e6)+" MB.")

		if(multicore):
			cores = mp.cpu_count()
			if(debug):
				print("Attempting to evaluate on "+str(cores)+" cores...")
				time = get_time()

			if(fortran_omp):
				try:
					if(ssaa):						
						if(debug):
							print("Computing with SSAAx"+str(ssfactor**2)+"...")
						grid = mandelfortran.mandel_calc_array_scaled_supersampled(grid,iters,depth,ssfactor,deltar,deltai)
					else:
						if(debug):
							print("Computing...")
						grid = mandelfortran.mandel_calc_array_scaled(grid,iters,depth)
				except MemoryError:
					print("Out of memory when sending work to Fortran.")
					grid = None
					return 0
				result = grid
				
				if(debug):
					time = get_time() - time
					print("Done in "+str(time)[:4]+" seconds.")			

			else:
				#Create a pool with the number of threads equal to the number of processor cores.
				if(debug):
					print("Creating thread pool...")
				pool = mp.Pool(processes=cores)
				#Warm up the pool
				pool.map(mandel_helper,np.ones((10,cores)))
				time = get_time() - time
				if(debug):
					print("Done in "+str(time)[:5]+" seconds.")
					print("Computing...")
					time = get_time()

				result = pool.map(mandel_helper,grid)
				
				if(debug):
					time = get_time() - time
					print("Done in "+str(time)[:4]+" seconds.")

		else:
			if(debug):
				print("Evaluating on a single core...")

			eval_type = "_singlecore"

			mandel_vector = np.vectorize(mandel_func)

			if(debug):
				print("Computing...")
				time = get_time()
			
			result = mandel_vector(grid,iters)
			
			if(debug):
				time = get_time() - time
				print("Done in "+str(time)[:4]+" seconds.")
		

		grid = None #Removes the grid of complex values from memory.
		
		#Extract the real part of the output of the fortran subroutine.
		#No information is lost since it only writes to the real part.
		result = np.real(result)

		if(saveimage):
			if(debug):
				print("Performing image manipulations...")
				time = get_time()
			
			gammaname = ""
			if(gamma != 1.):
				if(debug):
					print(" changing gamma...")
				#If the image is large enough we compute it multithreaded.
				#Don't do this. WAY slower than np.power. Gotta love vectorization.
				#if(im_eval_points > 5000):
				#	result = mandelfortran.multicore_pow(result,gamma)
				#else:
				try:
			 		result = np.power(result,gamma)
				except MemoryError:
					print("Out of memory when changing gamma.")
					result = None
					return 0

			if(blur):
				if(debug):
					print(" blurring...")
				try:
					blurred = np.zeros(np.shape(result))
					blurred = imagefortran.fastgauss(result,radius)
					result = blurred
				except MemoryError:
					print("Out of memory when blurring the image.")
					result = None
					blurred = None
					return 0

				blurred = None

			if(colorize):
				if(debug):
					print(" colouring...")
				try:
					colourized = np.zeros((np.concatenate((np.shape(result),np.array([3])))),order='F')
					colourized = imagefortran.fcolour(depth,colourized,np.real(result))
					result = colourized
				except MemoryError:
					print("Out of memory when colouring the image.")
					colourized = None
					result = None
					return 0

				colourized = None
			else:
				if(debug):
					print(" fitting to color depth...")
				#Scale up to 0-depth.
				result *= depth

				#Invert so that black is 0 and white is depth.
				#result -= depth 
				#result = np.abs(result) 
		
			#Convert to uints for image saving.
			result = result.astype(np.uint8)

			if(debug):
				print(" mirroring...")
			#Adds a flipped copy of the image to the top.
			try:
				result = np.concatenate((np.flip(result,axis=0),result))
			except AttributeError:
				try:
					result = np.concatenate((result[::-1],result))
				except MemoryError:
					print("Out of memory when mirroring image.")
					result = None
					return 0
			except MemoryError:
				print("Out of memory when mirroring image.")
				result = None
				return 0

			if(debug):
				time = get_time() - time
				print("Done in "+str(time)[:4]+" seconds.")

		return result
	
def write_image(fullname,image_file_ext,result,duration=1,debug=False):
	print("Writing image...")
	time = get_time()

	#Find the number of separate images in the result array.
	frames = np.shape(result)[0]

	#Write image to file.
	if(has_imageio):
		if(debug):
			print(" using imageio...")
		if(frames == 1):
			imageio.imwrite(fullname+image_file_ext,result[0])
		elif(frames > 1):
			image_file_ext = ".gif"
			imageio.mimwrite(fullname+image_file_ext,result,duration=duration/frames)
	elif(frames == 1):
		if(debug):
			print(" using PIL...")
			print("  converting to image object...")
		try:
			result = Image.fromarray(result[0])
		except OverflowError:
			print("  The image array is too large for PIL to handle. Try installing imageio.")
			result = None
			return 0

		if(debug):
			print("  saving...")
		result.save(fullname,optimize=True,quality=85)

	else:
		print("Can currrently not save gif with PIL.")
		return 0

	time = get_time() - time
	print("Done in "+str(time)[:4]+" seconds.")
	return 1	

def write_data(fullname,data_file_ext,result,debug=False):
	print("Writing raw data...")
	if(data_file_ext[-3:] == ".gz" and debug):
		print(" compressing...")
	if(debug):
		time = get_time()

	#Check if there is data from multiple frames to write.
	if(len(np.shape(result)) == 4):
		files = np.shape(result)[0]
		print("  WARNING: writing data from multiple frames.")
	else:
		files = 1

	#Write iteration data to file. If the file name ends in .gz numpy automatically compresses it.
	#Maybe in the future I'll be able to use this data to make an image.
		
	for i in range(files):
		index_name = "_"+str(i+1) if files > 1 else ""
		np.savetxt(fullname+index_name+data_file_ext,result,delimiter=' ')
	
	if(debug):
		time = get_time() - time
		print("Done in "+str(time)[:4]+" seconds.")
	return 1


if(__name__ == "__main__"):
	
	total_time = get_time()
	
	print("START")

	#Determines whether to use \ or / for file paths.
	pathdelim = "\\" if sys.platform == "win32" else "/"	
	
	#Generate the file name to save the image as.
	colorname = "_color" if colorize else "_bw"

	blurname = "_blur="+str(radius) if blur else ""

	eval_type = ""
	#eval_type = "_multicore" if multicore else "_singlecore" #Can be useful for debugging.

	gammaname = "_g="+str(gamma)[:4]

	ssaaname = "_ssaax"+str(ssfactor**2) if ssaa else ""

	#Avoid unneccesary work.
	if(frames > 1 and zoom == 1):
		frames = 1
		print("More than one frame has been requested, but at no zoom. Generating one image instead.")
	
	#The other alternatives (frames == 1 and zoom != 1) and (frames == 1 and zoom == 1)
	#are both allowed.

	#Always print details when making a single image.
	if(frames == 1):
		debug = True

	if(not saveresult and not saveimage):
		print("Note: program will produce no output.")

	result = []
	for i in range(frames):

		if(frames > 1):
			print("---Generating frame "+str(i+1)+"/"+str(frames)+"---")

		#Compute the region of the complex plane to plot.
		start,end = get_window(fractal_center,im_dist,aspect_ratio,1.+(zoom-1)/frames*i)
		
		#Generate an RGB matrix of the fractal.
		frame = mandelbrot(start,end,re_eval_points,im_eval_points,aspect_ratio,debug=debug)

		#If the result of the computation is 0 there was an error.
		if(type(frame) == int and result == 0):
			frame = None
			exit()

		#Save the generated image.
		result.append(frame)

		print("")

	#Clear up some memory.
	frame = None

	#Save the data as a file and not an image.
	if(saveresult):
		filename = path+pathdelim+"mandel"+colorname+eval_type
		success = write_data(filename,data_file_ext,result)	
		
		if(not success):
			exit()

	#Save an image.
	if(saveimage):
		filename = path+pathdelim+"mandelbrot_"+str(iters)+"_iters"+colorname+ssaaname+eval_type+blurname+gammaname
		success = write_image(filename,image_file_ext,result,duration)

		if(not success):
			exit()
		
	#Clear memory.
	result = None

	print("Total time consumption: "+str(get_time() - total_time)[:5]+" seconds.")
	print("END")