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
try:
	from numba import njit
	has_numba = True
except ImportError:
	has_numba = False

#DEFAULT SETTINGS:

#Color depth.
depth = 255 #Integer. Set to 255 if using colorize.

#The highest allowed number of iterations.
#Set to color depth if using the scaled algorithm and colorize.
iters = depth #Integer

#Sets the aspect ratio of the image.
#aspect_ratio = 21./9.
#aspect_ratio = 16./9.
aspect_ratio = 3./2.

#Defines parameters needed to determine the "window" to look at the fractal in.
im_dist = 8./3. #The distance along the imaginary axis where the fractal is.
fractal_center = -3./4.+0j #The point to center the view on.

#Number of points per axis to compute.
im_eval_points = 2160 #y-axis. Must be an even integer.
#im_eval_points = 1440
#im_eval_points = 1080

#Compute it multithreaded.
multicore = True

#The method of parallelization.
#If true it will be openmp in Fortran, otherwise multiprocessing in python.
#Fortran is ~100 to ~1000 times faster.
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
ssfactor = 3 #Integer. The computation will run slower by a factor of this number squared.

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
gamma = 1 #Closer to 0 means a darker image.
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

#Set to true if you want the fastest escaping point to always be blue.
colour_shift = False

#The default compression level for a png image. 0 if fast but large, 9 is slow but small.
compress_level = 6

#Define the different possible command line paramters.
parser = argparse.ArgumentParser(description="Computes and saves an image of the mandelbrot set.")
parser.add_argument("-cr","--Rcenter",required=False,type=float,default=np.real(fractal_center),help="Specify the real part of the point in the complex plane to center the image on. Defaults to "+str(np.real(fractal_center))+".")
parser.add_argument("-ci","--Icenter",required=False,type=float,default=np.imag(fractal_center),help="Specify the imaginary part of the point in the complex plane to center the image on. Defaults to "+str(np.imag(fractal_center))+".")
parser.add_argument("-y","--yresolution",required=False,type=int,default=im_eval_points,help="Specify the y-axis resolution of the image. Defaults to "+str(im_eval_points)+".")
parser.add_argument("-r","--aspectratio",required=False,type=float,default=aspect_ratio,help="Specifies the aspect ratio of the image. Defaults to "+str(aspect_ratio)+".")
parser.add_argument("-z","--zoom",required=False,type=float,default=1.,help="The zoom in factor to reach.")
parser.add_argument("-b","--blackwhite",required=False,action="store_false",help="Use this command if you want the image to be black and white.")
parser.add_argument("-C","--colourshift",required=False,action="store_true",help="Use this flag if you want to shift the colouring so that the fastest escaping point is always blue.")
parser.add_argument("-f","--frames",required=False,type=int,default=1,help="The number of frames per until the full zoom is achieved.")
parser.add_argument("-d","--duration",required=False,type=float,default=duration,help="The number of seconds the animation should take (applicable if frames and zoom > 1).")
parser.add_argument("-g","--gamma",required=False,type=float,default=gamma,help="Raises the output of the mandelbrot iterations to this number. Works as a gamma between 0.4 and 1.")
parser.add_argument("-s","--ssaafactor",required=False,type=int,default=ssfactor,help="Supersample each pixel this many times squared. Defaults to "+str(ssfactor)+". If 1, no SSAA will be computed.")
parser.add_argument("-e","--fileextension",required=False,default=image_file_ext,help="Set the file extension of the generated image. Defaults to "+image_file_ext[1:]+".")
parser.add_argument("-l","--compresslevel",required=False,type=int,default=compress_level,help="Specify the compression level of the image if the file extension is set to png. Defaults to "+str(compress_level)+".")
parser.add_argument("--saveresult",required=False,action="store_true",help="Use this argument if you want to save the results of the mandelbrot iterations to a "+data_file_ext+" file.")
parser.add_argument("--noimage",required=False,action="store_false",help="Use this argument if you do not want the program to output an image file.")
parser.add_argument("-v","--verbose",required=False,action="store_true",help="Use this argument if you want more detailed information on the computation.")

#Extract the given arguments.
args=vars(parser.parse_args())

#Put them in appropriate variables.
fractal_center = args["Rcenter"]-1j*args["Icenter"]
aspect_ratio = args["aspectratio"]
im_eval_points = args["yresolution"]
re_eval_points = int(round(aspect_ratio*im_eval_points))
gamma = args["gamma"]
ssfactor = args["ssaafactor"]
if ssfactor == 1:
	ssaa = False
zoom = args["zoom"]
colorize = args["blackwhite"]
colour_shift = args["colourshift"]
frames = args["frames"]
duration = args["duration"]
image_file_ext = args["fileextension"]
compress_level = args["compresslevel"]
saveresult = args["saveresult"]
saveimage = args["noimage"]
debug = args["verbose"]

if not has_imageio and frames > 1 and saveimage:
	print("Since the number of requested frames is larger than one, the output file type has been set to gif. But this computer does not have imageio installed, and PIL can not save gifs. Install imageio, or request only one image, and try again.")
	exit()

#start = -2.7-1.333j #Good for 3/2 aspect ratio.
#end = 1.3+1.333j

#start = -3.5 - 4/3j #Good for 16/9.
#end = 1.6 + 4/3j

#Multiplatform clock
get_time = time.perf_counter if sys.platform == "win32" else time.time

#Determines the working directory of the program.
path = os.path.dirname(os.path.abspath(__file__))
"""
if not fortran_omp:
	def mandel_func(c,maxiterations=iters,colordepth=float(depth)):
		#Takes a complex number and iterates the mandelbrot function on it until either its magnitude is larger than six (2 gives worse colour fade), or it has iterated enough.
		x = np.real(c)
		y = np.imag(c)
		y2 = y**2.
		#q=(x-.25)**2. + y2
		mag2=x**2. + y2

		#Filter out all points inside the main bulb and the period-2 bulb.
		if (x + 1.)**2 + y2 < 0.0625 or mag2*(8.*mag2-3.) <= .09375 - x:
			return 0.

		z = 0.+0.j
		iterations=0
		while(np.real(z)**2. + np.imag(z)**2. <= 36. and iterations < maxiterations):
			iterations += 1
			z = z**2. + c
		return (2+(maxiterations-iterations-2.)-4.*np.abs(z)**(-.4))/colordepth if iterations != maxiterations else 0.

	def mandel_helper(cs,maxiterations=iters):
		return [mandel_func(c,iters) for c in cs]
"""
#Makes exponentiation multicore if applicable on current machine
#and beneficial for the problem size.
cpu_cores = mp.cpu_count()
if has_numba and cpu_cores > 1 and im_eval_points >= 4000:
	@njit(parallel=True)
	def powah(array,power):
		return np.power(array,power)
else:
	def powah(array,power):
		return np.power(array,power)

#Takes in a quantity, shortens it and adds on the appropriate SI suffix.
def quantity_suffix(size):

	datasuffixes = {
		0: "",
		3: "k",
		6: "M",
		9: "G",
		12: "T" #Please don't ever need this.
	}

	exponent = 3*(int(np.log10(size))//3)

	return str(round(size/10**(exponent),1))[:5]+" "+datasuffixes.get(exponent)

def mandelbrot(fractal_center,im_dist,re_eval_points,im_eval_points,aspect_ratio,zoom,depth=depth,iters=iters,multicore=multicore,saveimage=saveimage,blur=blur,radius=radius,ssaa=ssaa,ssfactor=ssfactor,colorize=colorize,colour_shift=False,gamma=gamma,path=path,image_file_ext=image_file_ext,data_file_ext=data_file_ext,memory_debug=memory_debug,debug=False):

		im_dist *= 1/zoom
		re_dist = im_dist*aspect_ratio
		start = fractal_center - re_dist/2 - im_dist/2*1j
		end = fractal_center + re_dist/2 + im_dist/2*1j

		#We can get away with doing only half the work if we are centered on the real axis due to mirror symmetry.
		mirror = np.imag(fractal_center) == 0
		if mirror:
			if np.mod(im_eval_points,2)==0:
				im_eval_points = int(im_eval_points/2)
			else:
				print("Number of imaginary points must be even.")
				return 0

#----------------------------Grid generation section-------------------------------------
		
		if debug:
			print("Generating "+str(re_eval_points)+" by "+str(im_eval_points)+" grid...")
			time = get_time()
			
			grid_size = np.ones((re_eval_points,im_eval_points,3)).nbytes
			if mirror:
				grid_size *= 2


		#Generates two 1d arrays, one for the real parts and one for the imaginary parts
		#of all coordinates in the complex plane we want to compute the fractal for.
		
		#Real parts
		re_points= np.linspace(np.real(start),np.real(end),re_eval_points)
		deltar = re_points[1] - re_points[0]
		
		#Imaginary parts
		if mirror:
			im_points= np.linspace(0,np.imag(end),im_eval_points)
		else:
			im_points= np.linspace(np.imag(start),np.imag(end),im_eval_points)
		deltai = im_points[1] - im_points[0]

		
		if debug:
			time = get_time() - time
			
			grid_size += re_points.nbytes + im_points.nbytes

			if memory_debug:
				print(" fractal grid should take up roughly "+quantity_suffix(grid_size)+"B in memory.")
			
			print("Done in "+str(time)[:4]+" seconds.\n")
			

#----------------------------Fractal evaluation section----------------------------------

		cores = mp.cpu_count()
		if debug:
			if im_eval_points < 300: #The 300 limit is hard coded into the fortran code as of right now.
				print("Evaluating with a single core due to the image size...")
			else:
				print("Attempting to evaluate on "+str(cores)+" cores...")
			time = get_time()

		try:
			if ssaa :
				if debug:
					print("Computing with SSAAx"+str(ssfactor**2)+"...")
				result = mandelfortran.iterate_supersampled(re_points,im_points,iters,depth,ssfactor,deltar,deltai)
			else:
				if debug:
					print("Computing...")
				result = mandelfortran.iterate(re_points,im_points,iters,depth)
		except MemoryError:
			print("Out of memory when sending work to Fortran.")
			return 0
		
		if debug:
			time = get_time() - time
			print("Done in "+str(time)[:4]+" seconds.\n")			

#----------------------------Image post-processing section-------------------------------

		if saveimage:
			if debug:
				print("Performing image manipulations...")
				time = get_time()
			
			gammaname = ""
			if gamma != 1.:
				if debug:
					print(" changing gamma...")
				try:
			 		result = powah(result,gamma)
			 		#result = mandelfortran.multicore_pow(result,gamma)
				except MemoryError:
					print("Out of memory when changing gamma.")
					result = None
					return 0

			if blur:
				if debug:
					print(" blurring...")
				try:
					if debug:
						print("  generating target...")
					blurred = np.zeros(np.shape(result))
					if debug:
						print("  computing gaussian blur...")
					blurred = imagefortran.fastgauss(result,radius)
					result = blurred
				except MemoryError:
					print("Out of memory when blurring the image.")
					result = None
					blurred = None
					return 0

				blurred = None

			if colorize:
				if debug:
					print(" colouring...")
				try:
					#Shifts the colouring so that the fastest escaping point is blue.
					if colour_shift:
						if debug:
							print("  scaling...")
						result = np.multiply(result,.98/np.max(result))
					
					if debug:
						print("  generating target...")
					colourized = np.zeros((np.concatenate((np.shape(result),np.array([3])))),order='F')
					
					if debug:
						print("  computing colours...")
					colourized = imagefortran.fcolour(depth,colourized,np.real(result))
					result = colourized
				except MemoryError:
					print("Out of memory when colouring the image.")
					colourized = None
					result = None
					return 0

				colourized = None
			else:
				if colour_shift:
					if debug:
						print(" scaling")
					result = np.multiply(result,.98/np.max(result))

				if debug:
					print(" fitting to color depth...")
				#Scale up to 0-depth.
				result *= depth

				#Invert so that black is 0 and white is depth.
				#result -= depth 
				#result = np.abs(result) 
		
			#Convert to uints for image saving.
			result = result.astype(np.uint8)

			if mirror:
				if debug:
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

				if debug:
					time = get_time() - time
					#Need an extra digit of precision if the image is small.
					accuracy = 4 if im_eval_points > 200 else 5
					print("Done in "+str(time)[:accuracy]+" seconds.\n")

		return result
	
def write_image(fullname,image_file_ext,result,has_imageio,duration=1,debug=False,compress_level=compress_level):
	if debug:
		print("Writing image...")

	time = get_time()

	#Find the number of separate images in the result array.
	#frames = np.shape(result)[0] #This takes a weirldy large ammount of time.
	frames = len(result)
	
	#If we have generated multiple images we are making a gif.
	if frames > 1:
		image_file_ext = ".gif"

	#Write image to file.
	if has_imageio:
		if debug:
			print(" using imageio...")
			print("  saving "+fullname+image_file_ext+"...")
		if frames == 1:
			if image_file_ext == ".png":
				imageio.imwrite(fullname+image_file_ext,result[0],compress_level=compress_level) #Higher compress level oes not seem to do much for the file size.
			else:
				imageio.imwrite(fullname+image_file_ext,result[0])

		elif frames > 1:
			imageio.mimwrite(fullname+image_file_ext,result,duration=duration/frames)
	elif frames == 1:
		if debug:
			print(" using PIL...")
			print("  converting array to image object...")
		try:
			result = Image.fromarray(result[0])
		except OverflowError:
			print("  The image array is too large for PIL to handle. Try installing imageio.")
			result = None
			return 0

		if debug:
			print("  saving "+fullname+image_file_ext+"...")
		result.save(fullname+image_file_ext,optimize=True,quality=85)

	else:
		print("Can not use PIL to save a gif. Try instaling imageio.")
		return 0
		if debug:
			print(" using PIL...")
			print("  converting arrays to image objects...")
		try:
			result = [Image.fromarray(array) for array in result]
		except OverflowError:
			print("  The image arrays are too large for PIL to handle. Try installing imageio.")
			result = None
			return 0
		
		if debug:
			print("  saving "+fullname+image_file_ext+"...")
		result[0].save(fullname+image_file_ext,format='GIF',save_all=True,append_images=result[1:],duration=duration/frames)
		
		#print("Can currrently not save gif with PIL.")
		#return 0

	time = get_time() - time
	
	if debug:
		size = os.stat(fullname+image_file_ext).st_size
		print(" saved "+fullname+image_file_ext+" with size "+quantity_suffix(size)+"B.")
		print("Done in "+str(time)[:4]+" seconds.\n")
	
	return 1	

def write_data(fullname,data_file_ext,result,debug=False):
	if debug:
		print("Writing raw data...")

	if data_file_ext[-3:] == ".gz" and debug:
		print(" compressing...")
	if debug:
		time = get_time()

	#Check if there is data from multiple frames to write.
	if len(np.shape(result)) == 4:
		files = np.shape(result)[0]
		print("  WARNING: writing data from multiple frames.")
	else:
		files = 1

	#Write iteration data to file. If the file name ends in .gz numpy automatically compresses it.
	#Maybe in the future I'll be able to use this data to make an image.
		
	for i in range(files):
		index_name = "_"+str(i+1) if files > 1 else ""
		np.savetxt(fullname+index_name+data_file_ext,result,delimiter=' ')
	
	if debug:
		time = get_time() - time
		print("Done in "+str(time)[:4]+" seconds.\n")
	return 1


if __name__ == "__main__":
	
	total_time = get_time()
	
	#Determines whether to use \ or / for file paths.
	pathdelim = "\\" if sys.platform == "win32" else "/"	
	
	#Generate the file name to save the image as.
	colorname = "_colour" if colorize else "_bw"

	blurname = "_blur="+str(radius) if blur else ""

	eval_type = ""
	#eval_type = "_multicore" if multicore else "_singlecore" #Can be useful for debugging.

	gammaname = "_g="+str(gamma)[:4]

	ssaaname = "_ssaax"+str(ssfactor**2) if ssaa else ""

	#Avoid unneccesary work.
	if frames > 1 and zoom == 1:
		frames = 1
		print("More than one frame has been requested, but at no zoom. Generating one image instead.")
	
	#The other alternatives (frames == 1 and zoom != 1) and (frames == 1 and zoom == 1)
	#are both allowed.

	#Always print details when making a single image.
	#if frames == 1:
	#	debug = True

	if not saveresult and not saveimage:
		print("Note: program will produce no output.")

	result = []
	for i in range(frames):

		if frames > 1 and debug:
			print("---Generating frame "+str(i+1)+"/"+str(frames)+", "+str(100*float(i+1)/float(frames))[:5]+"%---")
		
		if frames == 1:
			z = zoom
		else:
			#Computes the ammount of zoom in this frame. 
			#z = 1.+(zoom-1)/frames*i #This slows down towards the end.
			#z = 1 - i/frames + zoom*i**2/frames**2 #This slows down less.
			#z = (zoom-1)/frames**2*i**2 + 1 #This is similar, but more efficient computation.
			#z = (zoom-1)/(np.exp(frames)-1)*(np.exp(i)-1) + 1 #Jump scare at the end.
			z = (zoom**(1./frames))**i #Constant relative speed.
		
		if debug and frames > 1:
			time = get_time()

		#Generate an RGB matrix of the fractal.
		frame = mandelbrot(fractal_center,im_dist,re_eval_points,im_eval_points,aspect_ratio,z,debug=debug,colour_shift=colour_shift)
		
		#If the result of the computation is 0 there was an error.
		if type(frame) == int and frame == 0:
			frame = None
			exit()

		if debug and frames > 1:
			print("---Frame finished in "+str(get_time() - time)[:4]+" seconds---\n")

		#Save the generated image.
		result.append(frame)

		

	#Clear up some memory.
	frame = None

	#Save the data as a file and not an image.
	if saveresult:
		filename = path+pathdelim+"mandel"+colorname+eval_type
		success = write_data(filename,data_file_ext,result,debug=debug)	
		
		if not success:
			exit()

	#Save an image.
	if saveimage:
		filename = path+pathdelim+"mandelbrot_"+str(iters)+"_iters"+colorname+ssaaname+eval_type+blurname+gammaname
		filename = "m"
		success = write_image(filename,image_file_ext,result,has_imageio,duration=duration,debug=debug)

		if not success:
			exit()
		
	#Clear memory.
	result = None

	if debug:
		print("Total time consumption: "+str(get_time() - total_time)[:5]+" seconds.")
