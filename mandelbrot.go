package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
	"runtime"
	"sync"
)

func main() {
	//Establish and parse command line arguments.
	verboseFlag := flag.Bool("v", false, "use this flag to print extra information")
	aspectRatioFlag := flag.Float64("r", 1.5, "set the aspect ratio of the image")
	realFractalCenterFlag := flag.Float64("cr", -0.75, "set the real part of the center of the fractal image")
	imagFractalCenterFlag := flag.Float64("ci", 0, "set the imaginary part of the center of the fractal image")
	resolutionFlag := flag.Int("y", 2016, "set the number of pixels along the y-axis")
	ssaaFlag := flag.Int("s", 3, "set the number of supersamples along one direction. Execution speed slows down with the square of this number")
	flag.Parse()

	//We use this flag to determine whether to print extra information.
	verbose := *verboseFlag

	//Work out size of plot in complex plane.
	aspectRatio := *aspectRatioFlag
	imagDistance := 8.0 / 3.0 //The distance covered along the imaginary axis in the image.
	realDistance := imagDistance * aspectRatio

	//Retrieve image center.
	realFractalCenter := *realFractalCenterFlag
	imagFractalCenter := *imagFractalCenterFlag

	//Set start and end points.
	realStart := realFractalCenter - realDistance/2.0
	realEnd := realFractalCenter + realDistance/2.0
	imagStart := imagFractalCenter - imagDistance/2.0
	imagEnd := imagFractalCenter + imagDistance/2.0

	//Colour and iteration settings.
	maxIterations := 255
	colorDepth := maxIterations

	//Data structure initialization.
	imagPointsLen := *resolutionFlag
	realPointsLen := int(float64(imagPointsLen) * aspectRatio)
	realPoints := make([]float64, realPointsLen)
	imagPoints := make([]float64, imagPointsLen)

	//Make an image object.
	upLeft := image.Point{0, 0}
	downRight := image.Point{realPointsLen, imagPointsLen}
	img := image.NewRGBA(image.Rectangle{upLeft, downRight})

	//Work out position of points in complex plane.
	linearSpace(realStart, realEnd, realPoints)
	realDelta := (realEnd - realStart) / float64(realPointsLen)
	linearSpace(imagStart, imagEnd, imagPoints)
	imagDelta := (imagEnd - imagStart) / float64(imagPointsLen)

	//Work out if we should mirror the fractal image from the top,
	//bottom, or not at all.
	//mirror := imagStart < 0 && imagEnd > 0
	//bottomHalfLarger := math.Abs(imagStart) >= math.Abs(imagEnd)
	//mirrorBottom := mirror && bottomHalfLarger
	//mirrorTop := mirror && !bottomHalfLarger

	//Establish a queue of jobs and fill it.
	rowJobs := make(chan mandelbrotRow, imagPointsLen)
	for i, cImag := range imagPoints {
		rowJobs <- newMandelbrotRow(img, cImag, i)

		//if cImag > 0 && mirrorBottom {
		//	break
		//}
	}
	close(rowJobs)

	//Parallel execution of cores separate threads that each compute
	//rows of the mandelbrot image taken from the jobs queue.
	if verbose {
		fmt.Println("Iterating...")
	}
	cores := runtime.NumCPU()
	var wg sync.WaitGroup
	for i := 0; i < cores; i++ {
		wg.Add(1) //Add one goroutine to the wait group
		go func() {
			mandelbrotRowWorker(rowJobs, realPoints, maxIterations, colorDepth, *ssaaFlag, realDelta, imagDelta)
			wg.Done()
		}() //make a gorountine that calls a mandelbrotworker, and notifies
		//the work group when it finishes.
	}
	wg.Wait() //Wait until they all finish.

	//if mirror {}

	//Encode image
	if verbose {
		fmt.Printf("Encoding fractal as png image...")
	}
	f, _ := os.Create("m.png")
	png.Encode(f, img)
}

func linearSpace(start float64, end float64, linSpace []float64) {
	length := len(linSpace)
	stepSize := (end - start) / float64(length-1)
	for i := 0; i < length; i++ {
		linSpace[i] = start + float64(i)*stepSize
	}
}

type mandelbrotRow struct {
	img    *image.RGBA
	cImag  float64
	yIndex int
}

func newMandelbrotRow(img *image.RGBA, cImag float64, yIndex int) mandelbrotRow {
	r := mandelbrotRow{img: img, cImag: cImag, yIndex: yIndex}
	return r
}

func mandelbrotRowWorker(rowJobs <-chan mandelbrotRow, cReal []float64, maxIterations int, colorDepth int, ssaa int, realDelta float64, imagDelta float64) {
	var escapeSpeed float64
	var columnOffset float64
	var rowOffset float64
	var total float64

	inverseFactor := 1.0 / float64(ssaa)

	for job := range rowJobs {
		if job.img.Bounds().Dx() != len(cReal) {
			fmt.Println("mandelbrotRowWorker: length of image row must be the same as data row. They are ", job.img.Bounds().Dx(), "and", len(cReal))
			os.Exit(1)
		}

		for j := 0; j < job.img.Bounds().Dx(); j++ {

			total = 0
			for k := 1; k <= int(math.Pow(float64(ssaa), 2.0)); k++ {
				//Computes offsets. These should range from -1/ssaa
				//to 1/ssaa with a 0 included if ssaa is odd.
				columnOffset = (float64(k%ssaa) - 1.0) * inverseFactor
				rowOffset = (float64(k-1)/float64(ssaa) - 1.0) * inverseFactor
				escapeSpeed = mandelbrotIteratePoint(cReal[j]+rowOffset*realDelta, job.cImag+columnOffset*imagDelta, maxIterations, colorDepth)
				total = total + escapeSpeed //math.Pow(escapeSpeed, 2.0) // makes the image brighter
			}

			total = total / float64(math.Pow(float64(ssaa), 2.0))

			job.img.Set(j, job.yIndex, color.RGBA{
				uint8(total * math.Pow(float64(colorDepth), 1.0-math.Pow(total, 45.0)*2.0)),
				uint8(total*70.0 - 880.0*math.Pow(total, 18.0) + 701.0*math.Pow(total, 9.0)),
				uint8(total*80.0 + math.Pow(total, 9.0)*float64(colorDepth) - 950.0*math.Pow(total, 99.0)),
				0xff})
		}
	}
}

func mandelbrotIteratePoint(cReal float64, cImag float64, maxIterations int, colorDepth int) float64 {

	cImagSqr := cImag * cImag
	magSqr := cReal*cReal + cImagSqr

	//Determine if the complex number is in the main cardioid or period two bulb.
	if math.Pow(cReal+1, 2)+cImagSqr <= 0.0625 || magSqr*(8.0*magSqr-3.0) <= 0.09375-cReal {
		//return maxIterations
		return 0.0
	}

	//Initialize variables
	zReal := 0.0
	zImag := 0.0
	iterations := 0
	zRealSqr := 0.0
	zImagSqr := 0.0

	//Iterates the mandelbrot function.
	//This loop has only three multiplications, which is the minimum.
	for iterations < maxIterations && zRealSqr+zImagSqr <= 36 {
		zImag = zReal * zImag
		zImag = zImag + zImag
		zImag = zImag + cImag
		zReal = zRealSqr - zImagSqr + cReal
		zRealSqr = zReal * zReal
		zImagSqr = zImag * zImag
		iterations++
	}

	if iterations == maxIterations {
		return 0.0
	}

	return (2.0 + float64(maxIterations-iterations-2) - 4.0*math.Pow(math.Sqrt(zRealSqr+zImagSqr), -0.4)) / float64(colorDepth)
}
