/*
    Create an LMS Noise Canceller for a noisy signal.  Choose the
    sinusoidal frequency.  Plot the impulse and frequency responses
	of the filter.  Filter the noisy signal using the canceller.

	Plot the LMS adaptive filter in complex(real, imag) from a file, one complex
	number  per line.
	r1 i1
	r2 i2
	...
	rn in

	The html/template package is used to generate the html sent to the client.
	Use CSS display grid to display a 300x300 grid of cells.
	Use CSS flexbox to display the labels on the x and y axes.

	Calculate the power spectral density (PSD) using Welch periodogram spectral
	estimation and plot it.  Average the periodograms using 50% overlap with a
	user-chosen window function such as Bartlett, Hanning, Hamming, or Welch.
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"net/http"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"
	"text/template"

	"github.com/mjibson/go-dsp/fft"
)

const (
	rows                     = 300                                // #rows in grid
	columns                  = 300                                // #columns in grid
	block                    = 512                                // size of buf1 and buf2, chunks of data to process
	tmpltime                 = "templates/plottimedata.html"      // html template address
	tmplfrequency            = "templates/plotfrequencydata.html" // html template address
	tmplplotresponse         = "templates/plotresponse.html"      // html template address
	tmplLmsNoiseCanceller    = "templates/LMSnoisecanceller.html" // html template address
	tmplfiltersignal         = "templates/filtersignal.html"      // html template address
	addr                     = "127.0.0.1:8080"                   // http server listen address
	patternFilterSignal      = "/filtersignal"                    // http handler pattern for filtering using the LMS Noise Canceller
	patternLmsNoiseCanceller = "/lmsnoisecanceller"               // http handler pattern for creating LMS noise canceller
	patternPlotResponse      = "/plotresponse"                    // http handler for plotting impulse or frequency responses
	xlabels                  = 11                                 // # labels on x axis
	ylabels                  = 11                                 // # labels on y axis
	dataDir                  = "data/"                            // directory for the data files
	deg2rad                  = math.Pi / 180.0                    // convert degrees to radians
	stateFile                = "filterstate.txt"                  // last M-1 filtered outputs from prev
	signal                   = "signal.txt"                       // signal to plot
	lmsFilterNoiseCanceller  = "lmsnoisecanceller.txt"            // LMS adaptive filter for Noise Canceller
	twoPi                    = 2.0 * math.Pi                      // Two PI ~ 6.28
)

// Type to contain all the HTML template actions
type PlotT struct {
	Grid        []string // plotting grid
	Status      string   // status of the plot
	Xlabel      []string // x-axis labels
	Ylabel      []string // y-axis labels
	Filename    string   // filename to plot
	SampleFreq  string   // data sampling rate in Hz
	FFTSegments string   // FFT segments, K
	FFTSize     string   // FFT size
	Samples     string   // complex samples in data file
	SNR         string   // signal-to-noise ratio
	SigFreq     string   // sine wave frequency in Hz
	PoleRad     string   // channel pole radius
	PoleAng     string   // channel pole angle in degrees
}

// Type to hold the minimum and maximum data values
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

// Window function type
type Window func(n int, m int) complex128

// previous sample block properties used for generating/filtering current block
type FilterState struct {
	lastFiltered    []float64 // last M-1 incomplete filtered samples from previous block
	firstSampleTime float64   // start time of current submit
	lastSampleTime  float64   // end time of currrent submit
}

// data to be supplied between the pipeline stages through the channel
type ComChan struct {
	desired float64 // primary source signal to the LMS algorithm
	in      float64 // input to the stage in the pipeline
}

// Data to construct the processing chain for display
type FilterSignal struct {
	toChan      chan float64 // synchronized channel to noisy channel from signal source generator
	fromChan    chan float64 // synchronized channel to adaptive filter from noisy channel
	wg          sync.WaitGroup
	samples     int       // total number of samples per submit
	sampleFreq  int       // sample frequency in Hz
	snr         int       // signal to noise ratio in dB
	FilterState           // used by current sample block from previous sample block
	filterCoeff []float64 // filter coefficients for the adaptive filter
	filterfile  string    // name of the adaptive filter file
	Endpoints             // embedded struct
	poleRad     float64   // noisy channel pole radius ->  r*cos(ang) + j r*sin(ang)
	poleAng     float64   // noisy channel pole angle in degrees
	signalFreq  int       // sine wave in Hz
	display     string    // what data to display
}

// LMS algorithm  data to create the noise canceller
type LMSAlgorithm struct {
	gain       float64   // gain of filter
	trials     int       // number of trials for algorithm
	order      int       // adaptive filter order
	samples    int       // number of samples
	samplerate int       // sample rate in Hz
	wEnsemble  []float64 // ensemble average coefficients
	wTrial     []float64 // trial coefficients
	wg         sync.WaitGroup
	toChan     chan float64 // synchronized channel to noisy channel from signal source generator
	fromChan   chan ComChan // synchronized channel to adaptive filter from noisy channel
	poleRad    float64      // noisy channel pole radius ->  r*cos(ang) + j r*sin(ang)
	poleAng    float64      // noisy channel pole angle in degrees
	signalFreq int          // pseudorandom (+/-)1 sine wave in Hz
	snr        int          // signal-to-noise ratio for noisy channel
}

var (
	timeTmpl              *template.Template
	freqTmpl              *template.Template
	lmsnoisecancellerTmpl *template.Template
	plotresponseTmpl      *template.Template
	filterSignalTmpl      *template.Template
	winType               = []string{"Bartlett", "Welch", "Hamming", "Hanning"}
)

// Bartlett window
func bartlett(n int, m int) complex128 {
	real := 1.0 - math.Abs((float64(n)-float64(m))/float64(m))
	return complex(real, 0)
}

// Welch window
func welch(n int, m int) complex128 {
	x := math.Abs((float64(n) - float64(m)) / float64(m))
	real := 1.0 - x*x
	return complex(real, 0)
}

// Hamming window
func hamming(n int, m int) complex128 {
	return complex(.54-.46*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Hanning window
func hanning(n int, m int) complex128 {
	return complex(.5-.5*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Rectangle window
func rectangle(n int, m int) complex128 {
	return 1.0
}

// init parses the html template files and is done only once at startup
func init() {
	timeTmpl = template.Must(template.ParseFiles(tmpltime))
	freqTmpl = template.Must(template.ParseFiles(tmplfrequency))
	plotresponseTmpl = template.Must(template.ParseFiles(tmplplotresponse))
	lmsnoisecancellerTmpl = template.Must(template.ParseFiles(tmplLmsNoiseCanceller))
	filterSignalTmpl = template.Must(template.ParseFiles(tmplfiltersignal))
}

// findEndpoints finds the minimum and maximum data values
func (ep *Endpoints) findEndpoints(input *bufio.Scanner, rad float64) {
	ep.xmax = -math.MaxFloat64
	ep.xmin = math.MaxFloat64
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64
	var (
		n      int = 0 // impulse response plot
		values []string
	)
	for input.Scan() {
		line := input.Text()
		// Each line has 1, 2 or 3 space-separated values, depending on if it is real or complex data:
		// real
		// time real
		// time real imaginary
		values = strings.Split(line, " ")
		var (
			x, y, t float64
			err     error
		)
		// no time data, just real value, used for impulse response plot
		if len(values) == 1 {
			if y, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}
			n++
			// real data
		} else if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				continue
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				continue
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				continue
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
		}

		// rotation
		if rad != 0.0 {
			xrot := x*math.Cos(rad) + y*math.Sin(rad)
			yrot := -x*math.Sin(rad) + y*math.Cos(rad)
			x = xrot
			y = yrot
		}

		if x > ep.xmax {
			ep.xmax = x
		}
		if x < ep.xmin {
			ep.xmin = x
		}

		if y > ep.ymax {
			ep.ymax = y
		}
		if y < ep.ymin {
			ep.ymin = y
		}
	}
	// impulse response plot
	if len(values) == 1 {
		ep.xmin = 0.0
		ep.xmax = float64(n - 1)
	}
}

// gridFill inserts the data points in the grid
func gridFill(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, rad float64, input *bufio.Scanner) error {
	var x float64 = -1
	for input.Scan() {
		line := input.Text()
		// Each line has 1, 2 or 3 space-separated values, depending on if it is real or complex data:
		// real
		// time real
		// time real imaginary
		values := strings.Split(line, " ")
		var (
			y, t float64
			err  error
		)
		// real data, no time
		if len(values) == 1 {
			x++
			if y, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}
			// real data
		} else if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				return err
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
		}

		// rotation
		if rad != 0.0 {
			xrot := x*math.Cos(rad) + y*math.Sin(rad)
			yrot := -x*math.Sin(rad) + y*math.Cos(rad)
			x = xrot
			y = yrot
		}

		// Check if inside the zoom values
		if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
			continue
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"
	}
	return nil
}

// gridFillInterp inserts the data points in the grid and draws a straight line between points
func gridFillInterp(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, rad float64, input *bufio.Scanner) error {

	var (
		x, y, t      float64
		prevX, prevY float64
		prevSet      bool = true
		err          error
	)

	const lessen = 1
	const increase = 10

	// Get first sample
	input.Scan()
	line := input.Text()
	// Each line has 1, 2 or 3 space-separated values, depending on if it is real or complex data:
	// real
	// time real
	// time real imaginary
	values := strings.Split(line, " ")
	// real data, no time
	if len(values) == 1 {
		x = 0
		if y, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}
		// real data
	} else if len(values) == 2 {
		if x, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}

		if y, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			return err
		}
		// complex data
	} else if len(values) == 3 {
		if t, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}

		if x, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			return err
		}

		if y, err = strconv.ParseFloat(values[2], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
			return err
		}

		// Calculate the modulus of the complex data becomes the y-axis
		// The time becomes the x-axis
		y = math.Sqrt(x*x + y*y)
		x = t
	}

	// rotation
	if rad != 0.0 {
		xrot := x*math.Cos(rad) + y*math.Sin(rad)
		yrot := -x*math.Sin(rad) + y*math.Cos(rad)
		x = xrot
		y = yrot
	}

	// Check if inside the zoom values
	if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
		prevSet = false
	} else {
		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		prevX = x
		prevY = y
	}

	// Scale factor to determine the number of interpolation points
	lenEP := math.Sqrt((endpoints.xmax-endpoints.xmin)*(endpoints.xmax-endpoints.xmin) +
		(endpoints.ymax-endpoints.ymin)*(endpoints.ymax-endpoints.ymin))

	// Continue with the rest of the points in the file
	for input.Scan() {
		line = input.Text()
		// Each line has 1, 2 or 3 space-separated values, depending on if it is real or complex data:
		// real
		// time real
		// time real imaginary
		values := strings.Split(line, " ")
		// real data, no time
		if len(values) == 1 {
			x++
			if y, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}
			// real data
		} else if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				return err
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
		}

		// rotation
		if rad != 0.0 {
			xrot := x*math.Cos(rad) + y*math.Sin(rad)
			yrot := -x*math.Sin(rad) + y*math.Cos(rad)
			x = xrot
			y = yrot
		}

		// Check if inside the zoom values
		if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
			continue
		} else if !prevSet {
			prevSet = true
			prevX = x
			prevY = y
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		// Interpolate the points between previous point and current point

		lenEdge := math.Sqrt((x-prevX)*(x-prevX) + (y-prevY)*(y-prevY))
		ncells := increase * int(columns*lenEdge/lenEP) / lessen // number of points to interpolate
		stepX := (x - prevX) / float64(ncells)
		stepY := (y - prevY) / float64(ncells)

		// loop to draw the points
		interpX := prevX
		interpY := prevY
		for i := 0; i < ncells; i++ {
			row := int((endpoints.ymax-interpY)*yscale + .5)
			col := int((interpX-endpoints.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
			interpX += stepX
			interpY += stepY
		}

		// Update the previous point with the current point
		prevX = x
		prevY = y
	}
	return nil
}

// processTimeDomain plots the time domain data from disk file
func processTimeDomain(w http.ResponseWriter, r *http.Request, filename string) error {

	// main data structure
	var (
		plot      PlotT
		xscale    float64
		yscale    float64
		endpoints Endpoints
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		// Mark the data x-y coordinate online at the corresponding
		// grid row/column.
		input := bufio.NewScanner(f)

		// Determine if rotate requested and perform the rotation of x and y
		rotate := r.FormValue("rotate")
		rad := 0.0
		if len(rotate) > 0 {
			deg, err := strconv.ParseFloat(rotate, 64)
			if err != nil {
				plot.Status = "Rotate degree conversion error"
				fmt.Printf("Rotate degree %v conversion error: %v\n", rotate, err)
			} else {
				rad = deg2rad * deg
			}
		}
		endpoints.findEndpoints(input, rad)

		f.Close()
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			input := bufio.NewScanner(f)

			// Determine if zoom requested and validate endpoints
			zoomxstart := r.FormValue("zoomxstart")
			zoomxend := r.FormValue("zoomxend")
			zoomystart := r.FormValue("zoomystart")
			zoomyend := r.FormValue("zoomyend")
			if len(zoomxstart) > 0 && len(zoomxend) > 0 &&
				len(zoomystart) > 0 && len(zoomyend) > 0 {
				x1, err1 := strconv.ParseFloat(zoomxstart, 64)
				x2, err2 := strconv.ParseFloat(zoomxend, 64)
				y1, err3 := strconv.ParseFloat(zoomystart, 64)
				y2, err4 := strconv.ParseFloat(zoomyend, 64)

				if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
					plot.Status = "Zoom x or y values are not numbers."
					fmt.Printf("Zoom error: x start error = %v, x end error = %v\n", err1, err2)
					fmt.Printf("Zoom error: y start error = %v, y end error = %v\n", err3, err4)
				} else {
					if (x1 < endpoints.xmin || x1 > endpoints.xmax) ||
						(x2 < endpoints.xmin || x2 > endpoints.xmax) || (x1 >= x2) {
						plot.Status = "Zoom values are not in x range."
						fmt.Printf("Zoom error: start or end value not in x range.\n")
					} else if (y1 < endpoints.ymin || y1 > endpoints.ymax) ||
						(y2 < endpoints.ymin || y2 > endpoints.ymax) || (y1 >= y2) {
						plot.Status = "Zoom values are not in y range."
						fmt.Printf("Zoom error: start or end value not in y range.\n")
						fmt.Printf("y1=%v, y2=%v, endpoints.ymin=%v, endpoints.ymax=%v\n",
							y1, y2, endpoints.ymin, endpoints.ymax)
					} else {
						// Valid Zoom endpoints, replace the previous min and max values
						endpoints.xmin = x1
						endpoints.xmax = x2
						endpoints.ymin = y1
						endpoints.ymax = y2
					}
				}
			}

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Check for interpolation and fill in the grid with the data points
			interp := r.FormValue("interpolate")
			if interp == "interpolate" {
				err = gridFillInterp(&plot, xscale, yscale, endpoints, rad, input)
			} else {
				err = gridFill(&plot, xscale, yscale, endpoints, rad, input)
			}
			if err != nil {
				return err
			}

			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := timeTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// processFrequencyDomain calculates the power spectral density (PSD) and plots it
func processFrequencyDomain(w http.ResponseWriter, r *http.Request, filename string) error {
	// Use complex128 for FFT computation
	// Get the number of complex samples nn, open file and count lines, close the file

	var (
		plot          PlotT // main data structure to execute with parsed html template
		endpoints     Endpoints
		N             int                                                        //  complex FFT size
		nn            int                                                        // number of complex samples in the data file
		K             int                                                        //  number of segments used in PSD with 50% overlap
		m             int                                                        // complex segment size
		win           string                                                     // FFT window type
		window        = make(map[string]Window, len(winType))                    // map of window functions
		sumWindow     float64                                                    // sum of squared window values for normalization
		normalizerPSD float64                                                    // normalizer for PSD
		PSD           []float64                                                  // power spectral density
		psdMax        float64                                 = -math.MaxFloat64 // maximum PSD value
		psdMin        float64                                 = math.MaxFloat64  // minimum PSD value
		xscale        float64                                                    // data to grid in x direction
		yscale        float64                                                    // data to grid in y direction
		samplingRate  float64                                                    // sampling rate in Hz
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Put the window functions in the map
	window["Bartlett"] = bartlett
	window["Welch"] = welch
	window["Hamming"] = hamming
	window["Hanning"] = hanning
	window["Rectangle"] = rectangle

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		input := bufio.NewScanner(f)
		// Number of real or complex data samples
		for input.Scan() {
			line := input.Text()
			if len(line) > 0 {
				nn++
			}
		}
		fmt.Printf("Data file %s has %d samples\n", filename, nn)
		// make even number of samples so if segments = 1, we won't
		// do the last FFT with one sample
		if nn%2 == 1 {
			nn++
		}

		f.Close()

		// Get number of segments from HTML form
		// Number of segments to average the periodograms to reduce the variance
		tmp := r.FormValue("fftsegments")
		if len(tmp) == 0 {
			return fmt.Errorf("enter number of FFT segments")
		}
		K, err = strconv.Atoi(tmp)
		if err != nil {
			fmt.Printf("FFT segments string convert error: %v\n", err)
			return fmt.Errorf("fft segments string convert error: %v", err)
		}

		// Require 1 <= K <= 20
		if K < 1 {
			K = 1
		} else if K > 20 {
			K = 20
		}

		// segment size complex samples
		m = nn / (K + 1)

		// Get window type:  Bartlett, Welch, Hanning, Hamming, etc
		// Multiply the samples by the window to reduce spectral leakage
		// caused by high sidelobes in rectangular window
		win = r.FormValue("fftwindow")
		if len(win) == 0 {
			return fmt.Errorf("enter FFT window type")
		}
		w, ok := window[win]
		if !ok {
			fmt.Printf("Invalid FFT window type: %v\n", win)
			return fmt.Errorf("invalid FFT window type: %v", win)
		}
		// sum the window values for PSD normalization due to windowing
		for i := 0; i < 2*m; i++ {
			x := cmplx.Abs(w(i, m))
			sumWindow += x * x
		}
		fmt.Printf("%s window sum = %.2f\n", win, sumWindow)

		// Get FFT size from HTML form
		// Check FFT Size >= 2*m, using 50%  overlap of segments
		// Check FFT Size is a power of 2:  2^n
		tmp = r.FormValue("fftsize")
		if len(tmp) == 0 {
			return fmt.Errorf("enter FFT size")
		}
		N, err = strconv.Atoi(tmp)
		if err != nil {
			return fmt.Errorf("fft size string convert error: %v", err)
		}

		if N < rows {
			fmt.Printf("FFT size < %d\n", rows)
			N = rows
		} else if N > rows*rows {
			fmt.Printf("FFT size > %d\n", rows*rows)
			N = rows * rows
		}
		// This rounds up to nearest FFT size that is a power of 2
		N = int(math.Exp2(float64(int(math.Log2(float64(N)) + .5))))
		fmt.Printf("N=%v\n", N)

		if N < 2*m {
			fmt.Printf("FFT Size %d not greater than 2%d\n", N, 2*m)
			return fmt.Errorf("fft Size %d not greater than 2*%d", N, 2*m)
		}

		// Power Spectral Density, PSD[N/2] is the Nyquist critical frequency
		// It is the (sampling frequency)/2, the highest non-aliased frequency
		PSD = make([]float64, N/2+1)

		// Reopen the data file
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			bufm := make([]complex128, m)
			bufN := make([]complex128, N)
			input := bufio.NewScanner(f)
			// Read in initial m samples into buf[m] to start the processing loop
			diff := 0.0
			prev := 0.0
			for k := 0; k < m; k++ {
				input.Scan()
				line := input.Text()
				// Each line has 1, 2 or 3 space-separated values, depending on if it is real or complex data:
				// real
				// time real
				// time real imaginary
				values := strings.Split(line, " ")
				if len(values) == 1 {
					var r float64
					if r, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if k != 0 {
						diff += .5
					}
					bufm[k] = complex(r, 0)
					// real data
				}
				if len(values) == 2 {
					// time real, calculate the sampling rate from the time steps
					var t, r float64

					if t, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if r, err = strconv.ParseFloat(values[1], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
						continue
					}

					if k == 0 {
						prev = t
					} else {
						diff += t - prev
						prev = t
					}
					bufm[k] = complex(r, 0)

					// complex data
				} else if len(values) == 3 {
					// time real imag
					var t, r, i float64

					// calculate the sampling rate from the time steps
					if t, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if r, err = strconv.ParseFloat(values[1], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
						continue
					}

					if i, err = strconv.ParseFloat(values[2], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
						continue
					}

					if k == 0 {
						prev = t
					} else {
						diff += t - prev
						prev = t
					}
					bufm[k] = complex(r, i)
				}
			}

			// Average the time steps and invert to get the sampling rate
			samplingRate = 1.0 / (diff / float64(m-1))
			fmt.Printf("sampling rate = %.2f\n", samplingRate)

			scanOK := true
			// loop over the rest of the file, reading in m samples at a time until EOF
			for {
				// Put the previous m samples in the front of the buffer N to
				// overlap segments
				copy(bufN, bufm)
				// Put the next m samples in back of the previous m samples
				kk := 0
				for k := 0; k < m; k++ {
					scanOK = input.Scan()
					if !scanOK {
						break
					}
					line := input.Text()
					// Each line has 1 - 3 values: [time], real, [imag]
					values := strings.Split(line, " ")
					if len(values) == 1 {
						var r float64

						if r, err = strconv.ParseFloat(values[0], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
							continue
						}
						// real data, no time or imag
						bufm[k] = complex(r, 0)
						// real data
					} else if len(values) == 2 {
						// time real, but don't need the time
						var r float64

						if r, err = strconv.ParseFloat(values[1], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
							continue
						}
						// real data, imaginary component is zero
						bufm[k] = complex(r, 0)

						// complex data
					} else if len(values) == 3 {
						// time real imag
						var r, i float64
						// Don't need the time
						if r, err = strconv.ParseFloat(values[1], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
							continue
						}

						if i, err = strconv.ParseFloat(values[2], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
							continue
						}

						bufm[k] = complex(r, i)
					}
					kk++
				}
				// Check for the normal EOF and the abnormal scan error
				// EOF does not give an error but is considered normal termination
				if !scanOK {
					if input.Err() != nil {
						fmt.Printf("Data file scan error: %v\n", input.Err().Error())
						return fmt.Errorf("data file scan error: %v", input.Err())
					}
				}
				// Put the next kk samples in back of the previous m
				copy(bufN[m:], bufm[:kk])

				// window the m + kk samples with chosen window
				for i := 0; i < m+kk; i++ {
					bufN[i] *= w(i, m)
				}

				// zero-pad N-m-kk samples in buf[N]
				for i := m + kk; i < N; i++ {
					bufN[i] = 0
				}

				// Perform N-point complex FFT and add squares to previous values in PSD[N/2+1]
				fourierN := fft.FFT(bufN)
				x := cmplx.Abs(fourierN[0])
				PSD[0] += x * x
				for i := 1; i < N/2; i++ {
					// Use positive and negative frequencies -> bufN[N-i] = bufN[-i]
					xi := cmplx.Abs(fourierN[i])
					xNi := cmplx.Abs(fourierN[N-i])
					PSD[i] += xi*xi + xNi*xNi
				}
				x = cmplx.Abs(fourierN[N/2])
				PSD[N/2] += x * x

				// part of K*Sum(w[i]*w[i]) PSD normalizer
				normalizerPSD += sumWindow

				// EOF reached
				if !scanOK {
					break
				}
			} // K segments done

			// Normalize the PSD using K*Sum(w[i]*w[i])
			// Use log plot for wide dynamic range
			if r.FormValue("plottype") == "linear" {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
				// log10 in dB
			} else {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					PSD[i] = 10.0 * math.Log10(PSD[i])
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
			}

			endpoints.xmin = 0
			endpoints.xmax = float64(N / 2) // equivalent to Nyquist critical frequency
			endpoints.ymin = psdMin
			endpoints.ymax = psdMax

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Store the PSD in the plot Grid
			for bin, pow := range PSD {
				row := int((endpoints.ymax-pow)*yscale + .5)
				col := int((float64(bin)-endpoints.xmin)*xscale + .5)
				plot.Grid[row*columns+col] = "online"
			}

			// Store in the form:  FFT Size, window type, number of samples nn, K segments, sampling frequency
			// Plot the PSD N/2 float64 values, execute the data on the plotfrequency.html template

			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Apply the  sampling rate in Hz to the x-axis using a scale factor
	// Convert the fft size to sampling rate/2, the Nyquist critical frequency
	sf := 0.5 * samplingRate / endpoints.xmax

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	format := "%.0f"
	if incr*sf < 1.0 {
		format = "%.2f"
	}
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf(format, x*sf)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Insert frequency domain parameters in the form
	plot.SampleFreq = fmt.Sprintf("%.0f", samplingRate)
	plot.FFTSegments = strconv.Itoa(K)
	plot.FFTSize = strconv.Itoa(N)
	plot.Samples = strconv.Itoa(nn)

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := freqTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// handlePlotResponse displays the impulse or frequency response of the
// LMS NoiseCanceller
func handlePlotResponse(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot PlotT
		err  error = nil
	)

	filename := r.FormValue("filename")
	// choose time or frequency domain processing
	if len(filename) > 0 {

		domain := r.FormValue("domain")
		switch domain {
		case "time":
			err = processTimeDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		case "frequency":
			err = processFrequencyDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		default:
			plot.Status = fmt.Sprintf("Invalid domain choice: %s", domain)
		}

		if err != nil {

			// Write to HTTP using template and grid``
			if err := plotresponseTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		plot.Status = "Missing filename to plot"
		// Write to HTTP using template and grid``
		if err := plotresponseTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// generateSignal generates a new set of input data and sends it to the toChannel
func (lms *LMSAlgorithm) generateSignal() error {

	// generate a new set of input data for this trial

	// increment wg
	lms.wg.Add(1)

	// launch a goroutine to generate input samples consisting of
	// a sine wave at unit amplitude and random phase
	go func() {
		defer func() {
			lms.wg.Done()
			close(lms.toChan)
		}()

		phase := twoPi * rand.Float64()
		omega := twoPi * float64(lms.signalFreq)
		delta := 1.0 / float64(lms.samplerate)
		// Send a sine wave with random phase to channel via the toChannel
		for i := 0; i < lms.samples; i++ {
			// new sample, randomly chosen
			for t := 0.0; t < float64(lms.samples)*delta; t += delta {
				// Send the sample through the channel
				lms.toChan <- math.Sin(omega*t + phase)
			}
		}
	}()
	return nil
}

// chanFilter adds Gaussian noise at the SNR to primary input
// and generates correlated noise for the reference input
func (lms *LMSAlgorithm) chanFilter() error {

	lms.wg.Add(1)
	go func() {

		defer func() {
			lms.wg.Done()
			close(lms.fromChan)
		}()

		// previous IIR outputs
		r := make([]float64, 2)

		// Convert degrees to radians
		theta := lms.poleAng * math.Pi / 180.0

		// IIR coefficients for correlated noise
		h := make([]float64, 2)
		h[0] = 2.0 * math.Cos(theta) * lms.poleRad
		h[1] = -lms.poleRad * lms.poleRad

		// Calculate the noise standard deviation using the SNR and signal amplitude=1
		noiseSD := math.Sqrt(math.Pow(10.0, -float64(2*lms.snr)/10.0))

		// range over the to channel to obtain the pure sine wave
		for sig := range lms.toChan {
			// generate a normal rv at specified standard deviation
			norm := noiseSD * rand.NormFloat64()
			// correlated noise for the reference LMS input
			correlated := norm
			for i, val := range r {
				correlated += h[i] * val
			}

			// Send output to next stage via the from channel
			// "desired" is the signal + normal noise, "in" is the correlated noise
			lms.fromChan <- ComChan{desired: sig + norm, in: correlated}
			// Update the recursion for the correlated noise
			r[0], r[1] = norm, r[0]
		}
	}()
	return nil
}

// runLms updates the adaptive weights using the LMS algorithm
func (lms *LMSAlgorithm) runLms() error {

	// increment wg
	lms.wg.Add(1)

	// launch a goroutine to generate the adaptive filter
	go func() {
		defer lms.wg.Done()
		L := lms.order + 1
		// holds the previous inputs
		x := make([]float64, L)
		i := 0
		// range over the channel containing the noisy signal
		for d := range lms.fromChan {
			x[i] = d.in
			y := 0.0
			k := i
			for j := 0; j < L; j++ {
				y += lms.wTrial[j] * x[k]
				k = (k - 1) % L
				if k < 0 {
					k = L + k
				}
			}
			// calculate the error
			e := d.desired - y
			k = i
			// update the trial adaptive weights
			for j := 0; j < L; j++ {
				// normalize the update by the average square-
				// input multiplied by the number of filter coefficients
				lms.wTrial[j] += lms.gain * e * x[k] / float64(L)
				// update the ensemble weight
				lms.wEnsemble[j] += lms.wTrial[j]
				k = (k - 1) % L
				if k < 0 {
					k = L + k
				}
			}
			// Increment the current input index for x slice
			i = (i + 1) % L
		}
	}()

	return nil
}

// handleLmsNoiseCanceller creates a Noise Canceller using the LMS algorithm
func handleLmsNoiseCanceller(w http.ResponseWriter, r *http.Request) {
	var plot PlotT

	// Get the number of samples to generate and the sample rate
	// The number of samples determines the number of iterations of the LMS algorithm
	samplestxt := r.FormValue("samples")
	sampleratetxt := r.FormValue("samplerate")
	// choose time or frequency domain processing
	if len(samplestxt) > 0 && len(sampleratetxt) > 0 {

		samples, err := strconv.Atoi(samplestxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Samples conversion to int error: %v", err.Error())
			fmt.Printf("Samples conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		samplerate, err := strconv.Atoi(sampleratetxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Sample rate conversion to int error: %v", err.Error())
			fmt.Printf("Sample rate conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt := r.FormValue("signalfreq")
		sigfreq, err := strconv.Atoi(txt)
		if err != nil {
			fmt.Printf("Signal frequency conversion error: %v\n", err)
			plot.Status = fmt.Sprintf("Signal frequency conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Verify signal frequency is less than Nyquist frequency = sample rate / 2
		if sigfreq > samplerate/2 {
			fmt.Printf("Signal frequency %v is greater than Nyquist frequency %v\n", sigfreq, samplerate/2)
			plot.Status = fmt.Sprintf("Signal frequency %v is greater than Nyquist frequency %v\n", sigfreq, samplerate/2)
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// adaptive filter order is the number of past samples
		// filter length = order + 1
		txt = r.FormValue("filtorder")
		order, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("Filter order conversion to int error: %v", err.Error())
			fmt.Printf("Filter order conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}
		// Make filter length odd, which means filter order is even;
		// if order is 7, it is changed to 8, and the length becomes 9.
		order = (order + 1) / 2 * 2

		// gain factor that regulates the speed and stability of adaption
		txt = r.FormValue("gain")
		gain, err := strconv.ParseFloat(txt, 64)
		if err != nil {
			plot.Status = fmt.Sprintf("Gain conversion to float64 error: %v", err.Error())
			fmt.Printf("Gain conversion to float64 error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// number of trials to perform the LMS algorithm to get ensemble average of the
		// weights
		txt = r.FormValue("trials")
		trials, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("Trials conversion to int error: %v", err.Error())
			fmt.Printf("Trials conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt = r.FormValue("snr")
		snr, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("SNR conversion to int error: %v", err.Error())
			fmt.Printf("SNR conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt = r.FormValue("poleang")
		poleAng, err := strconv.ParseFloat(txt, 64)
		if err != nil {
			plot.Status = fmt.Sprintf("Pole angle conversion to float64 error: %v", err.Error())
			fmt.Printf("Pole angle conversion to float64 error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Verify pole angle
		if poleAng < 0 || poleAng > 180 {
			plot.Status = fmt.Sprintf("Pole angle [%v] not between 0 and 180.\n", poleAng)
			fmt.Printf("Pole angle [%v] not between 0 and 180.\n", poleAng)
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt = r.FormValue("polerad")
		poleRad, err := strconv.ParseFloat(txt, 64)
		if err != nil {
			plot.Status = fmt.Sprintf("Pole radius conversion to float64 error: %v", err.Error())
			fmt.Printf("Pole radius conversion to float64 error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Verify pole radius
		if poleRad < .01 || poleRad > .99 {
			plot.Status = fmt.Sprintf("Pole radius [%v] not between .01 and .99.\n", poleRad)
			fmt.Printf("Pole radius [%v] not between .01 and .99.\n", poleRad)
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Construct object to hold LMS algorithm parameters
		lmsNoiseCanceller := LMSAlgorithm{
			gain:       gain,
			trials:     trials,
			order:      order,
			samples:    samples,
			samplerate: samplerate,
			wEnsemble:  make([]float64, order+1),
			wTrial:     make([]float64, order+1),
			snr:        snr,
			poleRad:    poleRad,
			poleAng:    poleAng,
			signalFreq: sigfreq,
		}

		// Run the Least-Mean-Square (LMS) algorithm to create the adaptive filter
		// Loop over the trials to generate the ensemble of filters which is averaged.
		for i := 0; i < lmsNoiseCanceller.trials; i++ {

			// Create new channels each trial since they are closed after each trial
			lmsNoiseCanceller.fromChan = make(chan ComChan)
			lmsNoiseCanceller.toChan = make(chan float64)

			// Generate a new set of input data consisting of a sine wave at unit amplitude
			err = lmsNoiseCanceller.generateSignal()
			if err != nil {
				plot.Status = fmt.Sprintf("generateSignal error: %v", err.Error())
				fmt.Printf("generateSignal error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Run the signal through the noisy channel with added Gaussian noise
			err = lmsNoiseCanceller.chanFilter()
			if err != nil {
				plot.Status = fmt.Sprintf("chanFilter error: %v", err.Error())
				fmt.Printf("chanFilter error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Run the LMS algorithm to find the adaptive filter coefficients
			err = lmsNoiseCanceller.runLms()
			if err != nil {
				plot.Status = fmt.Sprintf("generateLMSData error: %v", err.Error())
				fmt.Printf("runLmsNoiseCanceller error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			// Wait for this trial to finish
			lmsNoiseCanceller.wg.Wait()
		}

		// Save ensemble filter coefficients to disk (averaged coefficients)
		file, err := os.Create(path.Join(dataDir, lmsFilterNoiseCanceller))
		if err != nil {
			plot.Status = fmt.Sprintf("create %v error: %v", path.Join(dataDir, lmsFilterNoiseCanceller), err.Error())
			fmt.Printf("create %v error: %v\n", path.Join(dataDir, lmsFilterNoiseCanceller), err.Error())
			// Write to HTTP using template and grid
			if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}
		defer file.Close()

		// Average the weight ensemble to reduce variance
		for i := 0; i < lmsNoiseCanceller.order+1; i++ {
			fmt.Fprintf(file, "%v\n", lmsNoiseCanceller.wEnsemble[i]/(float64(trials*samples)))
		}

		plot.Status = fmt.Sprintf("Adaptive filter weights '%s' written to the data folder", lmsFilterNoiseCanceller)
		// Write to HTTP using template and grid
		if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}

	} else {
		plot.Status = "Enter samples, sample rate, SNR, etc."
		// Write to HTTP using template and grid
		if err := lmsnoisecancellerTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// gridFillInterp inserts the data points in the grid and draws a straight line between points
func (fs *FilterSignal) gridFillInterp(plot *PlotT) error {
	var (
		x            float64 = fs.firstSampleTime
		y            float64 = 0.0
		prevX, prevY float64
		err          error
		xscale       float64
		yscale       float64
		input        *bufio.Scanner
		timeStep     float64 = 1.0 / float64(fs.sampleFreq)
	)

	// Mark the data x-y coordinate online at the corresponding
	// grid row/column.

	fs.xmin = fs.firstSampleTime
	fs.xmax = fs.lastSampleTime

	// Calculate scale factors for x and y
	xscale = (columns - 1) / (fs.xmax - fs.xmin)
	yscale = (rows - 1) / (fs.ymax - fs.ymin)

	f, err := os.Open(path.Join(dataDir, signal))
	if err != nil {
		fmt.Printf("Error opening %s: %v\n", signal, err.Error())
		return err
	}
	defer f.Close()
	input = bufio.NewScanner(f)

	// Get first sample
	input.Scan()
	value := input.Text()

	if y, err = strconv.ParseFloat(value, 64); err != nil {
		fmt.Printf("gridFillInterp first sample string %s conversion to float error: %v\n", value, err)
		return err
	}

	plot.Grid = make([]string, rows*columns)

	// This cell location (row,col) is on the line
	row := int((fs.ymax-y)*yscale + .5)
	col := int((x-fs.xmin)*xscale + .5)
	plot.Grid[row*columns+col] = "online"

	prevX = x
	prevY = y

	// Scale factor to determine the number of interpolation points
	lenEPy := fs.ymax - fs.ymin
	lenEPx := fs.xmax - fs.xmin

	// Continue with the rest of the points in the file
	for input.Scan() {
		x += timeStep
		value = input.Text()
		if y, err = strconv.ParseFloat(value, 64); err != nil {
			fmt.Printf("gridFillInterp the rest of file string %s conversion to float error: %v\n", value, err)
			return err
		}

		// This cell location (row,col) is on the line
		row := int((fs.ymax-y)*yscale + .5)
		col := int((x-fs.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		// Interpolate the points between previous point and current point

		/* lenEdge := math.Sqrt((x-prevX)*(x-prevX) + (y-prevY)*(y-prevY)) */
		lenEdgeX := math.Abs((x - prevX))
		lenEdgeY := math.Abs(y - prevY)
		ncellsX := int(columns * lenEdgeX / lenEPx) // number of points to interpolate in x-dim
		ncellsY := int(rows * lenEdgeY / lenEPy)    // number of points to interpolate in y-dim
		// Choose the biggest
		ncells := ncellsX
		if ncellsY > ncells {
			ncells = ncellsY
		}

		stepX := (x - prevX) / float64(ncells)
		stepY := (y - prevY) / float64(ncells)

		// loop to draw the points
		interpX := prevX
		interpY := prevY
		for i := 0; i < ncells; i++ {
			row := int((fs.ymax-interpY)*yscale + .5)
			col := int((interpX-fs.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
			interpX += stepX
			interpY += stepY
		}

		// Update the previous point with the current point
		prevX = x
		prevY = y
	}
	return nil
}

// generateSignal creates a since  wave at unit amplitude
func (fs *FilterSignal) generateSignal() error {

	// determine where to send the signal

	// save the signal to a disk file
	if fs.display == "channelin" {
		f, _ := os.Create(path.Join(dataDir, signal))
		defer f.Close()
		r := fs.sampleFreq / fs.signalFreq
		current := 1.0

		// Save a signal sample to a disk file
		for i := 0; i < fs.samples; i++ {
			// new sample, randomly chosen
			if i%r == 0 {
				if rand.Float64() < 0.5 {
					current = -1.0
				} else {
					current = 1.0
				}
			}
			fmt.Fprintf(f, "%f\n", current)

			// find min/max of the signal as we go
			if current < fs.ymin {
				fs.ymin = current
			}
			if current > fs.ymax {
				fs.ymax = current
			}
		}
		// send the signal to the next stage
	} else {
		// increment wg
		fs.wg.Add(1)

		// launch a goroutine to generate input samples consisting of
		// a sine wave at unit amplitude
		go func() {
			defer func() {
				fs.wg.Done()
				close(fs.toChan)
			}()

			r := fs.sampleFreq / fs.signalFreq
			current := 1.0

			// Send a signal sample to noisy channel via the toChannel
			for i := 0; i < fs.samples; i++ {
				// new sample, randomly chosen
				if i%r == 0 {
					if rand.Float64() < 0.5 {
						current = -1.0
					} else {
						current = 1.0
					}
				}
				fs.toChan <- current
			}
		}()
	}
	return nil
}

// noiseCancellerFilter filters the signal with the LMS noiseCanceller
func (fs *FilterSignal) noiseCancellerFilter() error {
	// determine where to send the signal

	// this stage not used for these display choices
	if fs.display == "channelin" || fs.display == "channelout" {
		return nil
		// save the output of this stage to disk file
	} else {

		// Get noise canceller coefficients from running of LMS algorithm
		f, err := os.Open(path.Join(dataDir, fs.filterfile))
		if err != nil {
			fmt.Printf("Error opening %s: %v\n", fs.filterfile, err.Error())
			return err
		}
		input := bufio.NewScanner(f)
		for input.Scan() {
			line := input.Text()
			wt, err := strconv.ParseFloat(line, 64)
			if err != nil {
				fmt.Printf("strconv.ParseFloat of %s in %s error: %v\n", line, fs.filterfile, err)
				return fmt.Errorf("filter coefficient conversion error: %v", err)
			}
			fs.filterCoeff = append(fs.filterCoeff, wt)
		}
		f.Close()

		// Save output of this stage to disk file
		f, _ = os.Create(path.Join(dataDir, signal))
		// increment wg
		fs.wg.Add(1)
		// launch a goroutine to filter the signal
		go func() {
			defer func() {
				fs.wg.Done()
				f.Close()
			}()

			L := len(fs.filterCoeff)
			// x slice holds the previous inputs so we can do the convolution sum
			x := make([]float64, L)
			i := 0
			// range over the channel containing the noisy signal
			for d := range fs.fromChan {
				x[i] = d
				y := 0.0
				k := i
				for j := 0; j < L; j++ {
					y += fs.filterCoeff[j] * x[k]
					k = (k - 1) % L
					if k < 0 {
						k = L + k
					}
				}

				// Send output to disk file
				fmt.Fprintf(f, "%f\n", y)

				// find min/max of the signal as we go
				if y < fs.ymin {
					fs.ymin = y
				}
				if y > fs.ymax {
					fs.ymax = y
				}
				// Increment the current input index for x slice
				i = (i + 1) % L

			}
		}()
	}
	return nil
}

// chanFilter implements the noisy channel
func (fs *FilterSignal) chanFilter() error {
	// determine where to send the signal
	// this stage not used for this display choice
	if fs.display == "channelin" {
		return nil
		// save the output of this stage to disk file
	} else if fs.display == "channelout" {
		f, _ := os.Create(path.Join(dataDir, signal))
		fs.wg.Add(1)
		go func() {

			defer func() {
				fs.wg.Done()
				defer f.Close()
			}()

			// previous IIR outputs
			r := make([]float64, 2)

			// Convert degrees to radians
			theta := fs.poleAng * math.Pi / 180.0

			// IIR coefficients for channel impulse response
			h := make([]float64, 2)
			h[0] = 2.0 * math.Cos(theta) * fs.poleRad
			h[1] = -fs.poleRad * fs.poleRad

			// Calculate the noise standard deviation using the SNR and signal amplitude=1
			noiseSD := math.Sqrt(math.Pow(10.0, -float64(fs.snr)/10.0))

			// range over the input channel
			for d := range fs.toChan {
				temp := d + noiseSD*rand.NormFloat64()
				for i, val := range r {
					temp += h[i] * val
				}

				// Send output to disk file
				fmt.Fprintf(f, "%v\n", temp)

				// find min/max of the signal as we go
				if temp < fs.ymin {
					fs.ymin = temp
				}
				if temp > fs.ymax {
					fs.ymax = temp
				}

				// Update the recursion
				r[0], r[1] = d, r[0]
			}
		}()
		// send the output to the next stage
	} else {
		fs.wg.Add(1)
		go func() {

			defer func() {
				fs.wg.Done()
				close(fs.fromChan)
			}()

			// previous IIR outputs
			r := make([]float64, 2)

			// Convert degrees to radians
			theta := fs.poleAng * math.Pi / 180.0

			// IIR coefficients for channel impulse response
			h := make([]float64, 2)
			h[0] = 2.0 * math.Cos(theta) * fs.poleRad
			h[1] = -fs.poleRad * fs.poleRad

			// Calculate the noise standard deviation using the SNR and signal amplitude=1
			noiseSD := math.Sqrt(math.Pow(10.0, -float64(fs.snr)/10.0))

			// range over the input channel
			for d := range fs.toChan {
				temp := d + noiseSD*rand.NormFloat64()
				for i, val := range r {
					temp += h[i] * val
				}

				// Send output to next stage via the output channel
				fs.fromChan <- temp
				// Update the recursion
				r[0], r[1] = d, r[0]
			}
		}()
	}
	return nil
}

// label the plot and execute the PlotT on the HTML template
func (fs *FilterSignal) labelExec(w http.ResponseWriter, plot *PlotT) {

	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Construct x-axis labels
	incr := (fs.xmax - fs.xmin) / (xlabels - 1)
	x := fs.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (fs.ymax - fs.ymin) / (ylabels - 1)
	y := fs.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Fill in the form fields
	plot.Samples = strconv.Itoa(fs.samples)
	plot.SampleFreq = strconv.Itoa(fs.sampleFreq)
	plot.SNR = strconv.Itoa(fs.snr)
	plot.Filename = path.Base(fs.filterfile)
	plot.PoleAng = fmt.Sprintf("%f", fs.poleAng)
	plot.PoleRad = fmt.Sprintf("%f", fs.poleRad)
	plot.SigFreq = strconv.Itoa(fs.signalFreq)

	if len(plot.Status) == 0 {
		plot.Status = fmt.Sprintf("Signal consisting of a sine wave in noisy channel"+
			" with Gaussian white noise was filtered with %s", plot.Filename)
	}

	// Write to HTTP using template and grid
	if err := filterSignalTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}
}

// handleFilter filters a noisy signal using the LMS Noise Canceller
func handleFilterSignal(w http.ResponseWriter, r *http.Request) {
	var (
		plot       PlotT
		snr        int
		poleAng    float64
		poleRad    float64
		filterfile string
	)

	// need number of samples and sample frequency to continue
	temp := r.FormValue("samples")
	if len(temp) > 0 {
		samples, err := strconv.Atoi(temp)
		if err != nil {
			plot.Status = fmt.Sprintf("Samples conversion to int error: %v", err.Error())
			fmt.Printf("Samples conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		temp = r.FormValue("samplefreq")
		sampfreq, err := strconv.Atoi(temp)
		if err != nil {
			fmt.Printf("Sample frequency conversion error: %v\n", err)
			plot.Status = fmt.Sprintf("Sample frequency conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		temp = r.FormValue("signalfreq")
		sigfreq, err := strconv.Atoi(temp)
		if err != nil {
			fmt.Printf("Signal frequency conversion error: %v\n", err)
			plot.Status = fmt.Sprintf("Signal frequency conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Determine what else is needed to display the data
		display := r.FormValue("display")
		// Need pole data and snr
		if display == "channelout" || display == "noisecancellerout" {
			txt := r.FormValue("snr")
			snr, err = strconv.Atoi(txt)
			if err != nil {
				plot.Status = fmt.Sprintf("SNR conversion to int error: %v", err.Error())
				fmt.Printf("SNR conversion to int error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			txt = r.FormValue("poleang")
			poleAng, err = strconv.ParseFloat(txt, 64)
			if err != nil {
				plot.Status = fmt.Sprintf("Pole angle conversion to int error: %v", err.Error())
				fmt.Printf("Pole angle conversion to int error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Verify pole angle
			if poleAng < 0 || poleAng > 180 {
				plot.Status = fmt.Sprintf("Pole angle [%v] not between 0 and 180.\n", poleAng)
				fmt.Printf("Pole angle [%v] not between 0 and 180.\n", poleAng)
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			txt = r.FormValue("polerad")
			poleRad, err = strconv.ParseFloat(txt, 64)
			if err != nil {
				plot.Status = fmt.Sprintf("Pole radius conversion to float64 error: %v", err.Error())
				fmt.Printf("Pole radius conversion to float64 error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

		}
		// Need adaptive filter file
		if display == "noisecancellerout" {
			filterfile = r.FormValue("filterfile")
			if len(filterfile) == 0 {
				plot.Status = "Need Adaptive Filter File"
				fmt.Printf("Need Adaptive Filter File\n")
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
		}

		// Get FilterState from previous submit and store in FilterSignal
		var filterState FilterState
		// This file exists only after all the samples are generated and filtered
		f, err := os.Open(path.Join(dataDir, stateFile))
		if err == nil {
			defer f.Close()
			input := bufio.NewScanner(f)
			input.Scan()
			line := input.Text()
			// Last sample time from previous submit
			lst, err := strconv.ParseFloat(line, 64)
			if err != nil {
				fmt.Printf("From %s, sample time conversion error: %v\n", stateFile, err)
				plot.Status = fmt.Sprintf("Sample time conversion to float error: %v", err.Error())
				// Write to HTTP using template and grid
				if err := filterSignalTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			filterState.firstSampleTime = lst
			filterState.lastFiltered = make([]float64, 0)
			// Get the last incomplete filtered outputs from previous submit
			for input.Scan() {
				line := input.Text()
				fltOut, err := strconv.ParseFloat(line, 64)
				if err != nil {
					fmt.Printf("Sample time conversion error: %v\n", err)
					plot.Status = fmt.Sprintf("From %s, filtered output conversion to float error: %v", stateFile, err.Error())
					// Write to HTTP using template and grid
					if err := filterSignalTmpl.Execute(w, plot); err != nil {
						log.Fatalf("Write to HTTP output using template with error: %v\n", err)
					}
					return
				}
				filterState.lastFiltered = append(filterState.lastFiltered, fltOut)
			}
		} else {
			filterState = FilterState{
				firstSampleTime: 0.0,
				lastSampleTime:  float64(samples) / float64(sampfreq),
				lastFiltered:    make([]float64, 0),
			}
		}

		// create FilterSignal instance fs
		fs := FilterSignal{
			samples:     samples,
			sampleFreq:  sampfreq,
			signalFreq:  sigfreq,
			FilterState: filterState,
			filterCoeff: make([]float64, 0),
			filterfile:  filterfile,
			poleRad:     poleRad,
			poleAng:     poleAng,
			snr:         snr,
			fromChan:    make(chan float64),
			toChan:      make(chan float64),
			display:     display,
		}

		// Generate input data consisting of a sine wave at unit amplitude
		err = fs.generateSignal()
		if err != nil {
			plot.Status = fmt.Sprintf("generateSignal error: %v", err.Error())
			fmt.Printf("generateSignal error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Run the signal through the noisy channel with white Gaussian noise
		err = fs.chanFilter()
		if err != nil {
			plot.Status = fmt.Sprintf("chanFilter error: %v", err.Error())
			fmt.Printf("chanFilter error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Run the noisy channel output through the adaptive Noise Canceller
		err = fs.noiseCancellerFilter()
		if err != nil {
			plot.Status = fmt.Sprintf("noiseCancellerFilter error: %v", err.Error())
			fmt.Printf("noiseCancellerFilter error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// wait for the generator and filter goroutines to complete
		fs.wg.Wait()

		if err != nil {
			plot.Status = err.Error()
		}

		// Fill in the PlotT grid with the signal time vs amplitude
		err = fs.gridFillInterp(&plot)
		if err != nil {
			plot.Status = err.Error()
			fmt.Printf("gridFillInterp error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := filterSignalTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		//  generate x-labels, ylabels, status in PlotT and execute the data on the HTML template
		fs.labelExec(w, &plot)

	} else {
		// delete previous state file if initial connection (not a submit)
		if err := os.Remove(path.Join(dataDir, stateFile)); err != nil {
			// ignore error if file is not present
			fmt.Printf("Remove file %s error: %v\n", path.Join(dataDir, stateFile), err)
		}
		plot.Status = "Enter samples, sample frequency, SNR, etc."
		if err := filterSignalTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// executive program
func main() {

	// Setup http server with handler for filtering a noisy signal using the LMS Noise Canceller
	http.HandleFunc(patternFilterSignal, handleFilterSignal)

	// Setup http server with handler for creating a Noise Canceller using LMS
	http.HandleFunc(patternLmsNoiseCanceller, handleLmsNoiseCanceller)

	// Setup http server with handler for plotting impulse or frequency responses for the LMS Noise Canceller
	http.HandleFunc(patternPlotResponse, handlePlotResponse)

	fmt.Printf("LMS Noise Canceller Server listening on %v.\n", addr)

	http.ListenAndServe(addr, nil)
}
