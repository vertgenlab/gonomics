package main

import (
	"flag"
	"fmt"
	"log"
	"strconv"
)

func sizeToRatio(size float64) float64 {
	// TODO: where did this function come from?
	m := (1.0 - 0.4) / (200 - 800)
	b := 0.4 - m*800
	return m*size + b
}

func usage() {
	fmt.Print(
		"calcAmpureRatio - calculate a protocol to AMPure beads for a double size selection\n" +
			"Usage:\n" +
			" calcAmpureRatio minimumFragmentSize(bp) maximumFragmentSize(bp) startingSampleVolume(ul)\n" +
			"\n")
	flag.PrintDefaults()
}

func calcAmpureRatio(lowSize float64, highSize float64, volume float64) {
	lowRatio := sizeToRatio(lowSize)
	highRatio := sizeToRatio(highSize)

	fmt.Print(
		"Using AMPure beads to obtain DNA fragments in the size range of " + fmt.Sprintf("%f", lowSize) + " to " + fmt.Sprintf("%.2f", highSize) + " when starting with " + fmt.Sprintf("%.2f", volume) + " ul of sample\n" +
			"\n" +
			"For " + fmt.Sprintf("%.2f", volume) + " ul DNA:\n" +
			"    Add " + fmt.Sprintf("%.2f", highRatio*volume) + " ul AMPure beads (ratio=" + fmt.Sprintf("%.2f", highRatio) + ").\n" +
			"    Incubate for 1 minute and then place the tube near a magnet.\n" +
			"    The supernatant contains DNA<=" + fmt.Sprintf("%.2f", highSize) + ".\n" +
			"    Transfer supernatant to a new tube.\n" +
			"    Add " + fmt.Sprintf("%.2f", (lowRatio-highRatio)*volume) + "ul AMPure beads to the previous supernatant.\n" +
			"    Final ratio is " + fmt.Sprintf("%.2f", lowRatio) + ".\n" +
			"    Treat this as normal AMPure clean up:\n" +
			"        Place sample against magnet and incubate until clear.\n" +
			"        Discard the supernatant.\n" +
			"        Add 200 uL of 80% ethanol and incubate for 1 minute. Do not disturb beads.\n" +
			"        Discard supernatant and add 200 uL of 80% ethanol and incubate for 1 minute\n" +
			"        Remove supernatant. While the sample is still against magnet, incubate beads at RT for 3-5 minutes to air dry beads.\n" +
			"        Remove sample from magnet and resupend beads in 25 uL of elution buffer (EB).\n" +
			"        Incubate at room temperature for 5 minutes.\n" +
			"        Place samples against magnet, incubate until solution is clear, and transfer the supernatant to a new tube.\n" +
			"        Quantify by qubit or other methods.\n")
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	// TODO: I used 0 and 10kb, but what should be reasonable limits for the beads and the function?
	const minSizeForBeads float64 = 0
	const maxSizeForBeads float64 = 10000
	var lowSize, highSize, volume float64
	var err error
	var lowText = flag.Arg(0)
	var highText = flag.Arg(1)
	var volumeText = flag.Arg(2)

	// check to make sure input is a valid number and that the numbers are in ranges that make sense
	lowSize, err = strconv.ParseFloat(lowText, 64)
	if err != nil {
		log.Fatalf("Error: expecting a number for the low threshold, but got %s\n", lowText)
	}
	if lowSize <= minSizeForBeads || lowSize >= maxSizeForBeads {
		log.Fatalf("Error: expecting lower threshold to be between %f and %f, but got %f\n", minSizeForBeads, maxSizeForBeads, lowSize)
	}

	highSize, err = strconv.ParseFloat(highText, 64)
	if err != nil {
		log.Fatalf("Error: expecting a number for the high threshold, but got %s\n", highText)
	}
	if highSize <= minSizeForBeads || highSize >= maxSizeForBeads || highSize <= lowSize {
		log.Fatalf("Error: expecting higher threshold to be between %f and %f (and greater than the lower threshold), but got %f\n", minSizeForBeads, maxSizeForBeads, highSize)
	}

	volume, err = strconv.ParseFloat(volumeText, 64)
	if err != nil {
		log.Fatalf("Error: expecting a number for the volume, but got %s\n", volumeText)
	}

	// run program
	calcAmpureRatio(lowSize, highSize, volume)
}
