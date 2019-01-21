package main

import (
	"flag"
	"fmt"
	"log"
	"strconv"
)

func usage() {
	fmt.Print(
		"AMPure Calculator dev - b\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	m := (1.0 - 0.4) / (200 - 800)
	b := 0.4 - m*800
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	var high float64
	var low float64
	var volume float64
	//var volume = flag.Arg(2)
	if lowSize, err := strconv.ParseFloat(flag.Arg(0), 64); err == nil {
		high = m*lowSize + b
	}

	if highSize, err := strconv.ParseFloat(flag.Arg(1), 64); err == nil {
		low = m*highSize + b
	}
	//low := m*highSize + b
	//high := m*lowSize + b
	fmt.Println("To obtain DNA fragments in size range of", flag.Arg(0), "to", flag.Arg(1), "start with AMPure ratio of", fmt.Sprintf("%.2f", low))
	fmt.Println("Transfer supernatant from", fmt.Sprintf("%.2f", low), "to new tube.")
	fmt.Println("Add ", fmt.Sprintf("%.2f", high-low), " extra AMPure to hit target final ratio.")
	fmt.Println("Final ratio needs to be", fmt.Sprintf("%.2f", high))
	fmt.Println()
	fmt.Println("For", flag.Arg(2), "ul DNA:")
	if v, err := strconv.ParseFloat(flag.Arg(2), 64); err == nil {
		volume = v

	}
	fmt.Println("    Add ", fmt.Sprintf("%.2f", low*volume), " ul AMPure beads (ratio =", fmt.Sprintf("%.2f", low), ")")
	fmt.Println("    The supernatant contains DNA <=", flag.Arg(1), ".")
	fmt.Println("    Transfer supernatant new tube.")
	fmt.Println("    Add ", fmt.Sprintf("%.2f", (high-low)*volume), "ul AMPure beads to the previous supernatant.")
	fmt.Println("    Final ratio is", fmt.Sprintf("%.2f", high))
	fmt.Println("    Treat this as normal AMPure cleaning and take DNA off beads")
	fmt.Println("        Place sample against magnet and incubate until clear")
	fmt.Println("        Discard the supernatant")
	fmt.Println("        Add 200 uL of 80% ethanol and incubate for 1 minute. Do not disturb beads")
	fmt.Println("        Discard supernatant and add 200 uL of 80% ethanol and incubate for 1 minute")
	fmt.Println("        Remove supernatant. While samples is still against magnet, incubate beads at RT for 3-5 minutes to air dry beads")
	fmt.Println("        Remove sample from magnet and resupend beads in 25 uL of elution buffer (EB).")
	fmt.Println("        Incubate at room temperature for 5 minutes.")
	fmt.Println("        Place samples against magnet, incubate until solution is clear, and transfer the supernatant to a new tube.")
	fmt.Println("        Quantify by qubit")

}
