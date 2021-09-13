package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func compareBedDistanceBasedOnName (inputBed string, genomeBed string, outBed string) {
	out := fileio.EasyCreate(outBed)
	input := bed.Read(inputBed)
	genome := bed.Read(genomeBed)

	genomeMap := make(map[string][]bed.Bed)
	var i int

	for i = range genome {
		genomeMap[genome[i].Name] = append(genomeMap[genome[i].Name], genome[i])
		}

	var inputNameField []bed.Bed
	var found bool
	var distance int
	var lowestDistance int
	var err error
	var j int

	for i := 0; i < len(input); i++ {
		lowestDistance = 0
		inputNameField, found = genomeMap[input[i].Name]
		if found != true {
			log.Fatalf("Did not find %s", input[i].Name)
		}
		for j = range inputNameField {
			distance, err = bed.CompareDistance(inputNameField[j], input[i])
			if j == 0 {
				lowestDistance = distance
			} else {
				if distance < lowestDistance {
					lowestDistance = distance
				}
			}
		}

		if err != nil {
			log.Fatalf("Unable to compare distance, error message: %s", err)
		}
		input[i].Score = lowestDistance
		if input[i].FieldsInitialized < 5 {
			input[i].FieldsInitialized = 5
		}
		bed.WriteToFileHandle(out, input[i])
	}
		err = out.Close()
		if err != nil{
			log.Fatalf("File unable to close properly, error is: %s", err)
		}
}

func usage() {
	fmt.Print(
		"compareBedDistanceBasedOnName - Compares the Name fields between \n" +
			"two bed files, and if they match, the minimum distance between those\n" +
			"coordinates is calculated.\n" +
			"The first bed is iterated over line by line and the second bed is\n" +
			"searched for the entries that match the Name field.\n" +
			"Every entry in the first bed is expected to have at least one corresponding\n" +
			"bed entry in the second bed, where more than one matching bed entry is allowed\n" +
			"If there is more than one bed entry in the second bed that matches an entry in" +
			"the first bed, the distance between the coordinates will be calculated pairwise\n" +
			"and the lowest value will be reported (first bed paired with each matching entry from the second bed).\n" +
			"Thus, the second bed must include at least every Name included in the first bed.\n" +
			"For example the second bed might be the coordinates for every gene, genome wide (genomeBed)\n" +
			"while the first bed includes a list of bed coordinates you are interested in studying (inputBed).\n" +
			"The output file will be the inputBed file with the minimum distance reported in the Score column.\n" +
			"No headers allowed in input files.\n" +
			"This program was originally written for comparing eSTR coordinates to gene coordinates.\n" +
			"Usage:\n" +
			"input.bed genome.bed out.bed\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 3
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inputBed := flag.Arg(0)
	genomeBed := flag.Arg(1)
	outBed := flag.Arg(2)

	compareBedDistanceBasedOnName(inputBed, genomeBed, outBed)

}