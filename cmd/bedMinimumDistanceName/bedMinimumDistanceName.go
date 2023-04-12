package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

func bedMinimumDistanceName(inputBed string, genomeBed string, outBed string) {
	out := fileio.EasyCreate(outBed)
	input := bed.Read(inputBed)
	genome := bed.Read(genomeBed)

	genomeMap := make(map[string][]bed.Bed)
	var i int
	var ok bool

	//Create a map of all the "Names" in the genomeBed for easy reference.
	//Assumes all genomeBed entries have unique "Name" fields.
	for i = range genome {
		if _, ok = genomeMap[genome[i].Name]; ok {
			log.Fatalf("The following entry (Name field) is found twice in the bed: %s", genome[i].Name)
		}
		genomeMap[genome[i].Name] = append(genomeMap[genome[i].Name], genome[i])
	}

	var genomeBedMatchingNameField []bed.Bed
	var found bool
	var lowestDistance int
	var upstreamDownstreamIndicator bed.Strand
	var err error

	for i = range input { //range over the inputBed
		genomeBedMatchingNameField, found = genomeMap[input[i].Name] //grab the genomeBed entry that matches on the "Name" field.
		if found != true {
			log.Fatalf("Did not find genomeBed match for: %s", input[i].Name)
		}
		if len(genomeBedMatchingNameField) > 1 { //check that there is only one genomeBed entry that matches the current inputBed[i].Name
			log.Fatalf("There is more than one genomeBed.Name entry that matches %s", input[i])
		}
		lowestDistance, err = bed.MinimumDistance(genomeBedMatchingNameField[0], input[i])
		if err != nil {
			log.Fatalf("Unable to compare distance, error message: %s", err)
		}
		upstreamDownstreamIndicator = determineUpstreamDownstream(input[i], genomeBedMatchingNameField[0])

		input[i].Score = lowestDistance
		input[i].Strand = upstreamDownstreamIndicator
		if input[i].FieldsInitialized < 6 {
			input[i].FieldsInitialized = 6
		}
		bed.WriteToFileHandle(out, input[i])
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func determineUpstreamDownstream(inputBed bed.Bed, genomeBed bed.Bed) bed.Strand {
	var outputStrand bed.Strand
	if genomeBed.Strand == '+' {
		if (inputBed.ChromStart < genomeBed.ChromStart) || (inputBed.ChromStart == genomeBed.ChromStart) {
			outputStrand = '-'
		}
		if inputBed.ChromStart > genomeBed.ChromStart {
			outputStrand = '+'
		}
	} else if genomeBed.Strand == '-' {
		if (inputBed.ChromEnd > genomeBed.ChromEnd) || (inputBed.ChromEnd == genomeBed.ChromEnd) {
			outputStrand = '-'
		}
		if inputBed.ChromEnd < genomeBed.ChromEnd {
			outputStrand = '+'
		}
	} else {
		log.Fatalf("problem with genomeBed strand: %s", genomeBed)
	}
	return outputStrand
}

//TODO: update usage statement. Explain that genomeBed needs strand.
func usage() {
	fmt.Print(
		"bedMinimumDistanceName - For all entries in bed A (inputBed), look at the Name field and\n" +
			"find the matching Name field in bed B (genomeBed). Calculates the minimum distance between\n" +
			"the two bed entries and determines if the inputBed entry is upstream or downstream of the matching\n" +
			"genomeBed entry. Upstream is denoted as '-' and downstream is denoted as '+'. \n" +
			"The outputBed is the inputBed, with minimum distance reported in the 'Score' field and the \n" +
			"upstream/downstream labeling reported in the 'Strand' field. \n" +
			"Both 'Score' and 'Strand' fields will be overwritten from the inputBed if these \n" +
			"fields are not blank.\n" +
			"inputBed expectations: \n" +
			"\t - must have Chr, ChromStart, ChromEnd, and Name field\n" +
			"\t - optional to have Score and Strand but these will be rewritten\n" +
			"\t - all subsequent fields will be retained in the outputBed\n" +
			"\t - Name fields do not need to be unique" +
			"genomeBed expectations:\n" +
			"\t - must have Chr, ChromStart, ChromEnd, Name, and Strand\n" +
			"\t - all Name field entries must be unique.\n" +
			"\t - must include at least every Name from the inputBed, but it can have more.\n" +
			"No headers allowed in either input files.\n" +
			"This program was originally written for comparing eSTR coordinates to gene TSS coordinates.\n" +
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

	bedMinimumDistanceName(inputBed, genomeBed, outBed)

}
