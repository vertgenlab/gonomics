package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func compareBedOnName(eStrBed string, geneInfoBed string, outFile string) {
	out := fileio.EasyCreate(outFile)
	eStr := bed.Read(eStrBed)
	geneInfo := bed.Read(geneInfoBed)

	geneInfoMap := make(map[string][]bed.Bed)
	var i int

	for i = range geneInfo {
		geneInfoMap[geneInfo[i].Name] = append(geneInfoMap[geneInfo[i].Name], geneInfo[i])
		}

	var eStrGene []bed.Bed
	var found bool
	var distanceTss int
	var lowestDistanceTss int
	var err error
	var j int

	for i := 0; i < len(eStr); i++ {
		lowestDistanceTss = 0
		eStrGene, found = geneInfoMap[eStr[i].Name]
		if found != true {
			log.Fatalf("Did not find %s", eStr[i].Name)
		}
		for j = range eStrGene {
			log.Println(eStrGene[j])
			log.Printf("%v", j)
			distanceTss, err = bed.CompareDistance(eStrGene[j], eStr[i])
			log.Printf("%v", distanceTss)
			if j == 0 {
				lowestDistanceTss = distanceTss
				log.Printf("%v", lowestDistanceTss)
			} else {
				if distanceTss < lowestDistanceTss {
					lowestDistanceTss = distanceTss
				}
			}
		}

		if err != nil {
			log.Fatalf("Unable to compare distance, error message: %s", err)
		}
		eStr[i].Score = lowestDistanceTss
		bed.WriteToFileHandle(out, eStr[i])
	}
		err = out.Close()
		if err != nil{
			log.Fatalf("File unable to close properly, error is: %s", err)
		}
}

func usage() {
	fmt.Print(
		"compareBedDistanceBasedOnName - Compares the Name fields between two bed files, and if they match,\n" +
			"the minimum distance between those coordinates is calculated.\n" +
			"The second bed must include at least every Name included in the first bed. This second bed can \n" +
			"have more names than the ones that are being compared.\n" +
			"For example the second bed might be the coordinates for every gene, genome wide\n" +
			"while the first bed includes a list of bed coordinates you are interested in studying.\n" +
			"The output file will be the first input bed file where the minimum distance is reported in the Score column.\n" +
			"No headers allowed in input files.\n" +
			"Must have a zero in the score field for the first bed, for gonomics to be able to write the distance to this field.\n" +
			"This program was originally written for comparing eSTR coordinates to gene coordinates.\n" +
			"Usage:\n" +
			"eStr.bed geneInfo.bed out.bed\n")
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

	eStrBed := flag.Arg(0)
	geneInfoBed := flag.Arg(1)
	outFile := flag.Arg(2)

	compareBedOnName(eStrBed, geneInfoBed, outFile)

}