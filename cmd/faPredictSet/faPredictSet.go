// Command Group: Deep Learning

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

type Settings struct {
	Mode               string
	InFile             string
	OutFile            string
	WindowSize         int
	Stride             int
	WithRevComp        bool
	SatMutagenesisBeds string
	DivergentSitesFile string
	MidpointBed        bool
}

func PredictSetToBed(s Settings) {
	var err error
	var line string
	var midpoint int
	var words, nameFields []string
	var currBed bed.Bed
	var doneReading bool
	in := fileio.EasyOpen(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	for line, doneReading = fileio.EasyNextRealLine(in); !doneReading; line, doneReading = fileio.EasyNextRealLine(in) {
		words = strings.Split(line, "\t")
		if words[1] != "name" { //this skips the header line
			nameFields = strings.Split(words[1], ".")
			currBed = bed.Bed{Chrom: nameFields[0], ChromStart: common.StringToInt(nameFields[1]), ChromEnd: common.StringToInt(nameFields[2]), Name: words[2], FieldsInitialized: 4}
			if s.MidpointBed {
				midpoint = (currBed.ChromStart + currBed.ChromEnd) / 2
				currBed.ChromStart = midpoint
				currBed.ChromEnd = midpoint + 1
			}
			bed.WriteToFileHandle(out, currBed)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func PredictSetToEvolvabilityVector(s Settings) {
	var err error
	var currPos int
	var line, currMut, currOutFileName string
	var currValue, wildTypeValue float64
	var words, nameFields []string
	var doneReading, inMap bool
	var firstTime bool = true
	var divSites []vcf.Vcf
	var divSitesIntervals []interval.Interval = make([]interval.Interval, 0)
	var divSitesTree map[string]*interval.IntervalNode
	var overlappingDivIntervals []interval.Interval
	var currVcf vcf.Vcf
	var currRegion bed.Bed
	var currDivSites map[int]string //map position of divergent site to the mutation identity

	var out *fileio.EasyWriter

	in := fileio.EasyOpen(s.InFile)

	if s.DivergentSitesFile != "" {
		divSites, _ = vcf.Read(s.DivergentSitesFile)
		for i := range divSites {
			divSitesIntervals = append(divSitesIntervals, divSites[i])
		}
		divSitesTree = interval.BuildTree(divSitesIntervals)
	}

	for line, doneReading = fileio.EasyNextRealLine(in); !doneReading; line, doneReading = fileio.EasyNextRealLine(in) {
		words = strings.Split(line, "\t") //parse the first line of the file
		nameFields = strings.Split(words[1], ".")
		currValue = common.StringToFloat64(words[2])
		currPos = common.StringToInt(nameFields[3])
		currMut = nameFields[4]

		if firstTime && currMut != "WT" {
			log.Fatalf("Error. First line of predict file must specify the WT value for a region of interest.")
		}

		if currMut == "WT" {
			//close the old outFile
			if out != nil {
				err = out.Close()
				exception.PanicOnErr(err)
			}

			//make a new outfile for the next region
			currOutFileName = fmt.Sprintf("%s.%s.%s.evolvabilityVector.txt", nameFields[0], nameFields[1], nameFields[2])
			out = fileio.EasyCreate(currOutFileName)
			_, err = fmt.Fprintf(out, "Region\tPositionInRegion\tValue\tCategory\n")
			wildTypeValue = currValue

			currRegion = bed.Bed{Chrom: nameFields[0], ChromStart: common.StringToInt(nameFields[1]), ChromEnd: common.StringToInt(nameFields[2]), Name: fmt.Sprintf("%s.%s.%s.", nameFields[0], nameFields[1], nameFields[2])}

			//find overlapping divergentSites
			if s.DivergentSitesFile != "" {
				overlappingDivIntervals = interval.Query(divSitesTree, currRegion, "any")
				currDivSites = make(map[int]string)
				for i := range overlappingDivIntervals {
					currVcf = overlappingDivIntervals[i].(vcf.Vcf)
					currDivSites[currVcf.Pos-currRegion.ChromStart-1] = currVcf.Alt[0] //find relative position in the sequence. Go from 1-based vcf to 0-based bed coords.
				}
			}
		} else if s.DivergentSitesFile != "" {
			if _, inMap = currDivSites[currPos]; inMap { // if the current site is a divergent site.
				if currMut == currDivSites[currPos] { //if this is true, we have found the divergent mutation
					_, err = fmt.Fprintf(out, "%s\t%d\t%f\tDivSite\n", currRegion.Name, currPos, currValue-wildTypeValue)
				} else {
					_, err = fmt.Fprintf(out, "%s\t%d\t%f\tAllOthers\n", currRegion.Name, currPos, currValue-wildTypeValue)
				}
			} else {
				_, err = fmt.Fprintf(out, "%s\t%d\t%f\tAllOthers\n", currRegion.Name, currPos, currValue-wildTypeValue)
			}
		} else {
			_, err = fmt.Fprintf(out, "%s\t%d\t%f\tAllSites\n", currRegion.Name, currPos, currValue-wildTypeValue)
		}
		firstTime = false
	}
	if out != nil {
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

func PredictSetToHeatmap(s Settings) {
	var currPos, regionLength int
	var currMut, line, currOutFileName string
	var wildTypeValue, currValue float64
	var words, nameFields []string
	var doneReading bool
	var firstTime bool = true
	var currMatrix [][]float64 = make([][]float64, 0)

	in := fileio.EasyOpen(s.InFile)

	for line, doneReading = fileio.EasyNextRealLine(in); !doneReading; line, doneReading = fileio.EasyNextRealLine(in) {
		//parse the first line and set the value of the wild type allele
		words = strings.Split(line, "\t") //parse the first line of the file
		nameFields = strings.Split(words[1], ".")
		currMut = nameFields[4]
		currPos = common.StringToInt(nameFields[3])
		currValue = common.StringToFloat64(words[2])
		regionLength = common.StringToInt(nameFields[2]) - common.StringToInt(nameFields[1])

		if firstTime && currMut != "WT" {
			log.Fatalf("Error. First line of predict file must specify the WT value for a region of interest.")
		}

		if currMut == "WT" {
			wildTypeValue = currValue
			currOutFileName = fmt.Sprintf("%s.%s.%s.heatmap.txt", nameFields[0], nameFields[1], nameFields[2])
			if !firstTime {
				writeMatrixToFile(currMatrix, currOutFileName)
			}
			currMatrix = initializeNewHeatmap(regionLength)
		} else if currMut == "A" {
			currMatrix[currPos][0] = currValue - wildTypeValue
		} else if currMut == "C" {
			currMatrix[currPos][1] = currValue - wildTypeValue
		} else if currMut == "G" {
			currMatrix[currPos][2] = currValue - wildTypeValue
		} else if currMut == "T" {
			currMatrix[currPos][3] = currValue - wildTypeValue
		} else {
			log.Fatalf("Mut not recognized.")
		}

		firstTime = false
	}
	writeMatrixToFile(currMatrix, currOutFileName)
}

func writeMatrixToFile(answer [][]float64, outFileName string) {
	out := fileio.EasyCreate(outFileName)
	var i, j int
	var lineToWrite string
	var err error

	for i = range answer {
		lineToWrite = ""
		for j = range answer[i] {
			lineToWrite = fmt.Sprintf("%s\t%f", lineToWrite, answer[i][j])
		}
		_, err = fmt.Fprintf(out, "%s\n", lineToWrite)
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func initializeNewHeatmap(regionLength int) [][]float64 {
	var i int
	//initialize output matrix as zeros
	var answer = make([][]float64, regionLength) //length of matrix is length of region (typically 400)
	for i = range answer {
		answer[i] = make([]float64, 4) //matrix of depth 4, for the four base identities
	}

	return answer
}

func FaToSatMutagenesisPredictSet(s Settings) {
	var err error
	var i, j int
	var currIndex int
	var currFa fasta.Fasta
	var lineToWrite string
	var currMutant []dna.Base

	chroms := fasta.Read(s.InFile)
	out := fileio.EasyCreate(s.OutFile)
	beds := bed.Read(s.SatMutagenesisBeds)

	for i = range beds {
		currIndex = fasta.GetChromIndex(chroms, beds[i].Chrom)
		currFa = fasta.Extract(chroms[currIndex], beds[i].ChromStart, beds[i].ChromEnd, fmt.Sprintf("%s.%d.%d", beds[i].Chrom, beds[i].ChromStart, beds[i].ChromEnd))
		dna.AllToUpper(currFa.Seq)
		lineToWrite = fmt.Sprintf("%s.0.WT\t%s\n", currFa.Name, dna.BasesToString(currFa.Seq))
		_, err = fmt.Fprintf(out, lineToWrite)
		exception.PanicOnErr(err)
		currMutant = make([]dna.Base, len(currFa.Seq))
		for j = range currFa.Seq {
			if currFa.Seq[j] != dna.A {
				copy(currMutant, currFa.Seq)
				currMutant[j] = dna.A
				lineToWrite = fmt.Sprintf("%s.%d.A\t%s\n", currFa.Name, j, dna.BasesToString(currMutant))
				_, err = fmt.Fprintf(out, lineToWrite)
				exception.PanicOnErr(err)
			}
			if currFa.Seq[j] != dna.C {
				copy(currMutant, currFa.Seq)
				currMutant[j] = dna.C
				lineToWrite = fmt.Sprintf("%s.%d.C\t%s\n", currFa.Name, j, dna.BasesToString(currMutant))
				_, err = fmt.Fprintf(out, lineToWrite)
				exception.PanicOnErr(err)
			}
			if currFa.Seq[j] != dna.G {
				copy(currMutant, currFa.Seq)
				currMutant[j] = dna.G
				lineToWrite = fmt.Sprintf("%s.%d.G\t%s\n", currFa.Name, j, dna.BasesToString(currMutant))
				_, err = fmt.Fprintf(out, lineToWrite)
				exception.PanicOnErr(err)
			}
			if currFa.Seq[j] != dna.T {
				copy(currMutant, currFa.Seq)
				currMutant[j] = dna.T
				lineToWrite = fmt.Sprintf("%s.%d.T\t%s\n", currFa.Name, j, dna.BasesToString(currMutant))
				_, err = fmt.Fprintf(out, lineToWrite)
				exception.PanicOnErr(err)
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func FaToPredictSetComplete(s Settings) {
	var err error
	var i, j int
	records := fasta.Read(s.InFile)
	out := fileio.EasyCreate(s.OutFile)
	var currFa fasta.Fasta
	var lineToWrite string
	var revSeq []dna.Base = make([]dna.Base, s.WindowSize)

	for i = range records {
		for j = 0; j < len(records[i].Seq)-s.WindowSize; j += s.Stride {
			currFa = fasta.Extract(records[i], j, j+s.WindowSize, fmt.Sprintf("%s.%d.%d", records[i].Name, j, j+s.WindowSize))
			dna.AllToUpper(currFa.Seq)
			if s.WithRevComp {
				copy(revSeq, currFa.Seq)
				dna.ReverseComplement(revSeq)
				lineToWrite = fmt.Sprintf("%s\t%s\t%s\n", currFa.Name, dna.BasesToString(currFa.Seq), dna.BasesToString(revSeq))
			} else {
				lineToWrite = fmt.Sprintf("%s\t%s\n", currFa.Name, dna.BasesToString(currFa.Seq))
			}
			_, err = fmt.Fprintf(out, lineToWrite)
			exception.PanicOnErr(err)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"faPredictSet - Make and analyze deep learning prediction TSV files from input fasta format data.\n" +
			"Usage:\n" +
			"faPredictSet ToPredictSet input.fa output.txt\n" +
			"OR\n" +
			"faPredictSet PredictSetToHeatmaps predicted.txt\n" +
			"OR\n" +
			"faPredictSet PredictSetToEvolvabilityVector predicted.txt\n" +
			"OR\n" +
			"faPredictSet PredictSetToBed predicted.txt midpointPredictions.bed" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var windowSize *int = flag.Int("windowSize", 400, "Set the size of the sequence window in ToPredictSet.")
	var stride *int = flag.Int("stride", 1, "Sets the size of the stride between prediction windows in ToPredictSet.")
	var withRevComp *bool = flag.Bool("withRevComp", false, "Include the reverse complement sequence in the output file as an extra column in ToPredictSet.")
	var satMutagenesisBeds *string = flag.String("satMutagenesisBeds", "", "Specify a file to perform saturation mutagenesis on a set of genomic regions in ToPredictSet.")
	var DivergentSitesFile *string = flag.String("DivergentSitesFile", "", "Specify a VCF file of divergent sites to split the evolvability vector into divergent sites and all other mutations.")
	var midpointBed *bool = flag.Bool("midpointBed", false, "In mode: PredictSetToBed, return a bed of length 1 representing the midpoint of the prediction window.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) == 0 {
		flag.Usage()
		log.Fatalf("Error: Received no input arguments.")
	}

	mode := flag.Arg(0)
	var inFile, outFile string
	var s Settings

	if mode == "ToPredictSet" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		inFile = flag.Arg(1)
		outFile = flag.Arg(2)
		s = Settings{
			Mode:               mode,
			InFile:             inFile,
			OutFile:            outFile,
			WindowSize:         *windowSize,
			Stride:             *stride,
			WithRevComp:        *withRevComp,
			SatMutagenesisBeds: *satMutagenesisBeds,
		}
		if s.SatMutagenesisBeds != "" {
			FaToSatMutagenesisPredictSet(s)
		} else {
			FaToPredictSetComplete(s)
		}
	} else if mode == "PredictSetToHeatmaps" {
		expectedNumArgs = 2
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		inFile = flag.Arg(1)
		s = Settings{
			Mode:               mode,
			InFile:             inFile,
			WindowSize:         *windowSize,
			Stride:             *stride,
			WithRevComp:        *withRevComp,
			SatMutagenesisBeds: *satMutagenesisBeds,
		}
		PredictSetToHeatmap(s)
	} else if mode == "PredictSetToEvolvabilityVector" {
		expectedNumArgs = 2
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		inFile = flag.Arg(1)
		s = Settings{
			Mode:               mode,
			InFile:             inFile,
			WindowSize:         *windowSize,
			Stride:             *stride,
			WithRevComp:        *withRevComp,
			SatMutagenesisBeds: *satMutagenesisBeds,
			DivergentSitesFile: *DivergentSitesFile,
		}
		PredictSetToEvolvabilityVector(s)
	} else if mode == "PredictSetToBed" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		inFile = flag.Arg(1)
		outFile = flag.Arg(2)
		s = Settings{
			Mode:        mode,
			InFile:      inFile,
			OutFile:     outFile,
			MidpointBed: *midpointBed,
		}
		PredictSetToBed(s)
	} else {
		log.Fatalf("Error: Unrecognized mode. See usage.")
	}
}
