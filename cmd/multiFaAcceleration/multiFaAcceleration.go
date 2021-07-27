// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

type Settings struct {
	InFile                string
	ChromName             string
	VelOut                string
	AccelOut              string
	InitialVelOut         string
	SearchSpaceBed        string
	SearchSpaceProportion float64
	WindowSize            int
	Verbose               bool
}

func multiFaAcceleration(s Settings) {
	records := fasta.Read(s.InFile)
	velBed := fileio.EasyCreate(s.VelOut)
	accelBed := fileio.EasyCreate(s.AccelOut)
	initialVelBed := fileio.EasyCreate(s.InitialVelOut)
	var searchSpace []bed.Bed
	var referenceLength, i, j, threshold int
	var bitArray []byte

	//a couple of safety checks
	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration accepts a multiFa file with 4 records, found %v.", len(records))
	}

	if len(records[1].Seq) != len(records[0].Seq) || len(records[2].Seq) != len(records[0].Seq) || len(records[3].Seq) != len(records[0].Seq) {
		log.Fatalf("Error. All records must be of the same sequence length.")
	}

	if s.SearchSpaceBed != "" {
		searchSpace = bed.Read(s.SearchSpaceBed)
		referenceLength = fasta.AlnPosToRefPos(records[0], len(records[0].Seq)-1)
		bitArray = make([]byte, referenceLength)
		for i = range searchSpace {
			if searchSpace[i].Chrom == s.ChromName {
				for j = searchSpace[i].ChromStart; j < searchSpace[i].ChromEnd; j++ {
					bitArray[j] = 1
				}
			}
		}
		threshold = int(s.SearchSpaceProportion * float64(s.WindowSize))
	}

	//set up the system of equations as a blank matrix, with pairwise differences set by default to -1.
	var mat [][]float64 = [][]float64{{1, 1, 0, 0, 0, -1}, {1, 0, 1, 1, 0, -1}, {0, 1, 1, 1, 0, -1}, {1, 0, 1, 0, 1, -1}, {0, 1, 1, 0, 1, -1}, {0, 0, 0, 1, 1, -1}}

	var piS0S1, piS0S2, piS1S2, piS0S3, piS1S3, piS2S3 int
	var referenceCounter int = 0
	var reachedEnd bool = false
	var b1, b2 float64
	var solved [][]float64
	var currCount int

	for alignmentCounter := 0; reachedEnd == false; alignmentCounter++ {
		if s.Verbose && alignmentCounter%1000000 == 0 {
			fmt.Printf("alignmentCounter: %v\n", alignmentCounter)
		}
		if records[0].Seq[alignmentCounter] != dna.Gap && thresholdCheckPasses(s, currCount, threshold, bitArray, referenceCounter) {
			piS0S1, reachedEnd = countWindowDifference(records[0], records[1], alignmentCounter, s.WindowSize)
			piS0S2, _ = countWindowDifference(records[0], records[2], alignmentCounter, s.WindowSize)
			piS1S2, _ = countWindowDifference(records[1], records[2], alignmentCounter, s.WindowSize)
			piS0S3, _ = countWindowDifference(records[0], records[3], alignmentCounter, s.WindowSize)
			piS1S3, _ = countWindowDifference(records[1], records[3], alignmentCounter, s.WindowSize)
			piS2S3, _ = countWindowDifference(records[2], records[3], alignmentCounter, s.WindowSize)
			mat[0][5] = float64(piS0S1)
			mat[1][5] = float64(piS0S2)
			mat[2][5] = float64(piS1S2)
			mat[3][5] = float64(piS0S3)
			mat[4][5] = float64(piS1S3)
			mat[5][5] = float64(piS2S3)
			solved = numbers.Rref(mat)
			b1 = solved[0][5]
			b2 = solved[1][5]
			if !reachedEnd {
				bed.WriteBed(velBed, bed.Bed{Chrom: s.ChromName, ChromStart: referenceCounter, ChromEnd: referenceCounter + s.WindowSize, Name: fmt.Sprintf("%e", b1), FieldsInitialized: 4})
				bed.WriteBed(accelBed, bed.Bed{Chrom: s.ChromName, ChromStart: referenceCounter, ChromEnd: referenceCounter + s.WindowSize, Name: fmt.Sprintf("%e", b2-b1), FieldsInitialized: 4})
				bed.WriteBed(initialVelBed, bed.Bed{Chrom: s.ChromName, ChromStart: referenceCounter, ChromEnd: referenceCounter + s.WindowSize, Name: fmt.Sprintf("%e", b2), FieldsInitialized: 4})
				referenceCounter++
			}
		}
	}
	velBed.Close()
	accelBed.Close()
	initialVelBed.Close()
}

//bitArray is on reference coordinates, not alignment coordinates, so the window is simply equal to windowSize.
func thresholdCheckPasses(s Settings, currCount int, threshold int, bitArray []byte, referenceCounter int) bool {
	if s.SearchSpaceBed == "" { //no search space file, no need to look further
		return true
	}
	if referenceCounter == 0 {
		currCount = 0
		for i := 0; i < s.WindowSize; i++ {
			currCount += int(bitArray[i])
		}
	} else {
		if bitArray[referenceCounter-1] < 1 {
			currCount--
		}
		if bitArray[referenceCounter+s.WindowSize] > 0 {
			currCount++
		}
	}

	return currCount >= threshold
}

func countWindowDifference(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) (int, bool) {
	diff := 0
	baseCount := 0
	var seq1Indel bool = false
	var seq2Indel bool = false
	var reachedEnd bool = false
	var i int = 0

	for i = start; baseCount < windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
				baseCount++
			}
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			baseCount++
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			baseCount++
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}

	if baseCount != windowSize {
		reachedEnd = true
	}
	return diff, reachedEnd
}

func usage() {
	fmt.Print(
		"multiFaAcceleration - Performs velocity and acceleration on a four way multiple alignment in multiFa format." +
			"A four way multiple alignment must contain four species (index 0 to 3) in the topology that aln[0] is the most derived and species 1 to 3 are successive outgroups.\n" +
			"Three bed files are returned. The first produces the velocity score, the second returns the acceleration score, and the third returns the initial velocity score for each window of the genome for aln[0].\n" +
			"Usage:\n" +
			" multiFaAcceleration chromName in.fa velocity.bed acceleration.bed initialVelocity.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var searchSpaceBed *string = flag.String("searchSpaceBed", "", "Limits the generation of data to windows contained within these regions.")
	var searchSpaceProportion *float64 = flag.Float64("searchSpaceProportion", 0.5, "Proportion of window that must overlap search space in order to be evaluated.")
	var windowSize *int = flag.Int("windowSize", 500, "Set the size of the sliding window.")
	var verbose *bool = flag.Bool("verbose", false, "Enables debug prints.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	chromName := flag.Arg(0)
	inFile := flag.Arg(1)
	velOut := flag.Arg(2)
	accelOut := flag.Arg(3)
	initialVOut := flag.Arg(4)

	s := Settings{
		InFile:                inFile,
		ChromName:             chromName,
		VelOut:                velOut,
		AccelOut:              accelOut,
		InitialVelOut:         initialVOut,
		SearchSpaceBed:        *searchSpaceBed,
		SearchSpaceProportion: *searchSpaceProportion,
		WindowSize:            *windowSize,
		Verbose:               *verbose,
	}

	multiFaAcceleration(s)
}
