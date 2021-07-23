package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func multiFaAcceleration(inFile string, chromName string, velOut string, accelOut string, constraintOut string, windowSize int) {
	records := fasta.Read(inFile)

	//a couple of safety checks
	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration accepts a multiFa file with 4 records, found %v.", len(records))
	}

	if len(records[1].Seq) != len(records[0].Seq) || len(records[2].Seq) != len(records[0].Seq) || len(records[3].Seq) != len(records[0].Seq) {
		log.Fatalf("Error. All records must be of the same sequence length.")
	}


	velWig := wig.Wig{StepType: "fixedStep", Chrom: chromName, Start: 1, Step: 1, Values: make([]float64, len(records[0].Seq))}
	accelWig := wig.Wig{StepType: "fixedStep", Chrom: chromName, Start: 1, Step: 1, Values: make([]float64, len(records[0].Seq))}
	constraintWig := wig.Wig{StepType: "fixedStep", Chrom: chromName, Start: 1, Step: 1, Values: make([]float64, len(records[0].Seq))}

	//set up the system of equations as a blank matrix, with pairwise differences set by default to -1.
	var mat [][]float64 = [][]float64{{1, 1, 0, 0, 0, -1},{1, 0, 1, 1, 0, -1},{0, 1, 1, 1, 0, -1},{1, 0, 1, 0, 1, -1},{0, 1, 1, 0, 1, -1},{0, 0, 0, 1, 1, -1}}

	var piS0S1, piS0S2, piS1S2, piS0S3, piS1S3, piS2S3 int
	var referenceCounter int = 0
	var reachedEnd bool = false
	var b1, b2, b3, b4, b5 float64
	var solved [][]float64

	for alignmentCounter := 0; reachedEnd == false; alignmentCounter++ {
		if records[0].Seq[alignmentCounter] != dna.Gap {
			piS0S1, reachedEnd = countWindowDifference(records[0], records[1], alignmentCounter, windowSize)
			piS0S2, _ = countWindowDifference(records[0], records[2], alignmentCounter, windowSize)
			piS1S2, _ = countWindowDifference(records[1], records[2], alignmentCounter, windowSize)
			piS0S3, _ = countWindowDifference(records[0], records[3], alignmentCounter, windowSize)
			piS1S3, _ = countWindowDifference(records[1], records[3], alignmentCounter, windowSize)
			piS2S3, _ = countWindowDifference(records[2], records[3], alignmentCounter, windowSize)
			mat[0][5] = float64(piS0S1)
			mat[1][5] = float64(piS0S2)
			mat[2][5] = float64(piS1S2)
			mat[3][5] = float64(piS0S3)
			mat[4][5] = float64(piS1S3)
			mat[5][5] = float64(piS2S3)
			solved = numbers.Rref(mat)
			b1 = solved[0][5]
			b2 = solved[1][5]
			b3 = solved[2][5]
			b4 = solved[3][5]
			b5 = solved[4][5]
			if !reachedEnd {
				velWig.Values[referenceCounter+(windowSize / 2)] = b1
				accelWig.Values[referenceCounter+(windowSize / 2)] = b2 - b1
				constraintWig.Values[referenceCounter+(windowSize / 2)] = 1.0 / (b2 + b3 + b4 + b5)
				referenceCounter++
			}
		}
	}

	wig.Write(velOut, []wig.Wig{velWig})
	wig.Write(accelOut, []wig.Wig{accelWig})
	wig.Write(constraintOut, []wig.Wig{constraintWig})
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
			"Two wig files are returned. One produces the velocity score and the second returns the acceleration score for each window of the genome for aln[0].\n" +
			"Usage:\n" +
			" multiFaAcceleration chromName in.fa velocity.wig acceleration.wig constraint.wig\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var windowSize *int = flag.Int("windowSize", 500, "Set the size of the sliding window.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chromName := flag.Arg(1)
	velOut := flag.Arg(2)
	accelOut := flag.Arg(3)
	constraintOut := flag.Arg(4)

	multiFaAcceleration(inFile, chromName, velOut, accelOut, constraintOut, *windowSize)
}