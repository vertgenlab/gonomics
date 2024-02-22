package pFasta

import (
	"log"
	"fmt"
	"math"
	"math/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/roundSigFigs"
)

// checks if input pFasta has a sequence with chrom as name and returns its index
func checkIfChromInPfasta(input []PFasta, chrom string) int {
	chromInInput := false
	var answer int

	for inputIdx, inputpFa := range input {
		if inputpFa.Name == chrom {
			chromInInput = true
			answer = inputIdx
		}
	}

	if !chromInInput {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}

	return answer
}

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input []PFasta, start int, end int, outputName string, chrom string, takeCoords bool) PFasta {

	chromIdx := checkIfChromInPfasta(input, chrom)

	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input[chromIdx].Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var outName string
	if takeCoords {
		outName = fmt.Sprintf("%s:%v-%v", chrom, start, end)
	} else if len(outputName) > 0 {
		outName = outputName
	} else {
		outName = chrom
	}

	var answer = PFasta{Name: outName, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input[chromIdx].Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a pFa that has a list of subsequences of the input pFa
// defined by the regions in the bed region
// takeCoords specifies if name fields in output should be original names in region or identified by ChromStart and ChromEnd
func ExtractBed(input []PFasta, region []bed.Bed, takeCoords bool) []PFasta {
	answer := make([]PFasta, 0)
	for _, reg := range region {
		answer = append(answer, Extract(input, reg.ChromStart, reg.ChromEnd, "", reg.Chrom, takeCoords))
	}
	return answer
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution
func Sample(input []PFasta, chrom string) fasta.Fasta {
	chromIdx := checkIfChromInPfasta(input, chrom)

	var answer = fasta.Fasta{Name: input[chromIdx].Name, Seq: make([]dna.Base, len(input[chromIdx].Seq))}
	var currRand float32
	for inputIdx := range input[chromIdx].Seq {
		currRand = rand.Float32()
		if currRand < input[chromIdx].Seq[inputIdx].A {
			answer.Seq[inputIdx] = dna.A
		} else if currRand < (input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.C
		} else if currRand < (input[chromIdx].Seq[inputIdx].G + input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.G
		} else {
			answer.Seq[inputIdx] = dna.T
		}
	}

	return answer
}

// Browse produces command line visualizations of pFasta format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes noMask and lineLength formatting options as bools.
// If 0 sig figs, returns full probability
// func browsePfasta(infile string, outfile string, start int, end int, sigFigs int, lineLength int, seqName string) {
// 	//add 'end' as input in cmd (look at multfavisualiser)
// 	if !(start < end) {
// 		log.Fatalf("Error: Invalid arguments, start must be lower than end")
// 	}

// 	records := Read(infile) // for now, assuming only 1 seq in pfasta, will rewrite for multiple
// 	out := fileio.EasyCreate(outfile)
// 	lineA := make([]float32, lineLength)
// 	lineC := make([]float32, lineLength)
// 	lineG := make([]float32, lineLength)
// 	lineT := make([]float32, lineLength)
// 	setOfLinesIdx := 0

// 	for setOfLinesIdx = 0; setOfLinesIdx < (end-start)/sigFigs; setOfLinesIdx++ {
// 		printOneSetLines(lineLength, setOfLinesIdx, lineLength, lineA, lineC, lineG, lineT, start, records, out, sigFigs)
// 	}
// 	printOneSetLines(lineLength, setOfLinesIdx, (end-start)%sigFigs, lineA, lineC, lineG, lineT, start, records, out, sigFigs)
// }

// printOneSetLines prints from init_pos =(setOfLinesIdx*lineLength + start) in pFasta to (init_pos + numIters )
func printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, records []PFasta, out *fileio.EasyWriter, sigFigs int) {
	// add start position (look at mulfavisualiser), empty spacer lines
	recordIdx := setOfLinesIdx*lineLength + start
	lineIdx := 0
	for lineIdx = 0; lineIdx < numIters; lineIdx++ {
		fmt.Printf("sigfigs:%v", sigFigs)
		lineA[lineIdx], lineC[lineIdx], lineG[lineIdx], lineT[lineIdx] = getBaseProbsAtPos(records[0].Seq[recordIdx], sigFigs)
		recordIdx++
	}
	fmt.Fprintf(out, "%v", lineA[0:lineIdx])
	fmt.Fprintf(out, "%v", lineC[0:lineIdx])
	fmt.Fprintf(out, "%v", lineG[0:lineIdx])
	fmt.Fprintf(out, "%v", lineT[0:lineIdx])
}

// getBaseProbsAtPos returns the four probabilities rounded to sigFigs for a specified base
func getBaseProbsAtPos(base pDna.Float32Base, sigFigs int) (float32, float32, float32, float32) {
	// return roundToSigFigs(base.A, sigFigs), roundToSigFigs(base.C, sigFigs), roundToSigFigs(base.G, sigFigs), roundToSigFigs(base.T, sigFigs)
	return roundToSigFigs(float64(base.A), sigFigs), roundToSigFigs(float64(base.C), sigFigs), roundToSigFigs(float64(base.G), sigFigs), roundToSigFigs(float64(base.T), sigFigs)
}


