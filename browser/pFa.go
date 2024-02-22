package browser

import (
	"fmt"

	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/fasta/pfasta"
)

// PFaVisualizer produces command line visualizations of pFasta format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes noMask and lineLength formatting options as bools.
// If 0 sig figs, returns full probability
func PFaVisualizer(infile string, outfile string, start int, end int, sigFigs int, lineLength int, seqName string) {
	//add 'end' as input in cmd (look at multfavisualiser)
	if !(start < end) {
		log.Fatalf("Error: Invalid arguments, start must be lower than end")
	}

	records := Read(infile) // for now, assuming only 1 seq in pfasta, will rewrite for multiple
	out := fileio.EasyCreate(outfile)
	lineA := make([]float32, lineLength)
	lineC := make([]float32, lineLength)
	lineG := make([]float32, lineLength)
	lineT := make([]float32, lineLength)
	setOfLinesIdx := 0

	fmt.Fprintf(out, "%s", records[0].Name)
	for setOfLinesIdx = 0; setOfLinesIdx < (end-start)/sigFigs; setOfLinesIdx++ {
		printOneSetLines(lineLength, setOfLinesIdx, lineLength, lineA, lineC, lineG, lineT, start, records, out, sigFigs)
	}
	printOneSetLines(lineLength, setOfLinesIdx, (end-start)%sigFigs, lineA, lineC, lineG, lineT, start, records, out, sigFigs)
}

// printOneSetLines prints from init_pos =(setOfLinesIdx*lineLength + start) in pFasta to (init_pos + numIters )
func printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, records []PFasta, out *fileio.EasyWriter, sigFigs int, name string) {
	// add start position (look at mulfavisualiser), empty spacer lines
	recordIdx := setOfLinesIdx*lineLength + start
	lineIdx := 0
	for lineIdx = 0; lineIdx < numIters; lineIdx++ {
		fmt.Printf("sigfigs:%v", sigFigs)
		lineA[lineIdx], lineC[lineIdx], lineG[lineIdx], lineT[lineIdx] = getBaseProbsAtPos(records[0].Seq[recordIdx], sigFigs)
		recordIdx++
	}
	fmt.Fprintf(out, "Position: %d\n", recordIdx)
	fmt.Fprintf(out, "|%-*s | A | %v\n", lineA[0:lineIdx])
	fmt.Fprintf(out, "C | %v\n", lineC[0:lineIdx])
	fmt.Fprintf(out, "G | %v\n", lineG[0:lineIdx])
	fmt.Fprintf(out, "T | %v\n\n", lineT[0:lineIdx])
}

// getBaseProbsAtPos returns the four probabilities rounded to sigFigs for a specified base
func getBaseProbsAtPos(base pDna.Float32Base, sigFigs int) (float32, float32, float32, float32) {
	// return roundSigFigs(base.A, sigFigs), roundSigFigs(base.C, sigFigs), roundSigFigs(base.G, sigFigs), roundSigFigs(base.T, sigFigs)
	return numbers.RoundSigFigs(float64(base.A), sigFigs), numbers.RoundSigFigs(float64(base.C), sigFigs), numbers.RoundSigFigs(float64(base.G), sigFigs), numbers.RoundSigFigs(float64(base.T), sigFigs)
}

// calculateLongestName is a helper function of MultiFaVisualizer that returns the length of the longest name in a slice of fasta.Fasta structs.
func calculateLongestName(f []pFasta.PFasta) int {
	var ans int = 0
	var temp int
	for i := range f {
		temp = utf8.RuneCountInString(f[i].Name)
		if temp > ans {
			ans = temp
		}
	}
	return ans
}

