package browser

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"log"
	"strings"
	"unicode/utf8"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
)

// PFaVisualizer produces command line visualizations of pFasta format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes noMask and lineLength formatting options as bools.
// If 0 sig figs, returns full probability
func PFaVisualizer(infile string, outfile string, start int, end int, endOfAlignment bool, sigFigs int, lineLength int, seqName string) {
	//add 'end' as input in cmd (look at multfavisualiser)
	if !endOfAlignment && !(start < end) {
		log.Fatalf("Error: Invalid arguments, start must be lower than end")
	}

	records := pFasta.Read(infile) // for now, assuming only 1 seq in pfasta, will rewrite for multiple
	out := fileio.EasyCreate(outfile)
	lineA := make([]float32, lineLength)
	lineC := make([]float32, lineLength)
	lineG := make([]float32, lineLength)
	lineT := make([]float32, lineLength)
	setOfLinesIdx := 0
	if !endOfAlignment {
		fmt.Fprintf(out, "Start: %d. End: %d. SigFigs: %d.", start, end, sigFigs)
	} else {
		fmt.Fprintf(out, "Start: %d. End: to end. SigFigs: %d.", start, sigFigs)
	}

	if seqName == "" {
		longestName := calculateLongestNamePFa(records)
		// need to reconfigure this so that it can auto stop when it's "going off the chrom/approaching end"
		// WILL NOT WORK IF ALL RECORDS DO NOT HAVE SAME LENGTH
		for setOfLinesIdx = 0; setOfLinesIdx < (end-start)/lineLength; setOfLinesIdx++ {
			for seqIdx, seq := range records {
				fmt.Fprintf(out, "\n")
				printOneSetLines(lineLength, setOfLinesIdx, lineLength, lineA, lineC, lineG, lineT, start, seq, out, sigFigs, records[seqIdx].Name, longestName)
			}
		}

		for seqIdx, seq := range records {
			fmt.Fprintf(out, "\n")
			printOneSetLines(lineLength, setOfLinesIdx, (end-start)%lineLength, lineA, lineC, lineG, lineT, start, seq, out, sigFigs, records[seqIdx].Name, longestName)
		}
	} else {
		for desiredSeqIdx, desiredSeq := range records {
			if desiredSeq.Name == seqName {
				longestName := len(desiredSeq.Name)
				for setOfLinesIdx = 0; setOfLinesIdx < (end-start)/lineLength; setOfLinesIdx++ {
					fmt.Fprintf(out, "\n")
					printOneSetLines(lineLength, setOfLinesIdx, lineLength, lineA, lineC, lineG, lineT, start, desiredSeq, out, sigFigs, records[desiredSeqIdx].Name, longestName)
				}
				fmt.Fprintf(out, "\n")
				printOneSetLines(lineLength, setOfLinesIdx, (end-start)%lineLength, lineA, lineC, lineG, lineT, start, desiredSeq, out, sigFigs, records[desiredSeqIdx].Name, longestName)
				break
			}
		}
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

// printOneSetLines prints from init_pos =(setOfLinesIdx*lineLength + start) in pFasta to (init_pos + numIters )
func printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, record pFasta.PFasta, out *fileio.EasyWriter, sigFigs int, name string, longestName int) {
	recordIdx := setOfLinesIdx*lineLength + start
	lineIdx := 0
	fmt.Fprintf(out, "Position: %d\n", recordIdx)
	for lineIdx = 0; lineIdx < numIters; lineIdx++ {
		if sigFigs == 0 {
			lineA[lineIdx], lineC[lineIdx], lineG[lineIdx], lineT[lineIdx] = record.Seq[recordIdx].A, record.Seq[recordIdx].C, record.Seq[recordIdx].G, record.Seq[recordIdx].T
		} else {
			lineA[lineIdx], lineC[lineIdx], lineG[lineIdx], lineT[lineIdx] = getBaseProbsAtPos(record.Seq[recordIdx], sigFigs)
		}
		recordIdx++
	}
	fmt.Fprintf(out, ">%-*s\t|\tA\t|\t%v\n", longestName, record.Name, arrayToString(lineA[0:lineIdx], "\t"))
	fmt.Fprintf(out, "|%-*s\t|\tC\t|\t%v\n", longestName, "", arrayToString(lineC[0:lineIdx], "\t"))
	fmt.Fprintf(out, "|%-*s\t|\tG\t|\t%v\n", longestName, "", arrayToString(lineG[0:lineIdx], "\t"))
	fmt.Fprintf(out, "|%-*s\t|\tT\t|\t%v\n", longestName, "", arrayToString(lineT[0:lineIdx], "\t"))
}

// getBaseProbsAtPos returns the four probabilities rounded to sigFigs for a specified base
func getBaseProbsAtPos(base pDna.Float32Base, sigFigs int) (float32, float32, float32, float32) {
	return numbers.RoundSigFigs(float64(base.A), sigFigs), numbers.RoundSigFigs(float64(base.C), sigFigs), numbers.RoundSigFigs(float64(base.G), sigFigs), numbers.RoundSigFigs(float64(base.T), sigFigs)
}

// calculateLongestName is a helper function of MultiFaVisualizer that returns the length of the longest name in a slice of fasta.Fasta structs.
func calculateLongestNamePFa(f []pFasta.PFasta) int {
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

// https://stackoverflow.com/questions/37532255/one-liner-to-transform-int-into-string
// arrayToString converts an array of float32 to a string, separated by the given separator, and returns the string
func arrayToString(nums []float32, sep string) string {
	return strings.Trim(strings.Join(strings.Fields(fmt.Sprint(nums)), sep), "[]")
}
