package browser

import (
	"fmt"
	"log"
	"strings"
	"unicode/utf8"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
)

// PFaVisualizer produces command line visualizations of pFasta format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes lineLength formatting option as int.
// If 0 sig figs or 0 decimal places, returns full probability
// There will be slight floating point errors in the last (or last two) places
func PFaVisualizer(infile string, outfile string, start int, end int, startOfAlignment bool, endOfAlignment bool, sigFigs int, decimalPlaces int, lineLength int, seqName string) {
	if !startOfAlignment {
		if !endOfAlignment && !(start < end) {
			log.Fatalf("Error: Invalid arguments, start must be lower than end\n")
		} else if start < 0 {
			log.Fatalf("Error: Invalid arguments, start must be greater or equal to 0, or (case-insensitive) 'start'\n")
		}
	}

	records := pFasta.Read(infile)
	out := fileio.EasyCreate(outfile)

	if len(records) > 1 && seqName == "" {
		log.Fatalf("Error: Invalid arguments, must provide sequence name for file with multiple pFastas.\n")
	}

	var err error
	var formatting string
	var formatNum int
	if sigFigs == 0 {
		formatting = "Decimal Places"
		formatNum = decimalPlaces
	} else {
		formatting = "SigFigs"
		formatNum = sigFigs
	}

	if startOfAlignment {
		start = 0
	}

	if len(records) == 0 {
		log.Fatalf("Error: User provided empty pfasta file.\n")
	}
	if seqName == "" {
		// only accepts single pFa or a .pFa with multiple pFastas and a specified record name
		if len(records) > 1 {
			log.Fatalf("Error: User must specify sequence name for pFasta file with more than 1 sequence.\n")
		} else {
			// pfa with 1 entry

			if endOfAlignment {
				end = len(records[0].Seq)
			}

			_, err = fmt.Fprintf(out, "Start: %d. End: %d. %s: %d.", start, end, formatting, formatNum)
			exception.PanicOnErr(err)

			printAllSets(out, err, records[0], start, end, lineLength, sigFigs, decimalPlaces)
		}
	} else {
		// user can specify chrom if multiple entries
		seqPrinted := false
		for desiredSeqIdx, desiredSeq := range records {
			if desiredSeq.Name == seqName {
				seqPrinted = true

				if endOfAlignment {
					end = len(records[desiredSeqIdx].Seq)
				}

				_, err = fmt.Fprintf(out, "Start: %d. End: %d. %s: %d.", start, end, formatting, formatNum)
				exception.PanicOnErr(err)

				printAllSets(out, err, desiredSeq, start, end, lineLength, sigFigs, decimalPlaces)
				break
			}
		}
		if !seqPrinted {
			log.Fatalf("Error: User specified sequence not in input pfasta file.")
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// printAllSets prints probability distribution of bases from pos start to pos end in record
func printAllSets(out *fileio.EasyWriter, err error, record pFasta.PFasta, start int, end int, lineLength int, sigFigs int, decimalPlaces int) {
	if end == -1 {
		end = len(record.Seq)
	}
	lineA := make([]float32, lineLength)
	lineC := make([]float32, lineLength)
	lineG := make([]float32, lineLength)
	lineT := make([]float32, lineLength)
	setOfLinesIdx := 0

	for setOfLinesIdx = 0; setOfLinesIdx < (end-start)/lineLength; setOfLinesIdx++ {
		_, err = fmt.Fprintf(out, "\n")
		exception.PanicOnErr(err)
		printOneSetLines(lineLength, setOfLinesIdx, lineLength, lineA, lineC, lineG, lineT, start, record, out, sigFigs, decimalPlaces, len(record.Name))
	}
	_, err = fmt.Fprintf(out, "\n")
	exception.PanicOnErr(err)
	printOneSetLines(lineLength, setOfLinesIdx, (end-start)%lineLength, lineA, lineC, lineG, lineT, start, record, out, sigFigs, decimalPlaces, len(record.Name))
}

// printOneSetLines prints probability distribution of bases from init_pos =(setOfLinesIdx*lineLength + start) in record to (init_pos + numIters)
func printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, record pFasta.PFasta, out *fileio.EasyWriter, sigFigs int, decimalPlaces int, longestName int) {
	recordIdx := setOfLinesIdx*lineLength + start
	lineIdx := 0
	var err error
	_, err = fmt.Fprintf(out, "Position: %d\n", recordIdx)
	exception.PanicOnErr(err)
	for lineIdx = 0; lineIdx < numIters; lineIdx++ {
		lineA[lineIdx], lineC[lineIdx], lineG[lineIdx], lineT[lineIdx] = record.Seq[recordIdx+lineIdx].A, record.Seq[recordIdx+lineIdx].C, record.Seq[recordIdx+lineIdx].G, record.Seq[recordIdx+lineIdx].T
	}
	if sigFigs == 0 {
		printOneBaseDigits(lineA[0:lineIdx], "A", out, longestName, record.Name, decimalPlaces)
		printOneBaseDigits(lineC[0:lineIdx], "C", out, longestName, "", decimalPlaces)
		printOneBaseDigits(lineG[0:lineIdx], "G", out, longestName, "", decimalPlaces)
		printOneBaseDigits(lineT[0:lineIdx], "T", out, longestName, "", decimalPlaces)
	} else {
		printOneBaseScientific(lineA[0:lineIdx], "A", out, longestName, record.Name, sigFigs)
		printOneBaseScientific(lineC[0:lineIdx], "C", out, longestName, "", sigFigs)
		printOneBaseScientific(lineG[0:lineIdx], "G", out, longestName, "", sigFigs)
		printOneBaseScientific(lineT[0:lineIdx], "T", out, longestName, "", sigFigs)
	}
}

func printOneBaseScientific(line []float32, base string, out *fileio.EasyWriter, longestName int, name string, sigFigs int) {
	var answer strings.Builder
	for _, base := range line {
		answer.WriteString(fmt.Sprintf("\t%.*e", sigFigs-1, base))
	}
	var err error
	_, err = fmt.Fprintf(out, ">%-*s\t|\t%s\t|%s\n", longestName, name, base, answer.String())
	exception.PanicOnErr(err)
}

func printOneBaseDigits(line []float32, base string, out *fileio.EasyWriter, longestName int, name string, decimalPlaces int) {
	var answer strings.Builder
	for _, base := range line {
		answer.WriteString(fmt.Sprintf("\t%.*f", decimalPlaces, base))
	}
	var err error
	_, err = fmt.Fprintf(out, ">%-*s\t|\t%s\t|%s\n", longestName, name, base, answer.String())
	exception.PanicOnErr(err)
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

// PFaVisualizerTsv produces command line visualizations of pFasta format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes lineLength formatting option as int.
// If 0 sig figs or 0 decimal places, returns full probability
// There will be slight floating point errors in the last (or last two) places
func PFaVisualizerTsv(infile string, outfile string, start int, end int, startOfAlignment bool, endOfAlignment bool, sigFigs int, decimalPlaces int, lineLength int, seqName string) {
	if !startOfAlignment {
		if !endOfAlignment && !(start < end) {
			log.Fatalf("Error: Invalid arguments, start must be lower than end\n")
		} else if start < 0 {
			log.Fatalf("Error: Invalid arguments, start must be greater or equal to 0, or (case-insensitive) 'start'\n")
		}
	}
	records := pFasta.Read(infile)
	out := fileio.EasyCreate(outfile)

	if len(records) > 1 && seqName == "" {
		log.Fatalf("Error: Invalid arguments, must provide sequence name for file with multiple pFastas.\n")
	}

	var err error

	if startOfAlignment {
		start = 0
	}
	if len(records) == 0 {
		log.Fatalf("Error: User provided empty pfasta file.\n")
	}

	if seqName == "" {
		// only accepts single pFa or a .pFa with multiple pFastas and a specified record name
		if len(records) > 1 {
			log.Fatalf("Error: User must specify sequence name for pFasta file with more than 1 sequence.\n")
		} else {
			// pfa with 1 entry
			if endOfAlignment {
				end = len(records[0].Seq)
			} else {
				end += 1
			}

			_, err = fmt.Fprintf(out, "Position\tBase\tProbability\n")
			exception.PanicOnErr(err)
			printAllSetsTsv(out, err, records[0], start, end, lineLength, sigFigs, decimalPlaces)
		}
	} else {
		// user can specify chrom if multiple entries
		seqPrinted := false
		for desiredSeqIdx, desiredSeq := range records {
			if desiredSeq.Name == seqName {
				seqPrinted = true

				if endOfAlignment {
					end = len(records[desiredSeqIdx].Seq)
				}
				_, err = fmt.Fprintf(out, "Position\tBase\tProbability\n")
				exception.PanicOnErr(err)
				printAllSetsTsv(out, err, desiredSeq, start, end, lineLength, sigFigs, decimalPlaces)

				break
			}
		}
		if !seqPrinted {
			log.Fatalf("Error: User specified sequence not in input pfasta file.")
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// printAllSetsTsv prints probability distribution of bases from pos start to pos end in record
func printAllSetsTsv(out *fileio.EasyWriter, err error, record pFasta.PFasta, start int, end int, lineLength int, sigFigs int, decimalPlaces int) {
	idx := start
	
	if sigFigs == 0 {
		// gives decimal places
		for idx < end {
			base := record.Seq[idx]
			_, err = fmt.Fprintf(out, "%v\tA\t%.*f\n", idx, decimalPlaces, base.A)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tC\t%.*f\n", idx, decimalPlaces, base.C)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tG\t%.*f\n", idx, decimalPlaces, base.G)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tT\t%.*f\n", idx, decimalPlaces, base.T)
			exception.PanicOnErr(err)
			idx++
		}
	} else {
		// gives sigfigs
		for idx < end {
			base := record.Seq[idx]
			_, err = fmt.Fprintf(out, "%v\tA\t%.*f\n", idx, sigFigs-1, base.A)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tC\t%.*e\n", idx, sigFigs-1, base.C)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tG\t%.*e\n", idx, sigFigs-1, base.G)
			exception.PanicOnErr(err)
			_, err = fmt.Fprintf(out, "%v\tT\t%.*e\n", idx, sigFigs-1, base.T)
			exception.PanicOnErr(err)
			idx++
		}
	}
}
