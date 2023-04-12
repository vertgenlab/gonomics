// Package browser contains command line visualization tools for genomic information.
package browser

import (
	"fmt"
	"log"
	"unicode/utf8"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
)

// MultiFaVisualizer produces command line visualizations of multiFa format alignments from a specified start and end position.
// Can be written to a file or to standard out. Includes noMask and lineLength formatting options as bools.
func MultiFaVisualizer(infile string, outfile string, start int, end int, noMask bool, lineLength int, endOfAlignment bool) {
	if !(start < end) && !endOfAlignment {
		log.Fatalf("Invalid arguments, start must be lower than end")
	}
	var stop int
	records := fasta.Read(infile)
	if noMask {
		fasta.AllToUpper(records)
	}
	for i := 1; i < len(records); i++ {
		for j := range records[0].Seq {
			if records[i].Seq[j] == records[0].Seq[j] {
				records[i].Seq[j] = dna.Dot
			}
		}
	}
	long := calculateLongestName(records)
	var refCounter, startCounter, endCounter int = 0, 0, 0
	for t := 0; refCounter < start; t++ {
		startCounter++
		if t == len(records[0].Seq) {
			log.Fatalf("Ran out of chromosome")
		} else if records[0].Seq[t] != dna.Gap {
			refCounter++
		}
	}
	chromStart := refCounter

	out := fileio.EasyCreate(outfile)
	defer out.Close()

	fmt.Fprintf(out, "Start: %d. refCounter: %d. alignCounter: %d\n", start, refCounter, startCounter)

	refCounter = 0

	if endOfAlignment {
		endCounter = len(records[0].Seq)
	} else {
		for n := 0; refCounter < end; n++ {
			endCounter++
			if n == len(records[0].Seq) {
				log.Fatalf("Ran off the chromosome")
			} else if records[0].Seq[n] != dna.Gap {
				refCounter++
			}
		}
	}

	for k := startCounter; k < endCounter; k = k + lineLength {
		fmt.Fprintf(out, "Position: %d\n", chromStart)
		stop = numbers.Min(endCounter, k+lineLength)
		for m := range records {
			fmt.Fprintf(out, "|%-*s| %s\n", long, records[m].Name, dna.BasesToString(records[m].Seq[k:stop]))
		}
		fmt.Fprintf(out, "\n\n")
		chromStart = chromStart + lineLength - dna.CountGaps(records[0].Seq[k:stop])
	}
}

// calculateLongestName is a helper function of MultiFaVisualizer that returns the length of the longest name in a slice of fasta.Fasta structs.
func calculateLongestName(f []fasta.Fasta) int {
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
