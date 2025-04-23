package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sort"
)

// AssemblyStats takes the path to a fasta file and a flag for whether lower case letters
// should count as assembly gaps.  Six ints are returned, which encode:
// the N50 size, the L50 size, half the size of the genome, size of the genome, size of the largest contig, and the number of contigs.
func AssemblyStats(infile string, countLowerAsGaps bool) (int, int, int, int, int, int) {
	records := Read(infile)

	contigList := MakeContigList(records, countLowerAsGaps)

	genomeLength := calculateGenomeLength(contigList)

	sort.Ints(contigList)
	halfGenome := genomeLength / 2

	N50, L50 := CalculateN50L50(contigList, halfGenome)
	numContigs := len(contigList)
	largestContig := contigList[len(contigList)-1]

	return N50, L50, halfGenome, genomeLength, largestContig, numContigs
}

// calculateGenomeLength is a helper function that takes in a slice of int representing contigs in an assembly and returns the genome size (int)
func calculateGenomeLength(contigList []int) int {
	if len(contigList) == 0 {
		log.Fatalf("Error: Cannot calculate genome length -- contig list is empty")
	}
	var genomeLength int
	for i := 0; i < len(contigList); i++ {
		genomeLength += contigList[i]
	}
	return genomeLength
}

// CalculateN50L50 takes a sorted slice of contig lengths and the size of half the genome. It returns the N50 size and L50 size.
func CalculateN50L50(contigList []int, halfGenome int) (N50 int, L50 int) {
	var sum int = 0
	L50 = 0
	for i := len(contigList) - 1; i > -1; i-- {
		L50++
		sum += contigList[i]
		if sum >= halfGenome {
			return contigList[i], L50
		}
	}
	log.Fatalf("Unable to calculate N50 / L50.")
	return -1, -1
}

// MakeContigList takes a slice of fasta sequences and a flag for whether lower case
// letters should count as gaps.  A slice of contig sizes is the return value.
func MakeContigList(records []Fasta, countLowerAsGaps bool) []int {
	var contigLen int = 0
	var contig bool = false
	var contigList []int

	for i := 0; i < len(records); i++ {
		contig = false
		for k := 0; k < len(records[i].Seq); k++ {
			if contig {
				if countLowerAsGaps {
					if records[i].Seq[k] == dna.N || dna.IsLower(records[i].Seq[k]) {
						contig = false
						contigList = append(contigList, contigLen)
						contigLen = 0
					} else {
						contigLen++
					}
				} else {
					if records[i].Seq[k] == dna.N {
						contig = false
						contigList = append(contigList, contigLen)
						contigLen = 0
					} else {
						contigLen++
					}
				}
			} else {
				if !(records[i].Seq[k] == dna.N) {
					if countLowerAsGaps {
						if !(dna.IsLower(records[i].Seq[k])) {
							contig = true
							contigLen++
						}
					} else {
						contig = true
						contigLen++
					}
				}
			}
		}

		if contig {
			contigList = append(contigList, contigLen)
			contigLen = 0
		}
	}
	return contigList
}

// WriteAssemblyStats takes the name of an assembly, a path to an output file, and stats for:
// the N50 size, half the size of the genome, size of the genome, size of the largest contig, and the number of contigs.
// The stats, with some human-readable labels are written to the output file.
func WriteAssemblyStats(assemblyName string, outfile string, N50 int, L50 int, halfGenome int, genomeLength int, largestContig int, numContigs int) {
	file := fileio.EasyCreate(outfile)
	var err error
	_, err = fmt.Fprintf(file, "Assembly Name: %s\n", assemblyName)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "halfGenome: %d\n", halfGenome)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "genomeLength: %d\n", genomeLength)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "Number of contigs: %d\n", numContigs)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "Largest Contig: %d\n", largestContig)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "N50: %d\n", N50)
	exception.FatalOnErr(err)
	_, err = fmt.Fprintf(file, "L50: %d\n", L50)
	exception.FatalOnErr(err)

	err = file.Close()
	exception.FatalOnErr(err)
}
