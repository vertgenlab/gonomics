package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"os"
	"strconv"
	"strings"
)

type AlleleCount struct {
	Ref 	dna.Base
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
}

// Map structure: map[Chromosome]map[Position]*AlleleCount
type SampleMap map[string]map[int64]*AlleleCount


// Inputs a sam file and loops through while keeping a tally of each base present at each position. Stores in SampleMap
func CountAlleles(refFilename string, samFilename string, minMapQ int64) SampleMap {

	// Read in reference
	fmt.Printf("#Reading Reference\n")
	ref := fasta.Read(refFilename)
	fasta.AllToUpper(ref)


	var i, k int32
	var j int
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	var ChrSliceMatch int

	sam.ReadHeader(samFile)

	// Initialize empty map of chromosomes to map of positions
	AlleleMatrix := make(map[string]map[int64]*AlleleCount)

	var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		if progressMeter % 500000 == 0 {
			fmt.Printf("#Read %d Alignments\n", progressMeter)
		}
		progressMeter++

		// If mapping quality is less than the threshold then go to next alignment
		if aln.MapQ < minMapQ {
			continue
		}

		// Identify the position in reference fasta slice of current chromosome
		for j = 0; j < len(ref); j++ {
			if ref[j].Name == aln.RName {
				ChrSliceMatch = j
			}
		}

		if aln.Cigar[0].Op != '*' {

			//if the chromosome has already been added to the matrix, move along
			_, ok := AlleleMatrix[aln.RName]

			//if the chromosome is NOT in the matrix, initialize
			if ! ok {
				AlleleMatrix[aln.RName] = make(map[int64]*AlleleCount)
			}

			SeqIndex = 0
			RefIndex = aln.Pos - 1

			for i = 0; i < int32(len(aln.Cigar)); i++ {
				currentSeq = aln.Seq

				//Handle deletion relative to ref
				//Each position deleted is annotated with a Del read
				if aln.Cigar[i].Op == 'D' {
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

						//if the position has already been added to the matrix, move along
						_, ok := AlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if ! ok {
							AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, 0}
						}

						AlleleMatrix[aln.RName][RefIndex].Del++
						RefIndex++
					}

					//Handle insertion relative to ref
					//The base after the inserted sequence is annotated with an Ins read
				} else if aln.Cigar[i].Op == 'I' {

					//if the position has already been added to the matrix, move along
					_, ok := AlleleMatrix[aln.RName][RefIndex]

					//if the position is NOT in the matrix, add it
					if ! ok {
						AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, 0}
					}

					AlleleMatrix[aln.RName][RefIndex].Ins++
					SeqIndex = SeqIndex + aln.Cigar[i].RunLength

					//Handle matching pos relative to ref
				} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

						//if the position has already been added to the matrix, move along
						_, ok := AlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if ! ok {
							AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, 0}
						}

						switch currentSeq[SeqIndex] {
						case dna.A:
							AlleleMatrix[aln.RName][RefIndex].BaseA++
						case dna.T:
							AlleleMatrix[aln.RName][RefIndex].BaseT++
						case dna.G:
							AlleleMatrix[aln.RName][RefIndex].BaseG++
						case dna.C:
							AlleleMatrix[aln.RName][RefIndex].BaseC++
						}
						SeqIndex++
						RefIndex++
					}
				} else if aln.Cigar[i].Op != 'H' {
					SeqIndex = SeqIndex + aln.Cigar[i].RunLength
				}
			}
		}
	}

	return AlleleMatrix
}

// Removes positions with insufficient coverage from the map
func FilterAlleles(input SampleMap, coverageThreshold int32) SampleMap {
	for chrName, chr := range input {
		for pos, alleles := range chr {

			coverage := alleles.BaseA + alleles.BaseC + alleles.BaseG + alleles.BaseT + alleles.Del

			if coverage <= coverageThreshold {
				delete(input[chrName], pos)
			}


		}
	}
	return input
}

// Write SampleMap to file
func WriteAlleleCounts(input SampleMap, output string) {

	var outFile *os.File

	if output != "stdout" {
		fmt.Printf("#Creating Output File\n")
		outFile, _ = os.Create(output)
		defer outFile.Close()
		io.WriteString(outFile, "Chr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\n")
	} else {
		fmt.Printf("#No Output Specified. Printing to STDOUT.\n")
		fmt.Printf("Chr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\n")
	}

	var base string
	var progressMeter int

	for chrName, chr := range input {
		for pos, alleles := range chr {

			switch alleles.Ref {
			case dna.A:
				base = "A"
			case dna.C:
				base = "C"
			case dna.G:
				base = "G"
			case dna.T:
				base = "T"
			case dna.N:
				base = "N"
			case dna.Gap:
				base = "Gap"
			case dna.Dot:
				base = "Dot"
			default:
				base = "NA"
			}

			if output != "stdout" {

				if progressMeter % 500000 == 0 {
					fmt.Printf("Wrote %d Positions\n", progressMeter)
				}
				progressMeter++

				// Write to file
				fmt.Fprintf(outFile,
					"%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
					chrName,
					pos + 1,
					base,
					alleles.BaseA,
					alleles.BaseC,
					alleles.BaseG,
					alleles.BaseT,
					alleles.Ins,
					alleles.Del)
			} else {
				// Write to stdout
				fmt.Printf(
					"%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
					chrName,
					pos + 1,
					base,
					alleles.BaseA,
					alleles.BaseC,
					alleles.BaseG,
					alleles.BaseT,
					alleles.Ins,
					alleles.Del)
			}
		}
	}
}

// Reading file from WriteAlleleCounts and store as a SampleMap
func ReadAlleleCounts(inFilename string) SampleMap {
	var line string
	var answer SampleMap
	var doneReading = false

	file := fileio.EasyOpen(inFilename)
	defer file.Close()

	// Initialize SampleMap
	answer = make(map[string]map[int64]*AlleleCount)

	var progressMeter int

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {

		if progressMeter % 500000 == 0 {
			fmt.Printf("# Read %d Lines\n", progressMeter)
		}
		progressMeter++

		//order is: Chr, Pos, Ref, A, C, G, T, Ins, Del
		words := strings.Split(line, "\t")

		//exclude header line
		if words[1] != "Pos"{

			// Define each entry to be stored in the map
			var Ref dna.Base
			switch words[2] {
			case "A":
				Ref = dna.A
			case "C":
				Ref = dna.C
			case "G":
				Ref = dna.G
			case "T":
				Ref = dna.T
			case "Gap":
				Ref = dna.Gap
			case "Dot":
				Ref = dna.Dot
			}

			var Chr = words[0]
			var Pos, _ = strconv.ParseInt(words[1], 10, 64)
			var BaseA, _ = strconv.ParseInt(words[3], 10, 32)
			var BaseC, _ = strconv.ParseInt(words[4], 10, 32)
			var BaseG, _ = strconv.ParseInt(words[5], 10, 32)
			var BaseT, _ = strconv.ParseInt(words[6], 10, 32)
			var Ins, _ = strconv.ParseInt(words[7], 10, 32)
			var Del, _ = strconv.ParseInt(words[8], 10, 32)


			//if the chromosome has already been added to the matrix, move along
			_, ok := answer[Chr]

			//if the chromosome is NOT in the matrix, initialize
			if ! ok {
				answer[Chr] = make(map[int64]*AlleleCount)
			}

			answer[Chr][Pos-1] = &AlleleCount{
				Ref:	Ref,
				BaseA:	int32(BaseA),
				BaseC:	int32(BaseC),
				BaseG:	int32(BaseG),
				BaseT:	int32(BaseT),
				Ins:	int32(Ins),
				Del:	int32(Del)}
		}

	}

	return answer
}

// Find the allele with with the highest frequency within a subset of 5 alleles (helper for FindMinorAllele)
func MaxMinorAllele(allele1 int32, allele2 int32, allele3 int32, allele4 int32, allele5 int32) int32 {
	var minorAllele = allele1
	if allele2 > minorAllele {minorAllele = allele2}
	if allele3 > minorAllele {minorAllele = allele3}
	if allele4 > minorAllele {minorAllele = allele4}
	if allele5 > minorAllele {minorAllele = allele5}
	return minorAllele
}

// Find the allele with highest frequency
func FindMajorAllele(A int32, C int32, G int32, T int32, Ins int32, Del int32) int32{
	var majorAllele = A
	if C > majorAllele {majorAllele = C}
	if G > majorAllele {majorAllele = G}
	if T > majorAllele {majorAllele = T}
	if Ins > majorAllele {majorAllele = Ins}
	if Del > majorAllele {majorAllele = Del}

	return majorAllele
}

// Find the allele with the 2nd highest frequency
func FindMinorAllele(A int32, C int32, G int32, T int32, Ins int32, Del int32) int32 {

	majorAllele := FindMajorAllele(A, C, G, T, Ins, Del)

	var minorAllele = A
	if majorAllele == A {minorAllele = MaxMinorAllele(C, G, T, Ins, Del)}
	if majorAllele == C {minorAllele = MaxMinorAllele(A, G, T, Ins, Del)}
	if majorAllele == G {minorAllele = MaxMinorAllele(A, C, T, Ins, Del)}
	if majorAllele == T {minorAllele = MaxMinorAllele(A, C, G, Ins, Del)}
	if majorAllele == Ins {minorAllele = MaxMinorAllele(A, C, G, T, Del)}
	if majorAllele == Del {minorAllele = MaxMinorAllele(A, C, G, T, Ins)}

	return minorAllele
}