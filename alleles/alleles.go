package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strconv"
	"strings"
)
// TODO: for progress meter make it log.Printf instead of fmt.Printf to get timestamps
// TODO: look into map of struct {chr pos} as key instead of map of maps
// TODO: in vcf format multiple possible alleles into a single line
type AlleleCount struct {
	Ref    dna.Base
	Counts int32
	BaseA  int32
	BaseC  int32
	BaseG  int32
	BaseT  int32
	Ins    []Indel
	Del    []Indel
}

type Indel struct {
	Ref       []dna.Base
	Alt       []dna.Base
	RunLength int32
	Count     int32
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
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64
	var Match bool

	sam.ReadHeader(samFile)

	// Initialize empty map of chromosomes to map of positions
	AlleleMatrix := make(map[string]map[int64]*AlleleCount)

	var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		if progressMeter%500000 == 0 {
			log.Printf("#Read %d Alignments\n", progressMeter)
		}
		progressMeter++

		if aln.Cigar[0].Op != '*' {
		// If mapping quality is less than the threshold then go to next alignment
		if aln.MapQ < minMapQ {
			continue
		}

		// TODO: look into fasta map so this is unnecessary
		// Identify the position in reference fasta slice of current chromosome
		for j = 0; j < len(ref); j++ {
			if ref[j].Name == aln.RName {
				ChrSliceMatch = j
			}
		}

			//if the chromosome has already been added to the matrix, move along
			_, ok := AlleleMatrix[aln.RName]

			//if the chromosome is NOT in the matrix, initialize
			if !ok {
				AlleleMatrix[aln.RName] = make(map[int64]*AlleleCount)
			}

			SeqIndex = 0
			RefIndex = aln.Pos - 1

			for i = 0; i < int32(len(aln.Cigar)); i++ {
				currentSeq = aln.Seq

				//Handle deletion relative to ref
				//Each position deleted is annotated with a Del read
				if aln.Cigar[i].Op == 'D' {
					//fmt.Printf("DELETION\n")
					OrigRefIndex = RefIndex
					indelSeq = make([]dna.Base, 1)

					// First base in indel is the base prior to the indel sequence per VCF standard format
					indelSeq[0] = ref[ChrSliceMatch].Seq[OrigRefIndex-1]

					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

						//if the position has already been added to the matrix, move along
						_, ok := AlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if !ok {
							AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{
								ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0), make([]Indel, 0)}
						}

						// Keep track of deleted sequence
						indelSeq = append(indelSeq, ref[ChrSliceMatch].Seq[RefIndex])

						//AlleleMatrix[aln.RName][RefIndex].Del++
						AlleleMatrix[aln.RName][RefIndex].Counts++
						RefIndex++
					}

					Match = false
					for j = 0; j < len(AlleleMatrix[aln.RName][OrigRefIndex].Del); j++ {
						// If the deletion has already been seen before, increment the existing entry
						if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMatrix[aln.RName][OrigRefIndex].Del[j].Ref) == 0 {
							AlleleMatrix[aln.RName][OrigRefIndex].Del[j].Count++
							Match = true
							break
						}
					}

					// If the deletion has not been seen before, then append it to the Del slice
					if Match == false {
						var AltBase []dna.Base = make([]dna.Base, 1)
						AltBase[0] = indelSeq[0]
						currentIndel = Indel{indelSeq, AltBase, int32(len(indelSeq) - 1), 1}
						AlleleMatrix[aln.RName][OrigRefIndex].Del = append(AlleleMatrix[aln.RName][OrigRefIndex].Del, currentIndel)
					}

					//Handle insertion relative to ref
					//The base after the inserted sequence is annotated with an Ins read
				} else if aln.Cigar[i].Op == 'I' {
					//fmt.Printf("INSERTION\n")
					//if the position has already been added to the matrix, move along
					_, ok := AlleleMatrix[aln.RName][RefIndex]

					//if the position is NOT in the matrix, add it
					if !ok {
						AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{
							ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0), make([]Indel, 0)}
					}

					// Loop through read sequence and keep track of the inserted bases
					indelSeq = make([]dna.Base, 1)

					// First base in indel is the base prior to the indel sequence per VCF standard format
					indelSeq[0] = ref[ChrSliceMatch].Seq[RefIndex-1]

					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
						indelSeq = append(indelSeq, currentSeq[SeqIndex])
						SeqIndex++
					}

					Match = false
					for j = 0; j < len(AlleleMatrix[aln.RName][RefIndex].Ins); j++ {
						if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMatrix[aln.RName][RefIndex].Ins[j].Alt) == 0 {
							// If they the inserted sequence matches a previously inserted sequence, then increment the count
							AlleleMatrix[aln.RName][RefIndex].Ins[j].Count++
							Match = true
						}
					}

					if Match == false {
						var RefBase []dna.Base = make([]dna.Base, 1)
						RefBase[0] = indelSeq[0]
						currentIndel = Indel{RefBase, indelSeq, int32(len(indelSeq)), 1}
						AlleleMatrix[aln.RName][RefIndex].Ins = append(AlleleMatrix[aln.RName][RefIndex].Ins, currentIndel)
					}

					// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

					//Handle matching pos relative to ref
				} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {
					//fmt.Printf("BASE\n")
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

						//if the position has already been added to the matrix, move along
						_, ok := AlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if !ok {
							AlleleMatrix[aln.RName][RefIndex] = &AlleleCount{
								ref[ChrSliceMatch].Seq[RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0), make([]Indel, 0)}
						}

						switch currentSeq[SeqIndex] {
						case dna.A:
							AlleleMatrix[aln.RName][RefIndex].BaseA++
							AlleleMatrix[aln.RName][RefIndex].Counts++
						case dna.T:
							AlleleMatrix[aln.RName][RefIndex].BaseT++
							AlleleMatrix[aln.RName][RefIndex].Counts++
						case dna.G:
							AlleleMatrix[aln.RName][RefIndex].BaseG++
							AlleleMatrix[aln.RName][RefIndex].Counts++
						case dna.C:
							AlleleMatrix[aln.RName][RefIndex].BaseC++
							AlleleMatrix[aln.RName][RefIndex].Counts++
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

// Convert SampleMap to VCF
func AllelesToVcf(input SampleMap) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var current *vcf.Vcf
	var UknFmt string
	var base string
	var RefCount int32
	var i int
	var RefSeq string
	var AltSeq string

	for chrName, chr := range input {
		for pos, alleles := range chr {

			switch alleles.Ref {
			case dna.A:
				base = "A"
				RefCount = alleles.BaseA
			case dna.C:
				base = "C"
				RefCount = alleles.BaseC
			case dna.G:
				base = "G"
				RefCount = alleles.BaseG
			case dna.T:
				base = "T"
				RefCount = alleles.BaseT
			default:
				base = "NA"
				RefCount = 0
			}

			// Ref -> A
			UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseA, alleles.Counts)

			current = &vcf.Vcf{
				Chr:     chrName,
				Pos:     pos + 1,
				Id:      ".",
				Ref:     base,
				Alt:     "A",
				Qual:    1,
				Filter:  ".",
				Info:    ".",
				Format:  "RefCount:AltCount:Cov",
				Unknown: UknFmt}

			answer = append(answer, current)

			// Ref -> C
			UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseC, alleles.Counts)

			current = &vcf.Vcf{
				Chr:     chrName,
				Pos:     pos + 1,
				Id:      ".",
				Ref:     base,
				Alt:     "C",
				Qual:    1,
				Filter:  ".",
				Info:    ".",
				Format:  "RefCount:AltCount:Cov",
				Unknown: UknFmt}

			answer = append(answer, current)

			// Ref -> G
			UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseG, alleles.Counts)

			current = &vcf.Vcf{
				Chr:     chrName,
				Pos:     pos + 1,
				Id:      ".",
				Ref:     base,
				Alt:     "G",
				Qual:    1,
				Filter:  ".",
				Info:    ".",
				Format:  "RefCount:AltCount:Cov",
				Unknown: UknFmt}

			answer = append(answer, current)

			// Ref -> T
			UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseT, alleles.Counts)

			current = &vcf.Vcf{
				Chr:     chrName,
				Pos:     pos + 1,
				Id:      ".",
				Ref:     base,
				Alt:     "T",
				Qual:    1,
				Filter:  ".",
				Info:    ".",
				Format:  "RefCount:AltCount:Cov",
				Unknown: UknFmt}

			answer = append(answer, current)

			// Ref -> Ins
			for i = 0; i < len(alleles.Ins); i++ {
				UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.Ins[i].Count, alleles.Counts)

				RefSeq = dna.BaseToString(alleles.Ins[i].Ref[0])
				AltSeq = dna.BasesToString(alleles.Ins[i].Alt)

				current = &vcf.Vcf{
					Chr: chrName,
					// VCF format has insertions assigned to the prior base, so there is no pos + 1
					Pos:     pos,
					Id:      ".",
					Ref:     RefSeq,
					Alt:     AltSeq,
					Qual:    1,
					Filter:  ".",
					Info:    "INS",
					Format:  "RefCount:AltCount:Cov",
					Unknown: UknFmt}

				answer = append(answer, current)
			}

			// Ref -> Del
			for i = 0; i < len(alleles.Del); i++ {
				UknFmt = fmt.Sprintf("%d:%d:%d", RefCount, alleles.Del[i].Count, alleles.Counts)

				RefSeq = dna.BasesToString(alleles.Del[i].Ref)
				AltSeq = dna.BaseToString(alleles.Del[i].Alt[0])

				current = &vcf.Vcf{
					Chr: chrName,
					// VCF format has deletions assigned to the prior base, so there is no pos + 1
					Pos:     pos,
					Id:      ".",
					Ref:     RefSeq,
					Alt:     AltSeq,
					Qual:    1,
					Filter:  ".",
					Info:    "DEL",
					Format:  "RefCount:AltCount:Cov",
					Unknown: UknFmt}

				answer = append(answer, current)
			}
		}
	}
	return answer
}

// Removes positions with insufficient coverage from the map
func FilterAlleles(input SampleMap, coverageThreshold int32) SampleMap {
	for chrName, chr := range input {
		for pos, alleles := range chr {

			if alleles.Counts < coverageThreshold {
				delete(input[chrName], pos)
			}

		}
	}
	return input
}

// Readin VCF file from AllelesToVcf -> vcf.Write and store as a SampleMap
func ReadVcfToAlleleCounts(inFilename string) SampleMap {
	var line string
	var answer SampleMap
	var doneReading = false
	var RefSeq, AltSeq []dna.Base
	var UknFmt []string
	var AltCount, Counts int64
	var currentIndel Indel

	file := fileio.EasyOpen(inFilename)
	defer file.Close()

	// Initialize SampleMap
	answer = make(map[string]map[int64]*AlleleCount)

	var progressMeter int

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {

		if progressMeter%500000 == 0 {
			fmt.Printf("# Read %d Lines\n", progressMeter)
		}
		progressMeter++

		//order is: Chr, Pos, Id, Ref, Alt, Qual, Filter, Info, Format, Unknown (RefCount:AltCount:Coverage)
		words := strings.Split(line, "\t")

		//exclude header line
		if strings.HasPrefix(words[0], "#") == false {

			// Define each entry to be stored in the map
			var Chr string = words[0]

			//if the chromosome has already been added to the matrix, move along
			_, ok := answer[Chr]

			//if the chromosome is NOT in the matrix, initialize
			if !ok {
				answer[Chr] = make(map[int64]*AlleleCount)
			}

			// Convert strings to stored values
			RefSeq = dna.StringToBases(words[3])
			AltSeq = dna.StringToBases(words[4])
			UknFmt = strings.Split(words[9], ":")
			AltCount, _ = strconv.ParseInt(UknFmt[1], 10, 32)
			Counts, _ = strconv.ParseInt(UknFmt[2], 10, 32)

			// If point mutation
			if len(RefSeq) == 1 && len(AltSeq) == 1 {

				var Pos, _ = strconv.ParseInt(words[1], 10, 64)

				// If the position is in the map move along, else initialize
				// Subtract 1 from Pos for index 0
				_, okk := answer[Chr][Pos-1]
				if !okk {
					answer[Chr][Pos-1] = &AlleleCount{
						Ref:    0,
						Counts: 0,
						BaseA:  0,
						BaseC:  0,
						BaseG:  0,
						BaseT:  0,
						Ins:    make([]Indel, 0),
						Del:    make([]Indel, 0)}
				}

				answer[Chr][Pos-1].Ref = RefSeq[0]
				answer[Chr][Pos-1].Counts = int32(Counts)

				switch AltSeq[0] {
				case dna.A:
					answer[Chr][Pos-1].BaseA = int32(AltCount)
				case dna.C:
					answer[Chr][Pos-1].BaseC = int32(AltCount)
				case dna.G:
					answer[Chr][Pos-1].BaseG = int32(AltCount)
				case dna.T:
					answer[Chr][Pos-1].BaseT = int32(AltCount)
				}

				// If indel
			} else if len(RefSeq) != len(AltSeq) {

				var Pos, _ = strconv.ParseInt(words[1], 10, 64)

				// If the position is in the map move along, else initialize
				// VCF stores pos as base prior to indel so subtracting 1 for index 0 is unnecessary
				_, okk := answer[Chr][Pos]
				if !okk {
					answer[Chr][Pos] = &AlleleCount{
						Ref:    0,
						Counts: 0,
						BaseA:  0,
						BaseC:  0,
						BaseG:  0,
						BaseT:  0,
						Ins:    make([]Indel, 0),
						Del:    make([]Indel, 0)}
				}

				// If Ins
				if len(RefSeq) < len(AltSeq) {

					currentIndel = Indel{
						Ref:       RefSeq,
						Alt:       AltSeq,
						RunLength: int32(len(AltSeq[1:])),
						Count:     int32(AltCount)}

					answer[Chr][Pos].Ins = append(answer[Chr][Pos].Ins, currentIndel)
				}

				// If Del
				if len(RefSeq) > len(AltSeq) {

					currentIndel = Indel{
						Ref:       RefSeq,
						Alt:       AltSeq,
						RunLength: int32(len(RefSeq[1:])),
						Count:     int32(AltCount)}

					answer[Chr][Pos].Del = append(answer[Chr][Pos].Del, currentIndel)
				}
			}
		}
	}
	return answer
}

// Find the allele with with the highest frequency within a subset of 5 alleles (helper for FindMinorAllele)
func maxMinorAllele(allele1 int32, allele2 int32, allele3 int32, allele4 int32, allele5 int32) int32 {
	var minorAllele = allele1
	if allele2 > minorAllele {
		minorAllele = allele2
	}
	if allele3 > minorAllele {
		minorAllele = allele3
	}
	if allele4 > minorAllele {
		minorAllele = allele4
	}
	if allele5 > minorAllele {
		minorAllele = allele5
	}
	return minorAllele
}

// Find the allele with highest frequency
func FindMajorAllele(A int32, C int32, G int32, T int32, Ins int32, Del int32) int32 {
	var majorAllele = A
	if C > majorAllele {
		majorAllele = C
	}
	if G > majorAllele {
		majorAllele = G
	}
	if T > majorAllele {
		majorAllele = T
	}
	if Ins > majorAllele {
		majorAllele = Ins
	}
	if Del > majorAllele {
		majorAllele = Del
	}

	return majorAllele
}

// Find the allele with the 2nd highest frequency
func FindMinorAllele(A int32, C int32, G int32, T int32, Ins int32, Del int32) int32 {

	majorAllele := FindMajorAllele(A, C, G, T, Ins, Del)

	var minorAllele = A
	if majorAllele == A {
		minorAllele = maxMinorAllele(C, G, T, Ins, Del)
	}
	if majorAllele == C {
		minorAllele = maxMinorAllele(A, G, T, Ins, Del)
	}
	if majorAllele == G {
		minorAllele = maxMinorAllele(A, C, T, Ins, Del)
	}
	if majorAllele == T {
		minorAllele = maxMinorAllele(A, C, G, Ins, Del)
	}
	if majorAllele == Ins {
		minorAllele = maxMinorAllele(A, C, G, T, Del)
	}
	if majorAllele == Del {
		minorAllele = maxMinorAllele(A, C, G, T, Ins)
	}

	return minorAllele
}
