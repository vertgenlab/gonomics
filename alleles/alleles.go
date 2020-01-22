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
// TODO: in vcf format multiple possible alleles into a single line
type AlleleCount struct {
	Ref    	dna.Base
	Counts 	int32
	BaseA  	int32
	BaseC 	int32
	BaseG 	int32
	BaseT  	int32
	Indel	[]Indel
}

type Indel struct {
	Ref       []dna.Base
	Alt       []dna.Base
	Count     int32
}

type Location struct {
	Chr 	string
	Pos 	int64
}

// Map structure: map[Chromosome]map[Position]*AlleleCount
type SampleMap map[Location]*AlleleCount

// Temporary function until this is added to the fasta package
func refToMap(refFilename string) map[string][]dna.Base {
	inRef := fasta.Read(refFilename)
	fasta.AllToUpper(inRef)
	ref := make(map[string][]dna.Base)
	var curr *fasta.Fasta
	for i := 0; i < len(inRef); i++ {
		curr = inRef[i]
		_, ok := ref[curr.Name]
		if !ok {
			ref[curr.Name] = curr.Seq
		}
	}
	return ref
}

// Inputs a sam file and loops through while keeping a tally of each base present at each position. Stores in SampleMap
func CountAlleles(refFilename string, samFilename string, minMapQ int64) SampleMap {
	// Read in reference
	fmt.Printf("#Reading Reference\n")
	ref := refToMap(refFilename)

	var i, k int32
	var j int
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64
	var Match bool

	sam.ReadHeader(samFile)

	// Initialize empty map of chromosomes to map of positions
	AlleleMap := make(SampleMap)

	var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		if progressMeter%500000 == 0 {
			log.Printf("#Read %d Alignments\n", progressMeter)
		}
		progressMeter++

		// If read is unmapped then go to the next alignment
		if aln.Cigar[0].Op == '*' {
			continue
		}

		// If mapping quality is less than the threshold then go to next alignment
		if aln.MapQ < minMapQ {
			continue
		}

		SeqIndex = 0
		RefIndex = aln.Pos - 1

		for i = 0; i < int32(len(aln.Cigar)); i++ {
			currentSeq = aln.Seq

			//Handle deletion relative to ref
			//Each position deleted is annotated with counts + 1
			if aln.Cigar[i].Op == 'D' {
				OrigRefIndex = RefIndex
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				indelSeq[0] = ref[aln.RName][OrigRefIndex-1]

				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					// If the position has already been added to the map, move along
					_, ok := AlleleMap[Location{aln.RName, RefIndex}]

					// If the position is NOT in the map, initialize
					if !ok {
						AlleleMap[Location{aln.RName, RefIndex}] = &AlleleCount{
							ref[aln.RName][RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0)}
					}

					// Keep track of deleted sequence
					indelSeq = append(indelSeq, ref[aln.RName][RefIndex])

					AlleleMap[Location{aln.RName,RefIndex}].Counts++
					RefIndex++
				}

				Match = false
				for j = 0; j < len(AlleleMap[Location{aln.RName,OrigRefIndex}].Indel); j++ {
					// If the deletion has already been seen before, increment the existing entry
					// For a deletion the indelSeq should match the Ref
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Ref) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Alt) == 0{
						AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Count++
						Match = true
						break
					}
				}

				// If the deletion has not been seen before, then append it to the Del slice
				// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
				if Match == false {
					currentIndel = Indel{indelSeq, indelSeq[:1], 1}
					AlleleMap[Location{aln.RName, OrigRefIndex}].Indel = append(AlleleMap[Location{aln.RName, OrigRefIndex}].Indel, currentIndel)
				}

				//Handle insertion relative to ref
				//The base after the inserted sequence is annotated with an Ins read
			} else if aln.Cigar[i].Op == 'I' {

				// If the position has already been added to the map, move along
				_, ok := AlleleMap[Location{aln.RName, RefIndex}]

				// If the position is NOT in the map, initialize
				if !ok {
					AlleleMap[Location{aln.RName, RefIndex}] = &AlleleCount{
						ref[aln.RName][RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0)}
				}

				// Loop through read sequence and keep track of the inserted bases
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				indelSeq[0] = ref[aln.RName][RefIndex-1]

				// Keep track of inserted sequence by moving along the read
				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
					indelSeq = append(indelSeq, currentSeq[SeqIndex])
					SeqIndex++
				}

				Match = false
				for j = 0; j < len(AlleleMap[Location{aln.RName, RefIndex}].Indel); j++ {
					// If they the inserted sequence matches a previously inserted sequence, then increment the count
					// For an insertion, the indelSeq should match the Alt
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[1:], AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Ref) == 0 {
						AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Count++
						Match = true
						break
					}
				}

				if Match == false {
					currentIndel = Indel{indelSeq[:1], indelSeq, 1}
					AlleleMap[Location{aln.RName, RefIndex}].Indel = append(AlleleMap[Location{aln.RName, RefIndex}].Indel, currentIndel)
				}

				// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

				//Handle matching pos relative to ref
			} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					//if the position has already been added to the matrix, move along
					_, ok := AlleleMap[Location{aln.RName, RefIndex}]

					//if the position is NOT in the matrix, add it
					if !ok {
						AlleleMap[Location{aln.RName, RefIndex}] = &AlleleCount{
							ref[aln.RName][RefIndex], 0, 0, 0, 0, 0, make([]Indel, 0)}
					}

					switch currentSeq[SeqIndex] {
					case dna.A:
						AlleleMap[Location{aln.RName, RefIndex}].BaseA++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.T:
						AlleleMap[Location{aln.RName, RefIndex}].BaseT++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.G:
						AlleleMap[Location{aln.RName, RefIndex}].BaseG++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.C:
						AlleleMap[Location{aln.RName, RefIndex}].BaseC++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					}
					SeqIndex++
					RefIndex++
				}
			} else if aln.Cigar[i].Op != 'H' {
				SeqIndex = SeqIndex + aln.Cigar[i].RunLength
			}
		}
	}

	return AlleleMap
}

// Convert SampleMap to VCF
func AllelesToVcf(input SampleMap) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var current *vcf.Vcf
	var base string
	var RefCount int32
	var i int
	var RefSeq string
	var AltSeq string

	for loc, alleles := range input {

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

		current = &vcf.Vcf{
			Chr:     loc.Chr,
			Pos:     loc.Pos + 1,
			Id:      ".",
			Ref:     base,
			Qual:    1,
			Filter:  ".",
			Info:    ".",
			Format:  "RefCount:AltCount:Cov"}

		// Ref -> A
		current.Alt = "A"
		current.Unknown = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseA, alleles.Counts)
		answer = append(answer, current)

		// Ref -> C
		current.Alt = "C"
		current.Unknown = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseC, alleles.Counts)
		answer = append(answer, current)

		// Ref -> G
		current.Alt = "G"
		current.Unknown = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseG, alleles.Counts)
		answer = append(answer, current)

		// Ref -> T
		current.Alt = "T"
		current.Unknown = fmt.Sprintf("%d:%d:%d", RefCount, alleles.BaseT, alleles.Counts)
		answer = append(answer, current)

		// Ref -> Indel
		for i = 0; i < len(alleles.Indel); i++ {
			RefSeq = dna.BaseToString(alleles.Indel[i].Ref[0])
			AltSeq = dna.BasesToString(alleles.Indel[i].Alt)

			current.Unknown = fmt.Sprintf("%d:%d:%d", RefCount, alleles.Indel[i].Count, alleles.Counts)
			current.Pos = loc.Pos
			current.Ref = RefSeq
			current.Alt = AltSeq
			answer = append(answer, current)
		}
	}
	return answer
}

// Removes positions with insufficient coverage from the map
func FilterAlleles(input SampleMap, coverageThreshold int32) SampleMap {
	for loc, alleles := range input {
		if alleles.Counts < coverageThreshold {
			delete(input, loc)
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
	var AltCount, Counts, Pos int64
	var currentIndel Indel

	file := fileio.EasyOpen(inFilename)
	defer file.Close()

	// Initialize SampleMap
	answer = make(SampleMap)

	var progressMeter int

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {

		if progressMeter%500000 == 0 {
			log.Printf("# Read %d Lines\n", progressMeter)
		}
		progressMeter++

		//order is: Chr, Pos, Id, Ref, Alt, Qual, Filter, Info, Format, Unknown (RefCount:AltCount:Coverage)
		words := strings.Split(line, "\t")

		//exclude header line
		if strings.HasPrefix(words[0], "#") == false {

			// Define each entry to be stored in the map
			var Chr string = words[0]

			// Convert strings to stored values
			RefSeq = dna.StringToBases(words[3])
			AltSeq = dna.StringToBases(words[4])
			UknFmt = strings.Split(words[9], ":")
			AltCount, _ = strconv.ParseInt(UknFmt[1], 10, 32)
			Counts, _ = strconv.ParseInt(UknFmt[2], 10, 32)
			Pos, _ = strconv.ParseInt(words[1], 10, 64)

			// If point mutation
			if len(RefSeq) == 1 && len(AltSeq) == 1 {

				// If the position is in the map move along, else initialize
				// Subtract 1 from Pos for index 0
				_, ok := answer[Location{Chr, Pos-1}]
				if !ok {
					answer[Location{Chr, Pos-1}] = &AlleleCount{
						Ref:    0,
						Counts: 0,
						BaseA:  0,
						BaseC:  0,
						BaseG:  0,
						BaseT:  0,
						Indel:  make([]Indel, 0)}
				}

				answer[Location{Chr, Pos-1}].Ref = RefSeq[0]
				answer[Location{Chr, Pos-1}].Counts = int32(Counts)

				switch AltSeq[0] {
				case dna.A:
					answer[Location{Chr, Pos-1}].BaseA = int32(AltCount)
				case dna.C:
					answer[Location{Chr, Pos-1}].BaseC = int32(AltCount)
				case dna.G:
					answer[Location{Chr, Pos-1}].BaseG = int32(AltCount)
				case dna.T:
					answer[Location{Chr, Pos-1}].BaseT = int32(AltCount)
				}

				// If Indel
			} else if len(RefSeq) != len(AltSeq) {

				// If the position is in the map move along, else initialize
				// VCF stores pos as base prior to indel so subtracting 1 for index 0 is unnecessary
				_, ok := answer[Location{Chr, Pos}]
				if !ok {
					answer[Location{Chr, Pos}] = &AlleleCount{
						Ref:    0,
						Counts: 0,
						BaseA:  0,
						BaseC:  0,
						BaseG:  0,
						BaseT:  0,
						Indel:  make([]Indel, 0)}
				}

				currentIndel = Indel{
					Ref:       RefSeq,
					Alt:       AltSeq,
					Count:     int32(AltCount)}

				answer[Location{Chr, Pos}].Indel = append(answer[Location{Chr, Pos}].Indel, currentIndel)

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
