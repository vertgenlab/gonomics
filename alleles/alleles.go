package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strconv"
	"strings"
)

// TODO: in vcf format multiple possible alleles into a single line
type AlleleCount struct {
	Ref    dna.Base
	Pos	   int64
	Counts int32
	// Base counts are stored as slice with [0] = Total Count [1] = Forward Reads and [2] = Reverse Reads
	// Total Count is added for processing unpaired sequencing
	BaseAF int32
	BaseCF int32
	BaseGF int32
	BaseTF int32
	BaseAR int32
	BaseCR int32
	BaseGR int32
	BaseTR int32
	Indel []Indel
}

type Indel struct {
	Ref   	[]dna.Base
	Alt 	[]dna.Base
	CountF 	int32
	CountR	int32
}

type Location struct {
	Chr string
	Pos int64
}

// Map structure: map[Chromosome]map[Position]*AlleleCount
type SampleMap map[Location]*AlleleCount


// Convert SampleMap to VCF
func AllelesToVcf(input SampleMap) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var current *vcf.Vcf
	var base string
	var RefCount, RefCountF, RefCountR int32
	var i int

	var progressMeter int
	for loc, alleles := range input {
		progressMeter++
		if progressMeter%5000000 == 0 {
			log.Printf("processed %d positions", progressMeter)
		}
		switch alleles.Ref {
		case dna.A:
			base = "A"
			RefCount = alleles.BaseAF + alleles.BaseAR
			RefCountF = alleles.BaseAF
			RefCountR = alleles.BaseAR
		case dna.C:
			base = "C"
			RefCount = alleles.BaseCF + alleles.BaseCR
			RefCountF = alleles.BaseCF
			RefCountR = alleles.BaseCR
		case dna.G:
			base = "G"
			RefCount = alleles.BaseGF + alleles.BaseGR
			RefCountF = alleles.BaseGF
			RefCountR = alleles.BaseGR
		case dna.T:
			base = "T"
			RefCount = alleles.BaseTF + alleles.BaseGR
			RefCountF = alleles.BaseTF
			RefCountR = alleles.BaseTR
		default:
			base = "NA"
			RefCount = 0
			RefCountF = 0
			RefCountR = 0
		}

		// Ref -> A
		Sa := make([]string, 1)
		Sa[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseAF + alleles.BaseAR, alleles.BaseAF, alleles.BaseAR, alleles.Counts)

		current = &vcf.Vcf{
			Chr:    loc.Chr,
			Pos:    loc.Pos + 1,
			Id:     ".",
			Ref:    base,
			Alt:    "A",
			Qual:   1,
			Filter: ".",
			Info:   ".",
			Format: "RefCount,For,Rev:AltCount,For,Rev:Cov",
			Sample: Sa}

		answer = append(answer, current)

		// Ref -> C
		Sc := make([]string, 1)
		Sc[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseCF + alleles.BaseCR, alleles.BaseCF, alleles.BaseCR, alleles.Counts)

		current = &vcf.Vcf{
			Chr:    loc.Chr,
			Pos:    loc.Pos + 1,
			Id:     ".",
			Ref:    base,
			Alt:    "C",
			Qual:   1,
			Filter: ".",
			Info:   ".",
			Format: "RefCount,For,Rev:AltCount,For,Rev:Cov",
			Sample: Sc}

		answer = append(answer, current)

		// Ref -> G
		Sg := make([]string, 1)
		Sg[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseGF + alleles.BaseGR, alleles.BaseGF, alleles.BaseGR, alleles.Counts)

		current = &vcf.Vcf{
			Chr:    loc.Chr,
			Pos:    loc.Pos + 1,
			Id:     ".",
			Ref:    base,
			Alt:    "G",
			Qual:   1,
			Filter: ".",
			Info:   ".",
			Format: "RefCount,For,Rev:AltCount,For,Rev:Cov",
			Sample: Sg}

		answer = append(answer, current)

		// Ref -> T
		St := make([]string, 1)
		St[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseTF + alleles.BaseTR, alleles.BaseTF, alleles.BaseTR, alleles.Counts)

		current = &vcf.Vcf{
			Chr:    loc.Chr,
			Pos:    loc.Pos + 1,
			Id:     ".",
			Ref:    base,
			Alt:    "T",
			Qual:   1,
			Filter: ".",
			Info:   ".",
			Format: "RefCount,For,Rev:AltCount,For,Rev:Cov",
			Sample: St}

		answer = append(answer, current)

		// Ref -> Indel
		var RefSeq []string = make([]string, len(alleles.Indel))
		var AltSeq []string = make([]string, len(alleles.Indel))
		var Sindels [][]string = make([][]string, len(alleles.Indel))

		for i = 0; i < len(alleles.Indel); i++ {
			RefSeq[i] = dna.BasesToString(alleles.Indel[i].Ref)
			AltSeq[i] = dna.BasesToString(alleles.Indel[i].Alt)

			Sindel := make([]string, 1)
			Sindel[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.Indel[i].CountF + alleles.Indel[i].CountR, alleles.Indel[i].CountF, alleles.Indel[i].CountR, alleles.Counts)
			Sindels[i] = Sindel

			current = &vcf.Vcf{
				Chr:    loc.Chr,
				Pos:    loc.Pos,
				Id:     ".",
				Ref:    RefSeq[i],
				Alt:    AltSeq[i],
				Qual:   1,
				Filter: ".",
				Info:   ".",
				Format: "RefCount,For,Rev:AltCount,For,Rev:Cov",
				Sample: Sindels[i]}

			answer = append(answer, current)
		}
		delete(input, loc)
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
	var UknFmt, AltCountString []string
	var AltCountF, AltCountR, Counts, Pos int64
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
			AltCountString = strings.Split(UknFmt[1], ",")
			AltCountF, _ = strconv.ParseInt(AltCountString[1], 10, 32)
			AltCountR, _ = strconv.ParseInt(AltCountString[2], 10, 32)
			Counts, _ = strconv.ParseInt(UknFmt[2], 10, 32)
			Pos, _ = strconv.ParseInt(words[1], 10, 64)

			// If point mutation
			if len(RefSeq) == 1 && len(AltSeq) == 1 {

				// If the position is in the map move along, else initialize
				// Subtract 1 from Pos for index 0
				_, ok := answer[Location{Chr, Pos - 1}]
				if !ok {
					answer[Location{Chr, Pos - 1}] = &AlleleCount{
						Ref:    0,
						Counts: 0,
						BaseAF:  0,
						BaseCF:  0,
						BaseGF:  0,
						BaseTF:  0,
						BaseAR:  0,
						BaseCR:  0,
						BaseGR:  0,
						BaseTR:  0,
						Indel:  make([]Indel, 0)}
				}

				answer[Location{Chr, Pos - 1}].Ref = RefSeq[0]
				answer[Location{Chr, Pos - 1}].Counts = int32(Counts)

				switch AltSeq[0] {
				case dna.A:
					answer[Location{Chr, Pos - 1}].BaseAF = int32(AltCountF)
					answer[Location{Chr, Pos - 1}].BaseAR = int32(AltCountR)
				case dna.C:
					answer[Location{Chr, Pos - 1}].BaseCF = int32(AltCountF)
					answer[Location{Chr, Pos - 1}].BaseCR = int32(AltCountR)
				case dna.G:
					answer[Location{Chr, Pos - 1}].BaseGF = int32(AltCountF)
					answer[Location{Chr, Pos - 1}].BaseGR = int32(AltCountR)
				case dna.T:
					answer[Location{Chr, Pos - 1}].BaseTF = int32(AltCountF)
					answer[Location{Chr, Pos - 1}].BaseTR = int32(AltCountR)
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
						BaseAF:  0,
						BaseCF:  0,
						BaseGF:  0,
						BaseTF:  0,
						BaseAR:  0,
						BaseCR:  0,
						BaseGR:  0,
						BaseTR:  0,
						Indel:  make([]Indel, 0)}
				}

				currentIndel = Indel{
					Ref:   RefSeq,
					Alt:   AltSeq,
					CountF: int32(AltCountF),
					CountR: int32(AltCountR)}

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