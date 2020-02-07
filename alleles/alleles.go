package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/routines"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strconv"
	"strings"
	"sync"
)

// TODO: in vcf format multiple possible alleles into a single line
type AlleleCount struct {
	Ref    dna.Base
	Counts int32
	// Base counts are stored as slice with [0] = Total Count [1] = Forward Reads and [2] = Reverse Reads
	// Total Count is added for processing unpaired sequencing
	// TODO: instead of slice of ints, maybe make a second map to store the auxillary information that could be toggled with user inputs
	BaseA []int32
	BaseC []int32
	BaseG []int32
	BaseT []int32
	Indel []Indel
}

type Indel struct {
	Ref   []dna.Base
	Alt   []dna.Base
	Count []int32
}

type Location struct {
	Chr string
	Pos int64
}

// Map structure: map[Chromosome]map[Position]*AlleleCount
type SampleMap map[Location]*AlleleCount

type RefMap map[string][]dna.Base

// Temporary function until this is added to the fasta package
func refToMap(refFilename string) RefMap {
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

		if progressMeter%10000 == 0 {
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
							Ref: ref[aln.RName][RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					// Keep track of deleted sequence
					indelSeq = append(indelSeq, ref[aln.RName][RefIndex])

					AlleleMap[Location{aln.RName, RefIndex}].Counts++
					RefIndex++
				}

				Match = false
				for j = 0; j < len(AlleleMap[Location{aln.RName, OrigRefIndex}].Indel); j++ {
					// If the deletion has already been seen before, increment the existing entry
					// For a deletion the indelSeq should match the Ref
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Ref) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Alt) == 0 {
						AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, OrigRefIndex}].Indel[j].Count[2]++
						}
						Match = true
						break
					}
				}

				// If the deletion has not been seen before, then append it to the Del slice
				// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
				if Match == false {

					currentIndel = Indel{indelSeq, indelSeq[:1], make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == false {
						currentIndel.Count[2]++
					}
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
						Ref: ref[aln.RName][RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
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
					// If the inserted sequence matches a previously inserted sequence, then increment the count
					// For an insertion, the indelSeq should match the Alt
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Ref) == 0 {
						AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].Indel[j].Count[2]++
						}
						Match = true
						break
					}
				}

				if Match == false {
					currentIndel = Indel{indelSeq[:1], indelSeq, make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == true {
						currentIndel.Count[2]++
					}
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
							Ref: ref[aln.RName][RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					switch currentSeq[SeqIndex] {
					case dna.A:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseA[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseA[2]++
						}
						AlleleMap[Location{aln.RName, RefIndex}].BaseA[0]++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.T:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseT[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseT[2]++
						}
						AlleleMap[Location{aln.RName, RefIndex}].BaseT[0]++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.G:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseG[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseG[2]++
						}
						AlleleMap[Location{aln.RName, RefIndex}].BaseG[0]++
						AlleleMap[Location{aln.RName, RefIndex}].Counts++
					case dna.C:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseC[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[Location{aln.RName, RefIndex}].BaseC[2]++
						}
						AlleleMap[Location{aln.RName, RefIndex}].BaseC[0]++
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
			RefCount = alleles.BaseA[0]
			RefCountF = alleles.BaseA[1]
			RefCountR = alleles.BaseA[2]
		case dna.C:
			base = "C"
			RefCount = alleles.BaseC[0]
			RefCountF = alleles.BaseC[1]
			RefCountR = alleles.BaseC[2]
		case dna.G:
			base = "G"
			RefCount = alleles.BaseG[0]
			RefCountF = alleles.BaseG[1]
			RefCountR = alleles.BaseG[2]
		case dna.T:
			base = "T"
			RefCount = alleles.BaseT[0]
			RefCountF = alleles.BaseT[1]
			RefCountR = alleles.BaseT[2]
		default:
			base = "NA"
			RefCount = 0
			RefCountF = 0
			RefCountR = 0
		}

		// Ref -> A
		Sa := make([]string, 1)
		Sa[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseA[0], alleles.BaseA[1], alleles.BaseA[2], alleles.Counts)

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
		Sc[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseC[0], alleles.BaseC[1], alleles.BaseC[2], alleles.Counts)

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
		Sg[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseG[0], alleles.BaseG[1], alleles.BaseG[2], alleles.Counts)

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
		St[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.BaseT[0], alleles.BaseT[1], alleles.BaseT[2], alleles.Counts)

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
			Sindel[0] = fmt.Sprintf("%d,%d,%d:%d,%d,%d:%d", RefCount, RefCountF, RefCountR, alleles.Indel[i].Count[0], alleles.Indel[i].Count[1], alleles.Indel[i].Count[2], alleles.Counts)
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
	var AltCount, AltCountF, AltCountR, Counts, Pos int64
	var currentIndel Indel
	var Alt []int32

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
			AltCount, _ = strconv.ParseInt(AltCountString[0], 10, 32)
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
						BaseA:  make([]int32, 3),
						BaseC:  make([]int32, 3),
						BaseG:  make([]int32, 3),
						BaseT:  make([]int32, 3),
						Indel:  make([]Indel, 0)}
				}

				answer[Location{Chr, Pos - 1}].Ref = RefSeq[0]
				answer[Location{Chr, Pos - 1}].Counts = int32(Counts)

				Alt = make([]int32, 3)
				Alt[0] = int32(AltCount)
				Alt[1] = int32(AltCountF)
				Alt[2] = int32(AltCountR)

				switch AltSeq[0] {
				case dna.A:
					answer[Location{Chr, Pos - 1}].BaseA = Alt
				case dna.C:
					answer[Location{Chr, Pos - 1}].BaseC = Alt
				case dna.G:
					answer[Location{Chr, Pos - 1}].BaseG = Alt
				case dna.T:
					answer[Location{Chr, Pos - 1}].BaseT = Alt
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
						BaseA:  make([]int32, 3),
						BaseC:  make([]int32, 3),
						BaseG:  make([]int32, 3),
						BaseT:  make([]int32, 3),
						Indel:  make([]Indel, 0)}
				}

				currentIndel = Indel{
					Ref:   RefSeq,
					Alt:   AltSeq,
					Count: Alt}

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

// Count Alleles with Go Routines
type GoCountInput struct {
	matrix *SafeMap
	Ref    RefMap
	MapQT  int64
}

type SafeMap struct {
	m   SampleMap
	mux sync.Mutex
}

func getPos(matrix *SafeMap, key Location, ref RefMap) {
	//matrix.mux.Lock()
	_, ok := matrix.m[key]
	if !ok {
		current := &AlleleCount{
			Ref:    ref[key.Chr][key.Pos],
			Counts: 0,
			BaseA:  make([]int32, 3),
			BaseC:  make([]int32, 3),
			BaseG:  make([]int32, 3),
			BaseT:  make([]int32, 3),
			Indel:  make([]Indel, 0)}

		matrix.m[key] = current
	}
	//matrix.mux.Unlock()
}

func writeToMap(input chan SampleMap, output chan SampleMap, wg *sync.WaitGroup) {
	//start := time.Now()
	defer wg.Done()
	writeMap := make(SampleMap)

	for readMap := range input {

		for key, value := range readMap {
			_, ok := writeMap[key]
			if !ok {
				writeMap[key] = value
			} else {
				writeMap[key].Counts = writeMap[key].Counts + readMap[key].Counts
				writeMap[key].BaseA = addSlice(writeMap[key].BaseA, readMap[key].BaseA)
				writeMap[key].BaseC = addSlice(writeMap[key].BaseC, readMap[key].BaseC)
				writeMap[key].BaseG = addSlice(writeMap[key].BaseG, readMap[key].BaseG)
				writeMap[key].BaseT = addSlice(writeMap[key].BaseT, readMap[key].BaseT)
				addIndels(writeMap[key], readMap[key])
			}
		}

	}
	output <- writeMap
	//fmt.Println("Write took", time.Since(start))
}

func mergeMaps(a SampleMap, b SampleMap) SampleMap {
	for key, value := range b {
		_, ok := a[key]
		if !ok {
			a[key] = value
		} else {
			a[key].Counts = a[key].Counts + b[key].Counts
			a[key].BaseA = addSlice(a[key].BaseA, b[key].BaseA)
			a[key].BaseC = addSlice(a[key].BaseC, b[key].BaseC)
			a[key].BaseG = addSlice(a[key].BaseG, b[key].BaseG)
			a[key].BaseT = addSlice(a[key].BaseT, b[key].BaseT)
			addIndels(a[key], b[key])
		}
	}
	return a
}

func addIndels(write *AlleleCount, read *AlleleCount) {
	var Match bool
	var i, j int

	Match = false
	for i = 0; i < len(read.Indel); i++ {
		for j = 0; j < len(write.Indel); j++ {
			// If the deletion has already been seen before, increment the existing entry
			// For a deletion the indelSeq should match the Ref
			if dna.CompareSeqsIgnoreCase(read.Indel[i].Ref, write.Indel[j].Ref) == 0 &&
				dna.CompareSeqsIgnoreCase(read.Indel[i].Alt, write.Indel[j].Alt) == 0 {
				write.Indel[j].Count = addSlice(write.Indel[j].Count, read.Indel[i].Count)
				Match = true
				break
			}
		}

		if Match == false {
			write.Indel = append(write.Indel, read.Indel[i])
		}
	}
}

func addSlice(a []int32, b []int32) []int32 {
	c := make([]int32, len(a))
	for i := 0; i < len(a); i++ {
		c[i] = a[i] + b[i]
	}
	return c
}

func goCountAlleles(samRecord interface{}, input interface{}) interface{} {
	aln := samRecord.(*sam.SamAln)
	data := input.(*GoCountInput)
	ref := data.Ref
	minMapQ := data.MapQT

	answer := make(SampleMap)

	var i, k int32
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64

	// If read is unmapped then go to the next alignment
	if aln.Cigar[0].Op == '*' {
		return nil
	}

	// If mapping quality is less than the threshold then go to next alignment
	if aln.MapQ < minMapQ {
		return nil
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

				loc := Location{aln.RName, RefIndex}
				answer[loc] = &AlleleCount{
					Ref:    ref[aln.RName][RefIndex],
					Counts: 1,
					BaseA:  make([]int32, 3),
					BaseC:  make([]int32, 3),
					BaseG:  make([]int32, 3),
					BaseT:  make([]int32, 3),
					Indel:  make([]Indel, 0)}

				// Keep track of deleted sequence
				indelSeq = append(indelSeq, ref[aln.RName][RefIndex])
				RefIndex++
			}

			loc := Location{aln.RName, OrigRefIndex}

			currentIndel = Indel{indelSeq, indelSeq[:1], make([]int32, 3)}
			currentIndel.Count[0]++
			if sam.IsForwardRead(aln) == true {
				currentIndel.Count[1]++
			} else if sam.IsReverseRead(aln) == false {
				currentIndel.Count[2]++
			}

			answer[loc].Indel = append(answer[loc].Indel, currentIndel)

			//Handle insertion relative to ref
			//The base after the inserted sequence is annotated with an Ins read
		} else if aln.Cigar[i].Op == 'I' {

			loc := Location{aln.RName, RefIndex}
			answer[loc] = &AlleleCount{
				Ref:    ref[aln.RName][RefIndex],
				Counts: 0,
				BaseA:  make([]int32, 3),
				BaseC:  make([]int32, 3),
				BaseG:  make([]int32, 3),
				BaseT:  make([]int32, 3),
				Indel:  make([]Indel, 0)}

			// Loop through read sequence and keep track of the inserted bases
			indelSeq = make([]dna.Base, 1)

			// First base in indel is the base prior to the indel sequence per VCF standard format
			indelSeq[0] = ref[aln.RName][RefIndex-1]

			// Keep track of inserted sequence by moving along the read
			for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
				indelSeq = append(indelSeq, currentSeq[SeqIndex])
				SeqIndex++
			}

			currentIndel = Indel{indelSeq[:1], indelSeq, make([]int32, 3)}
			currentIndel.Count[0]++
			if sam.IsForwardRead(aln) == true {
				currentIndel.Count[1]++
			} else if sam.IsReverseRead(aln) == true {
				currentIndel.Count[2]++
			}

			answer[loc].Indel = append(answer[loc].Indel, currentIndel)

			// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

			//Handle matching pos relative to ref
		} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

			for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

				loc := Location{aln.RName, RefIndex}
				answer[loc] = &AlleleCount{
					Ref:    ref[aln.RName][RefIndex],
					Counts: 0,
					BaseA:  make([]int32, 3),
					BaseC:  make([]int32, 3),
					BaseG:  make([]int32, 3),
					BaseT:  make([]int32, 3),
					Indel:  make([]Indel, 0)}

				switch currentSeq[SeqIndex] {
				case dna.A:
					if sam.IsForwardRead(aln) == true {
						answer[loc].BaseA[1]++
					} else if sam.IsReverseRead(aln) == true {
						answer[loc].BaseA[2]++
					}
					answer[loc].BaseA[0]++
					answer[loc].Counts++
				case dna.T:
					if sam.IsForwardRead(aln) == true {
						answer[loc].BaseT[1]++
					} else if sam.IsReverseRead(aln) == true {
						answer[loc].BaseT[2]++
					}
					answer[loc].BaseT[0]++
					answer[loc].Counts++
				case dna.G:
					if sam.IsForwardRead(aln) == true {
						answer[loc].BaseG[1]++
					} else if sam.IsReverseRead(aln) == true {
						answer[loc].BaseG[2]++
					}
					answer[loc].BaseG[0]++
					answer[loc].Counts++
				case dna.C:
					if sam.IsForwardRead(aln) == true {
						answer[loc].BaseC[1]++
					} else if sam.IsReverseRead(aln) == true {
						answer[loc].BaseC[2]++
					}
					answer[loc].BaseC[0]++
					answer[loc].Counts++
				}

				SeqIndex++
				RefIndex++
			}
		} else if aln.Cigar[i].Op != 'H' {
			SeqIndex = SeqIndex + aln.Cigar[i].RunLength
		}
	}
	return answer
}

func GoCountAlleles(refFilename string, samFilename string, minMapQ int64, threads int) SampleMap {

	// Read in reference
	fmt.Printf("#Reading Reference\n")
	ref := refToMap(refFilename)

	// Initialize empty map of chromosomes to map of positions
	AlleleMap := make(SampleMap)

	input := &GoCountInput{&SafeMap{m: AlleleMap}, ref, minMapQ}

	channel, _ := routines.GoWorkOnLine(samFilename, input, goCountAlleles, routines.NextSamLine, threads)

	interMaps := 20

	sendMap := make(chan SampleMap)
	recieveMap := make(chan SampleMap, interMaps)

	var wg sync.WaitGroup
	for i := 0; i < interMaps; i++ {
		wg.Add(1)
		go writeToMap(sendMap, recieveMap, &wg)
	}

	go func() {
		wg.Wait()
		close(recieveMap)
	}()

	for j := range channel {
		if j != nil {
			sendMap <- j.(SampleMap)
		}
	}
	close(sendMap)

	wg.Wait()
	for k := range recieveMap {
		AlleleMap = mergeMaps(AlleleMap, k)
	}

	return AlleleMap
}

func safeMapToMap(mSafe sync.Map) SampleMap {
	var answer SampleMap
	answer = make(SampleMap)
	mSafe.Range(func(key interface{}, value interface{}) bool {
		answer[key.(Location)] = value.(*AlleleCount)
		return true
	})
	return answer
}
