package alleles

import (
	"fmt"
	"strconv"
	"strings"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"os"
)


type RefAlleleCount struct {
	Chr		string
	Pos 	int32
	Ref 	dna.Base
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
}

type AlleleCount struct {
	Chr		string
	Pos 	int32
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
}





func InitializeAlleleMatrix(refFilename string) (map[string][]*RefAlleleCount, []*fasta.Fasta) {
	//fmt.Printf("Reading reference\n")
	ref := fasta.Read(refFilename)
	fasta.AllToUpper(ref)

	//fmt.Printf("Initializing allele matrix\n")
	alleleMatrix := make(map[string][]*RefAlleleCount, len(ref))
	var i, k int32
	for i = 0; i < int32(len(ref)); i++ {
		alleleMatrix[ref[i].Name] = make([]*RefAlleleCount, len(ref[i].Seq))
		for k = 0; k < int32(len(ref[i].Seq)); k++ {
			alleleMatrix[ref[i].Name][k] = &RefAlleleCount{ref[i].Name, k, ref[i].Seq[k], 0, 0, 0, 0, 0, 0}
		}
	}
	return alleleMatrix, ref
}






func WriteVariantsToFile(input []*AlleleCount, outFilename string) {
	outFile, _ := os.Create(outFilename)
	defer outFile.Close()
	io.WriteString(outFile, "Chr\tPos\tA\tC\tG\tT\tIns\tDel\n")

	var i int32

	for i = 0; i < int32(len(input)); i++ {
		fmt.Fprintf(outFile,
			"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			input[i].Chr,
			input[i].Pos,
			input[i].BaseA,
			input[i].BaseC,
			input[i].BaseG,
			input[i].BaseT,
			input[i].Ins,
			input[i].Del)
	}
}


func WriteVariantsToTerm(input []*AlleleCount) {
	fmt.Printf("Chr\tPos\tA\tC\tG\tT\tIns\tDel\n")

	var i int32

	for i = 0; i < int32(len(input)); i++ {
		fmt.Printf(
			"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			input[i].Chr,
			input[i].Pos,
			input[i].BaseA,
			input[i].BaseC,
			input[i].BaseG,
			input[i].BaseT,
			input[i].Ins,
			input[i].Del)
	}
}






func WholeGenomeAlleleCount(samFilename string, refFilename string) []*RefAlleleCount {
	var answer []*RefAlleleCount
	var current *RefAlleleCount
	var i, k int32

	alleleMatrix, ref := InitializeAlleleMatrix(refFilename)


	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false //var done bool = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln

	//fmt.Printf("Allele matrix initialized. Looping through sam. \n")
	sam.ReadHeader(samFile)

	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		//fmt.Printf("Mark \n")
		if aln.Cigar[0].Op != '*' {
			SeqIndex = 0
			RefIndex = aln.Pos - 1
			for i = 0; i < int32(len(aln.Cigar)); i++ {
				currentSeq = aln.Seq

				//Handle deletion relative to ref
				//Each position deleted is annotated with a Del read
				if aln.Cigar[i].Op == 'D' {
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
						alleleMatrix[aln.RName][RefIndex].Del++
						RefIndex++
					}

				//Handle insertion relative to ref
				//The base after the inserted sequence is annotated with an Ins read
				} else if aln.Cigar[i].Op == 'I' {
					alleleMatrix[aln.RName][RefIndex].Ins++
					SeqIndex = SeqIndex + aln.Cigar[i].RunLength

				//Handle matching pos relative to ref
				} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
						switch currentSeq[SeqIndex] {
						case dna.A:
							alleleMatrix[aln.RName][RefIndex].BaseA++
						case dna.T:
							alleleMatrix[aln.RName][RefIndex].BaseT++
						case dna.G:
							alleleMatrix[aln.RName][RefIndex].BaseG++
						case dna.C:
							alleleMatrix[aln.RName][RefIndex].BaseC++
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

	for i = 0; i < int32(len(ref[i].Name)); i++ {
		for k = 0; k < int32(len(ref[i].Seq)); k++ {
			/*coverage := alleleMatrix[ref[i].Name][k].BaseA +
				alleleMatrix[ref[i].Name][k].BaseC +
				alleleMatrix[ref[i].Name][k].BaseG +
				alleleMatrix[ref[i].Name][k].BaseT +
				alleleMatrix[ref[i].Name][k].Del +
				alleleMatrix[ref[i].Name][k].Ins

			if coverage > 50 {

			 */
				current = &RefAlleleCount{
					ref[i].Name,
					alleleMatrix[ref[i].Name][k].Pos+1,
					alleleMatrix[ref[i].Name][k].Ref,
					alleleMatrix[ref[i].Name][k].BaseA,
					alleleMatrix[ref[i].Name][k].BaseC,
					alleleMatrix[ref[i].Name][k].BaseG,
					alleleMatrix[ref[i].Name][k].BaseT,
					alleleMatrix[ref[i].Name][k].Del,
					alleleMatrix[ref[i].Name][k].Ins}

				answer = append(answer, current)
			//}
			}
		}

	return answer
}






func MaxMinorAllele(allele1 int32, allele2 int32, allele3 int32, allele4 int32, allele5 int32) int32 {
	var minorAllele = allele1
	if allele2 > minorAllele {minorAllele = allele2}
	if allele3 > minorAllele {minorAllele = allele3}
	if allele4 > minorAllele {minorAllele = allele4}
	if allele5 > minorAllele {minorAllele = allele5}
	return minorAllele
}






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





func FindMajorAllele(A int32, C int32, G int32, T int32, Ins int32, Del int32) int32{
	var majorAllele = A
	if C > majorAllele {majorAllele = C}
	if G > majorAllele {majorAllele = G}
	if T > majorAllele {majorAllele = T}
	if Ins > majorAllele {majorAllele = Ins}
	if Del > majorAllele {majorAllele = Del}

	return majorAllele
}






func SortAlleles(input []*AlleleCount) []*AlleleCount {
	var filteredAlleles []*AlleleCount
	var current *AlleleCount
	var i int32
	for i = 0; i < int32(len(input)); i++ {
		//TODO: make coverage an input parameter
		//coverage := input[i].BaseA + input[i].BaseC + input[i].BaseG + input[i].BaseT + input[i].Ins + input[i].Del

		//var minorAllele = FindMinorAllele(input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del)

		//must be a minimum of 50 reads total
		//TODO: make coverage an input parameter
		//if coverage > 50 {
			//minor allele must have greater than 5 supporting reads
			//if minorAllele > 5{
				//minor allele must have a frequency > 1%
				//if float64(minorAllele)/float64(coverage) > 0.01 {
					current = &AlleleCount{input[i].Chr, input[i].Pos, input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del}
					filteredAlleles = append(filteredAlleles, current)
				//}
			//}
		//}
	}
	return filteredAlleles
}






func SortAllelesWithRef(input []*RefAlleleCount) []*AlleleCount {
	var filteredAlleles []*AlleleCount
	var current *AlleleCount
	var i int32
	for i = 0; i < int32(len(input)); i++ {
		var coverage = input[i].BaseA + input[i].BaseC + input[i].BaseG + input[i].BaseT + input[i].Ins + input[i].Del

		var minorAllele = FindMinorAllele(input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del)

		//convert allele count int into a dna.Base for comparisons
		var majorAlleleCount = FindMajorAllele(input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del)
		var majorAllele dna.Base

		switch majorAlleleCount {
		case input[i].BaseA:
			majorAllele = dna.A

		case input[i].BaseC:
			majorAllele = dna.C

		case input[i].BaseG:
			majorAllele = dna.G

		case input[i].BaseT:
			majorAllele = dna.T

		default:
			majorAllele = input[i].Ref
		}

		//if the major alleles is different from the reference allele, include an entry
		if majorAllele != input[i].Ref {
			current = &AlleleCount{input[i].Chr, input[i].Pos, input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del}
			filteredAlleles = append(filteredAlleles, current)
		}

		//if the major allele does match the reference, only include if it passes the following filters
		if majorAllele == input[i].Ref {
			//must be a minimum of 50 reads total
			if coverage > 50 {
				//minor allele must have greater than 5 supporting reads
				if minorAllele > 5{
					//minor allele must have a frequency > 1%
					if float64(minorAllele)/float64(coverage) > 0.01 {
						current = &AlleleCount{input[i].Chr, input[i].Pos, input[i].BaseA, input[i].BaseC, input[i].BaseG, input[i].BaseT, input[i].Ins, input[i].Del}
						filteredAlleles = append(filteredAlleles, current)
					}
				}
			}
		}
	}
	return filteredAlleles
}








func DynamicAlleleCount (samFilename string) []*AlleleCount {
	var i, k int32
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false //var done bool = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln

	sam.ReadHeader(samFile)

	dynamicAlleleMatrix := make(map[string]map[int64]*AlleleCount)

	//var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		//fmt.Println("Alignment no", progressMeter)
		//progressMeter++

		if aln.Cigar[0].Op != '*' {

			//if the chromosome has already been added to the matrix, move along
			_, ok := dynamicAlleleMatrix[aln.RName]

			//if the chromosome is NOT in the matrix, initialize
			if ! ok {
				dynamicAlleleMatrix[aln.RName] = make(map[int64]*AlleleCount)
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
						_, ok := dynamicAlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if ! ok {
							dynamicAlleleMatrix[aln.RName][RefIndex] = &AlleleCount{aln.RName, int32(RefIndex+1), 0, 0, 0, 0, 0, 0}
						}

						dynamicAlleleMatrix[aln.RName][RefIndex].Del++
						RefIndex++
					}

					//Handle insertion relative to ref
					//The base after the inserted sequence is annotated with an Ins read
				} else if aln.Cigar[i].Op == 'I' {

					//if the position has already been added to the matrix, move along
					_, ok := dynamicAlleleMatrix[aln.RName][RefIndex]

					//if the position is NOT in the matrix, add it
					if ! ok {
						dynamicAlleleMatrix[aln.RName][RefIndex] = &AlleleCount{aln.RName, int32(RefIndex+1), 0, 0, 0, 0, 0, 0}
					}

					dynamicAlleleMatrix[aln.RName][RefIndex].Ins++
					SeqIndex = SeqIndex + aln.Cigar[i].RunLength

					//Handle matching pos relative to ref
				} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {
					for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

						//if the position has already been added to the matrix, move along
						_, ok := dynamicAlleleMatrix[aln.RName][RefIndex]

						//if the position is NOT in the matrix, add it
						if ! ok {
							dynamicAlleleMatrix[aln.RName][RefIndex] = &AlleleCount{aln.RName, int32(RefIndex+1), 0, 0, 0, 0, 0, 0}
						}

						switch currentSeq[SeqIndex] {
						case dna.A:
							dynamicAlleleMatrix[aln.RName][RefIndex].BaseA++
						case dna.T:
							dynamicAlleleMatrix[aln.RName][RefIndex].BaseT++
						case dna.G:
							dynamicAlleleMatrix[aln.RName][RefIndex].BaseG++
						case dna.C:
							dynamicAlleleMatrix[aln.RName][RefIndex].BaseC++
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
	//extract values from the map and store as slice for answer

	var alleles []*AlleleCount
	for _, chr := range dynamicAlleleMatrix {
		for _, allele := range chr {
			alleles = append(alleles, allele)
		}
	}
	return alleles
}







func ExtractAllelesWholeGenome(inFilename string, refFilename string) []*AlleleCount {
	alleles := WholeGenomeAlleleCount(inFilename, refFilename)
	//return alleles
	filteredAlleles := SortAllelesWithRef(alleles)
	return filteredAlleles
	//WriteVariantsToFile(filteredAlleles, outFilename)
}

func ExtractAlleles(samFilename string) []*AlleleCount {
	alleles := DynamicAlleleCount(samFilename)
	//return alleles
	filteredAlleles := SortAlleles(alleles)
	return filteredAlleles
	//WriteVariantsToFile(filteredAlleles, outFilename)
}






func ReadAlleleCounts(inFilename string) []*AlleleCount {
	var line string
	var answer []*AlleleCount
	var current *AlleleCount
	var doneReading = false

	file := fileio.EasyOpen(inFilename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		//order is: Chr, Pos, A, C, G, T, Ins, Del
		words := strings.Split(line, "\t")

		//exclude header line
		if words[1] != "Pos"{
			var Chr string = words[0]
			var Pos, _ = strconv.ParseInt(words[1], 10, 32)
			var BaseA, _ = strconv.ParseInt(words[2], 10, 32)
			var BaseC, _ = strconv.ParseInt(words[3], 10, 32)
			var BaseG, _ = strconv.ParseInt(words[4], 10, 32)
			var BaseT, _ = strconv.ParseInt(words[5], 10, 32)
			var Ins, _ = strconv.ParseInt(words[6], 10, 32)
			var Del, _ = strconv.ParseInt(words[7], 10, 32)

			current = &AlleleCount{
				Chr:	Chr,
				Pos:	int32(Pos),
				BaseA:	int32(BaseA),
				BaseC:	int32(BaseC),
				BaseG:	int32(BaseG),
				BaseT:	int32(BaseT),
				Ins:	int32(Ins),
				Del:	int32(Del)}

			answer = append(answer, current)
		}

	}
	return answer
}



/*
func main() {
	inFile := os.Args[1]
	outFile := os.Args[2]
	alleles := DynamicAlleleCount(inFile)
	filteredAlleles := SortAlleles(alleles)
	WriteVariantsToFile(filteredAlleles, outFile)
}
 */