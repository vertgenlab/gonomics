package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"io"
	"io/ioutil"
	"os"
)

type BatchAlleleCount struct {
	Sample	string
	Ref		dna.Base
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
}

// Map structure: map[Chromosome]map[Position][Sample]*BatchAlleleCount
type BatchSampleMap map[string]map[int64][]*BatchAlleleCount



// Input a directory filled with AlleleCount (.ac) files. Merges the files into a nested map (chr->pos->[sample]data)
func CreateSampleMap(inDirectory string) BatchSampleMap {

	files, _ := ioutil.ReadDir(inDirectory)
	var SampleName string
	var current *BatchAlleleCount
	var fileCount int = 1

	SampleMap := make(map[string]map[int64][]*BatchAlleleCount)

	for _, file := range files {

		SampleName = file.Name()
		SamplePath := fmt.Sprintf("%s%s", inDirectory, SampleName)
		fmt.Printf("#\n# File %d\n", fileCount)
		AlleleCounts := ReadAlleleCounts(SamplePath)
		fileCount++

		// Loop through input map and deposit info into output map
		for chrName, chr := range AlleleCounts {
			for pos, alleles := range chr {

				//if the chromosome has already been added to the matrix, move along
				_, ok := SampleMap[chrName]

				//if the chromosome is NOT in the matrix, initialize
				if ! ok {
					SampleMap[chrName] = make(map[int64][]*BatchAlleleCount)
				}

				// Check if position is in the map, if not then add the background struct
				_, okk := SampleMap[chrName][pos]
				if ! okk {
					current = &BatchAlleleCount{"Background", alleles.Ref, 0, 0, 0, 0, 0, 0}
					SampleMap[chrName][pos] = append(SampleMap[chrName][pos], current)
				}

				//add in and increment background as entry zero in batch allele count slice
				SampleMap[chrName][pos][0].BaseA = SampleMap[chrName][pos][0].BaseA + alleles.BaseA
				SampleMap[chrName][pos][0].BaseC = SampleMap[chrName][pos][0].BaseC + alleles.BaseC
				SampleMap[chrName][pos][0].BaseG = SampleMap[chrName][pos][0].BaseG + alleles.BaseG
				SampleMap[chrName][pos][0].BaseT = SampleMap[chrName][pos][0].BaseT + alleles.BaseT
				SampleMap[chrName][pos][0].Ins = SampleMap[chrName][pos][0].Ins + alleles.Ins
				SampleMap[chrName][pos][0].Del = SampleMap[chrName][pos][0].Del + alleles.Del




				current = &BatchAlleleCount{
					Sample: SampleName,
					Ref:	alleles.Ref,
					BaseA:  alleles.BaseA,
					BaseC:  alleles.BaseC,
					BaseG:  alleles.BaseG,
					BaseT:  alleles.BaseT,
					Ins:    alleles.Ins,
					Del:    alleles.Del}

				SampleMap[chrName][pos] = append(SampleMap[chrName][pos], current)
			}
		}
	}

	return SampleMap
}


func WriteBatchAlleleCounts(input BatchSampleMap, output string) {

	var outFile *os.File

	if output != "stdout" {
		fmt.Printf("#Creating Output File\n")
		outFile, _ = os.Create(output)
		defer outFile.Close()
		io.WriteString(outFile, "Sample\tChr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\n")
	} else {
		fmt.Printf("#No Output Specified. Printing to STDOUT.\n")
		fmt.Printf("Sample\tChr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\n")
	}

	var base string
	var progressMeter int
	var i int

	for chrName, chr := range input {
		for pos, alleles := range chr {
			//TODO: change i = 0 to i = 1 after troubleshooting so that the write function ignores the background struct
			for i = 0; i < len(alleles); i++ {

			switch alleles[i].Ref {
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
					"%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
					alleles[i].Sample,
					chrName,
					pos + 1,
					base,
					alleles[i].BaseA,
					alleles[i].BaseC,
					alleles[i].BaseG,
					alleles[i].BaseT,
					alleles[i].Ins,
					alleles[i].Del)
			} else {
				// Write to stdout
				fmt.Printf(
					"%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
					alleles[i].Sample,
					chrName,
					pos+1,
					base,
					alleles[i].BaseA,
					alleles[i].BaseC,
					alleles[i].BaseG,
					alleles[i].BaseT,
					alleles[i].Ins,
					alleles[i].Del)
				}
			}
		}
	}
}