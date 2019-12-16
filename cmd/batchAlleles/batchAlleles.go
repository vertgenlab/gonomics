package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"io"
	"log"
	"os"
)


func usage() {
	fmt.Print(
		"callVariants - Inputs a directory of allele count files and ouputs a normalized file with Z scores for each alternative allele.\n" +
			"Usage:\n" +
			" callVariants [options] inputDirectory/ \n" +
			"options:\n")
	flag.PrintDefaults()
}


func WriteZscoreToFile(input map[string]map[int32][]*alleles.BatchZScores, outFilename string) {
	outFile, _ := os.Create(outFilename)
	defer outFile.Close()
	io.WriteString(outFile, "Sample\tChr\tPos\tCoverage\tA%\tC%\tG%\tT%\tIns%\tDel%\tZscoreA\tZscoreC\tZscoreG\tZscoreT\tZscoreIns\tZscoreDel\n")

	var i int

	for _, pos := range input {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {
				fmt.Fprintf(outFile,
					"%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					data[i].Sample,
					data[i].Chr,
					data[i].Pos,
					data[i].Coverage,
					data[i].NormBaseA,
					data[i].NormBaseC,
					data[i].NormBaseG,
					data[i].NormBaseT,
					data[i].NormIns,
					data[i].NormDel,
					data[i].ZscoreA,
					data[i].ZscoreC,
					data[i].ZscoreG,
					data[i].ZscoreT,
					data[i].ZscoreIns,
					data[i].ZscoreDel)
			}
		}
	}
}


func WriteSortedAllelesToFile(input map[string]map[int32][]*alleles.MultiplexSortedAlleles, outFilename string) {
	outFile, _ := os.Create(outFilename)
	defer outFile.Close()
	io.WriteString(outFile, "Sample\tChr\tPos\tA\tC\tG\tT\tIns\tDel\tBkgdA\tBkgdC\tBkgdG\tBkgdT\tBkgdIns\tBkgdDel\n")

	var i int

	for _, pos := range input {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {
				fmt.Fprintf(outFile,
					"%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					data[i].Sample,
					data[i].Chr,
					data[i].Pos,
					data[i].BaseA,
					data[i].BaseC,
					data[i].BaseG,
					data[i].BaseT,
					data[i].Ins,
					data[i].Del,
					data[i].BkgdA,
					data[i].BkgdC,
					data[i].BkgdG,
					data[i].BkgdT,
					data[i].BkgdIns,
					data[i].BkgdDel)
			}
		}
	}
}


func WriteSortedAllelesToTerm(input map[string]map[int32][]*alleles.MultiplexSortedAlleles) {

	fmt.Printf("Sample\tChr\tPos\tA\tC\tG\tT\tIns\tDel\tBkgdA\tBkgdC\tBkgdG\tBkgdT\tBkgdIns\tBkgdDel\n")

	var i int

	for _, pos := range input {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {
				fmt.Printf(
					"%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					data[i].Sample,
					data[i].Chr,
					data[i].Pos,
					data[i].BaseA,
					data[i].BaseC,
					data[i].BaseG,
					data[i].BaseT,
					data[i].Ins,
					data[i].Del,
					data[i].BkgdA,
					data[i].BkgdC,
					data[i].BkgdG,
					data[i].BkgdT,
					data[i].BkgdIns,
					data[i].BkgdDel)
			}
		}
	}
}

/*
func WriteFEToFile(input map[string]map[int32][]*batch.BatchFishersE, outFilename string) {
	outFile, _ := os.Create(outFilename)
	defer outFile.Close()
	io.WriteString(outFile, "Sample\tChr\tPos\tA\tC\tG\tT\tIns\tDel\tFishers Exact A\tFishers Exact C\tFishers Exact G\tFishers Exact T\tFishers Exact Ins\tFishers ExactDel\n")

	var i int

	for _, pos := range input {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {
				fmt.Fprintf(outFile,
					"%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n",
					data[i].Sample,
					data[i].Chr,
					data[i].Pos,
					data[i].BaseA,
					data[i].BaseC,
					data[i].BaseG,
					data[i].BaseT,
					data[i].Ins,
					data[i].Del,
					data[i].FetA,
					data[i].FetC,
					data[i].FetG,
					data[i].FetT,
					data[i].FetIns,
					data[i].FetDel)
			}
		}
	}
}
 */


func main() {
	var expectedNumArgs int=1
	var outFile *string = flag.String("o", "stdout", "Write output to a file")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	flag.Parse()
	inDirectory := flag.Arg(0)

	fmt.Printf("Merging Samples\n")
	SampleMap := alleles.CreateSampleMap(inDirectory)
	//fmt.Printf("Normalizing Coverage\n")
	//NormalizedSampleMap := batch.NormalizeCoverage(SampleMap)
	//fmt.Printf("Calculating Z Scores\n")
	//Zscores := batch.CalculateZscores(NormalizedSampleMap)
	//fmt.Printf("Filtering Variants\n")
	//output := batch.FilterAF(Zscores)
	//WriteVariantsToFile(Zscores, "test")

	fmt.Printf("Sorting Alleles\n")
	output := alleles.BatchSortAlleles(SampleMap)

	if *outFile == "stdout" {
		WriteSortedAllelesToTerm(output)
	} else {
		fmt.Printf("Writing to File\n")
		WriteSortedAllelesToFile(output, *outFile)
	}
}