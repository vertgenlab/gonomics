package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

type ScoreSettings struct {
	ScoreType  string
	InFileList string
	OutFile    string
}

func scoreUsage(scoreFlags *flag.FlagSet) {
	fmt.Print("samAssembler score - Validate assembly accuracy from five-way alignment including known divergent sequences.\n" +
		"The user must specify a list of filenames, corresponding to each multiFa to be used in joint scoring (usually one file per chromosome).\n" +
		"'scoreType' may be chosen from the following strings: baseMatrix, baseMatrixByRefBase\n" +
		"Usage:\n" +
		"samAssembler score scoreType inFileList outFile.txt\n" +
		"Options:\n")
	scoreFlags.PrintDefaults()
}

func parseScoreArgs() {
	var err error
	var expectedNumArgs int = 3
	scoreFlags := flag.NewFlagSet("score", flag.ExitOnError)

	err = scoreFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	scoreFlags.Usage = func() { scoreUsage(scoreFlags) }

	if len(scoreFlags.Args()) != expectedNumArgs {
		scoreFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(scoreFlags.Args()))
	}
	scoreType := scoreFlags.Arg(0)
	inFileList := scoreFlags.Arg(1)
	outFile := scoreFlags.Arg(2)
	s := ScoreSettings{
		ScoreType:  scoreType,
		InFileList: inFileList,
		OutFile:    outFile,
	}
	samAssemblerScore(s)
}

// samAssemblerScore provides tools for scoring diploid assembly results from samAssembler.
// Currently, two modes are supported: baseMatrix, and baseMatrixByRefBase. The details of
// these modes are described in their respective helper functions below.
func samAssemblerScore(s ScoreSettings) {
	switch s.ScoreType {
	case "baseMatrix":
		baseMatrixByRefBase(s.InFileList, s.OutFile, false)
	case "baseMatrixByRefBase":
		baseMatrixByRefBase(s.InFileList, s.OutFile, true)
	default:
		log.Fatalf("Error: Unrecognized score mode: %s.\n", s.ScoreType)
	}
}

// baseMatrixByRefBase returns a confusion matrix for diploid genotype multiclass classification validation
// of samAssembler results. The inFileList represents a new line delimited list of multiFa files to consider.
// the outFile should specify the ".txt" destination file location.
// When byRefBase is true, this program produces 4 output matrices, one for each reference base state (A, C, G, or T).
func baseMatrixByRefBase(inFileList string, outFile string, byRefBase bool) {
	var i, alnPos int
	var currHeader string
	var err error
	var actualGeno, predGeno sam.DiploidBase
	var records []fasta.Fasta
	//initialize answer matrices
	var answerA = make([][]int, 10)
	for i = range answerA {
		answerA[i] = make([]int, 10)
	}
	var answerC = make([][]int, 10)
	for i = range answerC {
		answerC[i] = make([]int, 10)
	}
	var answerG = make([][]int, 10)
	for i = range answerG {
		answerG[i] = make([]int, 10)
	}
	var answerT = make([][]int, 10)
	for i = range answerT {
		answerT[i] = make([]int, 10)
	}
	var answerMerged = make([][]int, 10)
	for i = range answerMerged {
		answerMerged[i] = make([]int, 10)
	}

	//fill matrix from fasta files
	inFiles := fileio.Read(inFileList)
	for i = range inFiles {
		records = fasta.Read(inFiles[i])
		if !validFasta(records) {
			log.Fatalf("Fasta files must have five entries of the same sequence length. Error in: %v.\n", records)
		}
		for alnPos = range records[0].Seq {
			if records[0].Seq[alnPos] < 4 && records[1].Seq[alnPos] < 5 && records[2].Seq[alnPos] < 5 && records[3].Seq[alnPos] < 5 && records[4].Seq[alnPos] < 5 { // if we don't have a gap, but rather have A, C, G, T, or N in estimates, and no N in ref
				actualGeno = basesToDiploidBase(records[1].Seq[alnPos], records[2].Seq[alnPos])
				predGeno = basesToDiploidBase(records[3].Seq[alnPos], records[4].Seq[alnPos])
				answerMerged[refPhasedIndexMapping(predGeno, records[0].Seq[alnPos])][refPhasedIndexMapping(actualGeno, records[0].Seq[alnPos])]++
				switch records[0].Seq[alnPos] {
				case dna.A:
					answerA[refPhasedIndexMapping(predGeno, dna.A)][refPhasedIndexMapping(actualGeno, dna.A)]++
				case dna.C:
					answerC[refPhasedIndexMapping(predGeno, dna.C)][refPhasedIndexMapping(actualGeno, dna.C)]++
				case dna.G:
					answerG[refPhasedIndexMapping(predGeno, dna.G)][refPhasedIndexMapping(actualGeno, dna.G)]++
				case dna.T:
					answerT[refPhasedIndexMapping(predGeno, dna.T)][refPhasedIndexMapping(actualGeno, dna.T)]++
				default:
					log.Fatalf("Error: Unrecognized refBase: %v.\n", records[0].Seq[alnPos])
				}
			}
		}
	}

	var rowNames []string = []string{"HomoRef", "HetRefTs", "HetRefTv1", "HetRefTv2", "HomoTs", "HetTsTv1", "HetTsTv2", "HomoTv1", "HetTv1Tv2", "HomoTv2"}

	//write answer
	out := fileio.EasyCreate(outFile)
	if !byRefBase {
		currHeader = "X\tHomoRef\tHetRefTs\tHetRefTv1\tHetRefTv2\tHomoTs\tHetTsTv1\tHetTsTv2\tHomoTv1\tHetTv1Tv2\tHomoTv2\n"
		writeMatrix(out, answerMerged, rowNames, currHeader)
	} else {
		currHeader = "Ref:A\tHomoRef\tHetRefTs\tHetRefTv1\tHetRefTv2\tHomoTs\tHetTsTv1\tHetTsTv2\tHomoTv1\tHetTv1Tv2\tHomoTv2\n"
		writeMatrix(out, answerA, rowNames, currHeader)
		currHeader = "Ref:C\tHomoRef\tHetRefTs\tHetRefTv1\tHetRefTv2\tHomoTs\tHetTsTv1\tHetTsTv2\tHomoTv1\tHetTv1Tv2\tHomoTv2\n"
		writeMatrix(out, answerC, rowNames, currHeader)
		currHeader = "Ref:G\tHomoRef\tHetRefTs\tHetRefTv1\tHetRefTv2\tHomoTs\tHetTsTv1\tHetTsTv2\tHomoTv1\tHetTv1Tv2\tHomoTv2\n"
		writeMatrix(out, answerG, rowNames, currHeader)
		currHeader = "Ref:T\tHomoRef\tHetRefTs\tHetRefTv1\tHetRefTv2\tHomoTs\tHetTsTv1\tHetTsTv2\tHomoTv1\tHetTv1Tv2\tHomoTv2\n"
		writeMatrix(out, answerT, rowNames, currHeader)
	}

	writeSummaryStatistics(out, answerMerged)

	err = out.Close()
	exception.PanicOnErr(err)
}

// writeSummaryStatistics is a helpr function of baseMatrixByRefBase that writes some summary statistics to the output file.
func writeSummaryStatistics(out *fileio.EasyWriter, dataMerged [][]int) {
	var row, column int
	var matrixTotal int = 0
	for row = 0; row < len(dataMerged); row++ {
		for column = 0; column < len(dataMerged[row]); column++ {
			matrixTotal += dataMerged[row][column]
		}
	}
	// the number of correct predictions is the sum of diagonal elements
	var correctTotal int = 0
	for row = 0; row < len(dataMerged); row++ {
		correctTotal += dataMerged[row][row]
	}

	// trueNegative is the number of true homoRef predicted as homoRef
	var trueNegative = dataMerged[0][0]

	// falseNegative is the number of true mutations misclassified as homoRef
	var falseNegative int
	for column = 1; column < len(dataMerged); column++ {
		falseNegative += dataMerged[0][column]
	}

	// falsePositive is the number of trueHomoRef classified as a mutation
	var falsePositive int
	for row = 1; row < len(dataMerged); row++ {
		falsePositive += dataMerged[row][0]
	}

	// misclassified is the number of mutations (!homoRef) misclassified as a different mutation
	var misclassified int
	for row = 1; row < len(dataMerged); row++ {
		for column = 1; column < len(dataMerged); column++ {
			if row != column {
				misclassified += dataMerged[row][column]
			}
		}
	}

	// actual negative is all positions with homoref genome, or first column
	var actualNegative int
	for row = 0; row < len(dataMerged); row++ {
		actualNegative += dataMerged[row][0]
	}
	var actualPositive = matrixTotal - actualNegative

	// predicted negative is all positions with homoref prediction, or first row
	var predNegative int
	for column = 0; column < len(dataMerged); column++ {
		predNegative += dataMerged[0][column]
	}
	var predPositive = matrixTotal - predNegative
	var accuracy = float64(correctTotal) / float64(matrixTotal)
	var recall = float64(correctTotal-dataMerged[0][0]) / float64(actualPositive)
	var precision = float64(correctTotal-dataMerged[0][0]) / float64(predPositive)
	var misclassificationRate = float64(misclassified) / float64(misclassified+correctTotal-dataMerged[0][0])
	var falseVariantRate = float64(falsePositive) / float64(falsePositive+dataMerged[0][0])
	var falseReferenceRate = float64(falseNegative) / float64(falseNegative+dataMerged[0][0])

	var err error
	_, err = fmt.Fprintf(out, "\nmatrixTotal\tcorrectTotal\ttrueNegative\tfalseNegative\tfalsePositive\tmisclassified\tAccuracy\tInaccuracy\tRecall\tPrecision\tMisclassificaitonRate\tFalseVariantRate\tfalseReferenceRate\n")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "%v\t%v\t%v\t%v\t%v\t%v\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", matrixTotal, correctTotal, trueNegative, falseNegative, falsePositive, misclassified, accuracy, 1-accuracy, recall, precision, misclassificationRate, falseVariantRate, falseReferenceRate)
	exception.PanicOnErr(err)
}

// writeMatrix writes a two-dimensional data matrix of type [][]float64 to a file, with
// a header specified by an input string and rowNames specified by a []string. The data matrix must be of dimensions
// 10x10, and the program thus expects 10 rowNames.
func writeMatrix(out *fileio.EasyWriter, data [][]int, rowNames []string, header string) {
	var i, j int
	var err error

	if len(data) != 10 {
		log.Fatalf("Error: writeMatrix expects an input data matrix with 10 rows. Found: %v.\n", len(data))
	}
	if len(rowNames) != 10 {
		log.Fatalf("Error: writeMatrix expects 10 rownames. Found %v.\n", len(rowNames))
	}

	_, err = fmt.Fprintf(out, header)
	exception.PanicOnErr(err)
	for i = 0; i < 10; i++ {
		if len(data[i]) != 10 {
			log.Fatalf("Error: writeMatrix expects an input data matrix with 10 columns. Found %v columns on row %v.\n", len(data[i]), i)
		}
		_, err = fmt.Fprintf(out, "%s\t", rowNames[i])
		exception.PanicOnErr(err)
		for j = 0; j < 10; j++ {
			_, err = fmt.Fprintf(out, "%v\t", data[i][j])
			exception.PanicOnErr(err)
		}
		_, err = fmt.Fprintf(out, "\n")
		exception.PanicOnErr(err)
	}
}

// refPhasedIndexMapping takes a reference base (A,C,G,or T) and a DiploidBase genotype (AA, AC, AG, ...)
// and returns an index between 0 and 9. These indices correspond to the following mutation types in order:
// 0 - homozygous reference
// 1 - heterozygous reference with transition
// 2 - heterozygous reference with transversion 1 (defined alphanumerically, so for refBase A, tv1 = C and tv2 = T)
// 3 - heterozygous reference with transversion 2
// 4 - homozygous transition
// 5 - heterozygous transition / transversion 1
// 6 - heterozygous transition/ transversion 2
// 7 - heterozygous transversion 1
// 8 - heterozygous transversion 1 / transversion 2
// 9 - homozygous transversion 2.
func refPhasedIndexMapping(geno sam.DiploidBase, refBase dna.Base) int {
	switch refBase {
	case dna.A:
		switch geno {
		case sam.AA:
			return 0
		case sam.AG:
			return 1
		case sam.AC:
			return 2
		case sam.AT:
			return 3
		case sam.GG:
			return 4
		case sam.CG:
			return 5
		case sam.GT:
			return 6
		case sam.CC:
			return 7
		case sam.CT:
			return 8
		case sam.TT:
			return 9
		default:
			log.Fatalf("Error: Unrecognized genotype: %v.\n", geno)
			return -1
		}
	case dna.C:
		switch geno {
		case sam.CC:
			return 0
		case sam.CT:
			return 1
		case sam.AC:
			return 2
		case sam.CG:
			return 3
		case sam.TT:
			return 4
		case sam.AT:
			return 5
		case sam.GT:
			return 6
		case sam.AA:
			return 7
		case sam.AG:
			return 8
		case sam.GG:
			return 9
		default:
			log.Fatalf("Unrecognized genotype: %v.\n", geno)
			return -1
		}
	case dna.G:
		switch geno {
		case sam.GG:
			return 0
		case sam.AG:
			return 1
		case sam.CG:
			return 2
		case sam.GT:
			return 3
		case sam.AA:
			return 4
		case sam.AC:
			return 5
		case sam.AT:
			return 6
		case sam.CC:
			return 7
		case sam.CT:
			return 8
		case sam.TT:
			return 9
		default:
			log.Fatalf("Unrecognized genotype: %v.\n", geno)
			return -1
		}
	case dna.T:
		switch geno {
		case sam.TT:
			return 0
		case sam.CT:
			return 1
		case sam.AT:
			return 2
		case sam.GT:
			return 3
		case sam.CC:
			return 4
		case sam.AC:
			return 5
		case sam.CG:
			return 6
		case sam.AA:
			return 7
		case sam.AG:
			return 8
		case sam.GG:
			return 9
		default:
			log.Fatalf("Error: Unrecognized genotype: %v.\n", geno)
			return -1
		}
	default:
		log.Fatalf("Error: Unrecognized refBase: %v.\n", refBase)
		return -1
	}
}

// basesToDiploidBase converts two input dna.Base variables into a sam.DiploidBase variable.
func basesToDiploidBase(base1 dna.Base, base2 dna.Base) sam.DiploidBase {
	switch base1 {
	case dna.A:
		switch base2 {
		case dna.A:
			return sam.AA
		case dna.C:
			return sam.AC
		case dna.G:
			return sam.AG
		case dna.T:
			return sam.AT
		case dna.N:
			return sam.NN
		case dna.Gap:
			return sam.NN
		default:
			log.Fatalf("Error: Unrecognized base: %v.\n", base2)
		}
	case dna.C:
		switch base2 {
		case dna.A:
			return sam.AC
		case dna.C:
			return sam.CC
		case dna.G:
			return sam.CG
		case dna.T:
			return sam.CT
		case dna.N:
			return sam.NN
		case dna.Gap:
			return sam.NN
		default:
			log.Fatalf("Error: Unrecognized base: %v.\n", base2)
		}
	case dna.G:
		switch base2 {
		case dna.A:
			return sam.AG
		case dna.C:
			return sam.CG
		case dna.G:
			return sam.GG
		case dna.T:
			return sam.GT
		case dna.N:
			return sam.NN
		case dna.Gap:
			return sam.NN
		default:
			log.Fatalf("Error: Unrecognized base: %v.\n", base2)
		}
	case dna.T:
		switch base2 {
		case dna.A:
			return sam.AT
		case dna.C:
			return sam.CT
		case dna.G:
			return sam.GT
		case dna.T:
			return sam.TT
		case dna.N:
			return sam.NN
		case dna.Gap:
			return sam.NN
		}
	case dna.N:
		return sam.NN
	case dna.Gap:
		return sam.NN
	default:
		log.Fatalf("Error: Unrecognized base: %v.\n", base1)
		return sam.NN
	}
	log.Fatalf("Error: Something went wrong.")
	return sam.NN
}

// validFasta checks if an input multiFa has five records of equal length. This is required for scoring.
func validFasta(records []fasta.Fasta) bool {
	if len(records) != 5 {
		return false
	}
	lenRef := len(records[0].Seq)
	if len(records[1].Seq) != lenRef {
		return false
	}
	if len(records[2].Seq) != lenRef {
		return false
	}
	if len(records[3].Seq) != lenRef {
		return false
	}
	if len(records[4].Seq) != lenRef {
		return false
	}
	return true
}
