package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

func samAssemblerScore(scoreType string, inFileList string, outFile string) {
	switch scoreType {
	case "smallBaseMatrix":
		smallBaseMatrix(inFileList, outFile)
	case "baseMatrixByRefBase":
		baseMatrixByRefBase(inFileList, outFile)
	}
}

func smallBaseMatrix(inFileList string, outFile string) {
	var err error
	var i, j, alnPos int
	var records []fasta.Fasta
	var actualGeno, predGeno sam.DiploidBase
	//initialize answer
	var answer = make([][]float64, 11)
	for i = range answer {
		answer[i] = make([]float64, 11)
	}

	//fill matrix from fasta files
	inFiles := fileio.Read(inFileList)
	for i = range inFiles {
		records = fasta.Read(inFiles[i])
		if !validFasta(records) {
			log.Fatalf("Fasta files must have five entries of the same sequence length. Error in: %v.\n", records)
		}
		for alnPos = range records[0].Seq {
			if records[0].Seq[alnPos] < 5 && records[1].Seq[alnPos] < 5 && records[2].Seq[alnPos] < 5 && records[3].Seq[alnPos] < 5 && records[4].Seq[alnPos] < 5 { // if we don't have a gap, but rather have A, C, G, T, or N
				actualGeno = basesToDiploidBase(records[1].Seq[alnPos], records[2].Seq[alnPos])
				predGeno = basesToDiploidBase(records[3].Seq[alnPos], records[4].Seq[alnPos])
				answer[predGeno][actualGeno]++
			}
		}
	}

	//write answer to file
	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\tNN\n")
	exception.PanicOnErr(err)

	for i = 0; i < 11; i++ {
		_, err = fmt.Fprintf(out, "%s\t", sam.DiploidBaseString(sam.DiploidBase(i)))
		exception.PanicOnErr(err)
		for j = 0; j < 11; j++ {
			_, err = fmt.Fprintf(out, "%v\t", answer[i][j])
			exception.PanicOnErr(err)
		}
		_, err = fmt.Fprintf(out, "\n")
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

//baseMatrixByRefBase will make a matrix with elements Aij where:
// i = the actual base, with index i = OutputDiploidGeno * (1+InputDiploidGeno)
func baseMatrixByRefBase(inFileList string, outFile string) {
	var i, j, alnPos, refBase int
	var err error
	var actualGeno, predGeno sam.DiploidBase
	var records []fasta.Fasta
	//initialize answer
	var answer = make([][]float64, 44)
	for i = range answer {
		answer[i] = make([]float64, 44)
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
				refBase = 11 * int(records[0].Seq[alnPos])
				actualGeno = basesToDiploidBase(records[1].Seq[alnPos], records[2].Seq[alnPos])
				predGeno = basesToDiploidBase(records[3].Seq[alnPos], records[4].Seq[alnPos])
				answer[int(predGeno)+refBase][int(actualGeno)+refBase]++
			}
		}
	}

	//write answer
	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "\t")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, strings.Repeat("\tA", 11))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, strings.Repeat("\tC", 11))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, strings.Repeat("\tG", 11))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, strings.Repeat("\tT", 11))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "\n")
	exception.PanicOnErr(err)

	_, err = fmt.Fprintf(out, "\t")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, strings.Repeat("\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\tNN", 4))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "\n")
	exception.PanicOnErr(err)

	for i = range answer { //for each row of answer
		if i/11 <= 0 {
			_, err = fmt.Fprintf(out, "A\t%s\t", sam.DiploidBaseString(sam.DiploidBase(i)))
		} else if i/22 <= 0 {
			_, err = fmt.Fprintf(out, "C\t%s\t", sam.DiploidBaseString(sam.DiploidBase(i-11)))
		} else if i/33 <= 0 {
			_, err = fmt.Fprintf(out, "G\t%s\t", sam.DiploidBaseString(sam.DiploidBase(i-22)))
		} else {
			_, err = fmt.Fprintf(out, "T\t%s\t", sam.DiploidBaseString(sam.DiploidBase(i-33)))
		}
		for j = range answer[i] {
			_, err = fmt.Fprintf(out, "%v\t", answer[i][j])
			exception.PanicOnErr(err)
		}
		_, err = fmt.Fprintf(out, "\n")
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

//TODO: Return N is a temp feature
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
			log.Fatalf("Unrecognized base: %v.\n", base2)
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
			log.Fatalf("Unrecognized base: %v.\n", base2)
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
			log.Fatalf("Unrecognized base: %v.\n", base2)
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
		log.Fatalf("Unrecognized base: %v.\n", base1)
		return sam.NN
	}
	log.Fatalf("Something went wrong.")
	return sam.NN
}

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
