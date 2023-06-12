package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func mergeMultiFa(inAFile string, inBFile string, outFile string) {
	var inA []fasta.Fasta = fasta.Read(inAFile)
	var inB []fasta.Fasta = fasta.Read(inBFile)
	var i int

	if inA[0].Name != inB[0].Name {
		log.Fatalf("The first file reference name (%s) does not match the second file reference name (%s).", inA[0].Name, inB[0].Name)
	}

	if len(inA) < 2 {
		log.Fatalf("The first multiFa file has less than two entries, and is therefore not a valid multiFa.")
	}
	if len(inB) < 2 {
		log.Fatalf("The second multiFa file has less than two entries, and is therefore not a valid multiFa.")
	}
	for i = 1; i < len(inA); i++ {
		if len(inA[0].Seq) != len(inA[i].Seq) {
			log.Fatalf("In the first multiFa file, entry at index %v is not the same length as the reference sequence. This file is therefore not a valid multiFa.", i)
		}
	}

	for i = 1; i < len(inB); i++ {
		if len(inB[0].Seq) != len(inB[i].Seq) {
			log.Fatalf("In second first multiFa file, entry at index %v is not the same length as the reference sequence. This file is therefore not a valid multiFa.", i)
		}
	}

	//we'll build the empty output sequence
	var answerLen = len(inA) + len(inB) - 1 //the answer length is the sum of the number of entries in the two files, minus 1, for the shared reference.
	var answer = make([]fasta.Fasta, answerLen)
	for i = range inA {
		answer[i].Name = inA[i].Name
	}
	for i = 1; i < len(inB); i++ {
		answer[i+len(inA)-1].Name = inB[i].Name
	}

	var alnPosA, alnPosB int = 0, 0

	for alnPosA < len(inA[0].Seq) && alnPosB < len(inB[0].Seq) {
		if (inA[0].Seq[alnPosA] != dna.Gap && inB[0].Seq[alnPosB] != dna.Gap) || (inA[0].Seq[alnPosA] == dna.Gap && inB[0].Seq[alnPosB] == dna.Gap) {
			if inA[0].Seq[alnPosA] != inB[0].Seq[alnPosB] {
				log.Fatalf("Error in mergeMultiFa. Reference sequences at alignment position %v of the first file showed different bases between the two files. Found in first file: %s. Found in second file: %s.\n", alnPosA, dna.BaseToString(inA[0].Seq[alnPosA]), dna.BaseToString(inB[0].Seq[alnPosB]))
			}
			for i = range inA {
				answer[i].Seq = append(answer[i].Seq, inA[i].Seq[alnPosA])
			}
			for i = 1; i < len(inB); i++ {
				answer[i+len(inA)-1].Seq = append(answer[i+len(inA)-1].Seq, inB[i].Seq[alnPosB])
			}
			alnPosA++
			alnPosB++
		} else if inB[0].Seq[alnPosB] == dna.Gap {
			for i = range inA {
				answer[i].Seq = append(answer[i].Seq, dna.Gap)
			}
			for i = 1; i < len(inB); i++ {
				answer[i+len(inA)-1].Seq = append(answer[i+len(inA)-1].Seq, inB[i].Seq[alnPosB])
			}
			alnPosB++
		} else if inA[0].Seq[alnPosA] == dna.Gap {
			for i = range inA {
				answer[i].Seq = append(answer[i].Seq, inA[i].Seq[alnPosA])
			}
			for i = 1; i < len(inB); i++ {
				answer[i+len(inA)-1].Seq = append(answer[i+len(inA)-1].Seq, dna.Gap)
			}
			alnPosA++
		} else {
			log.Fatalf("Something went wrong.")
		}
	}

	fasta.Write(outFile, answer)
}

func usage() {
	fmt.Print(
		"mergeMultiFa - Merge two multiFa files on a shared reference. Does not perform local realignment on INDELs.\n" +
			"Usage:\n" +
			" mergeMultiFa input1.fa input2.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inAFile := flag.Arg(0)
	inBFile := flag.Arg(1)
	outFile := flag.Arg(2)

	mergeMultiFa(inAFile, inBFile, outFile)
}
