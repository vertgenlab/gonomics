// Command Group: "Sequence Evolution & Reconstruction"

// Reconstruct ancient sequences using extant genomes and a newick tree with branch lengths
package main

import (
	"flag"
	"fmt"
	"log"
	"bufio"
	"os"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/wig"
)

type IlsReconstructSeqSettings struct {
	PostProbFiles string
	ReconFiles string
	OutDir string
	Precision float32
}

// IlsReconstructSeqUsage defines the usage statement for the ilsReconstructSeq command
func IlsReconstructSeqUsage(IlsReconstructSeqFlags *flag.FlagSet) {
	fmt.Print(
		"ilsReconstructSeq accounts for incomplete lineage sorting in ancestral sequence reconstruction, by weighting the input sequence reconstructions by the posterior probability at that position of the associated topology (out of four). Memory scales by pDNA size, so if working with large genomes, make sure to split by chromosome." + 
			"This program takes as input a txt file containing the addresses of n wig files, and a txt file containing the addresses of n corresponding pFastas."
			"This program returns a pfasta file containing a sequence representing the ----------------" + 
			"Usage:\n" + 
			"ilsReconstructSeq <posterior-probabilities>.txt <recons>.txt outDir\n" + 
			"options:\n")
	IlsReconstructSeqFlags.PrintDefaults()
}

// parseIlsReconstructSeqArgs is the main function of the IlsReconstructSeq command. It parses options and runs the ilsReconstructSeq function.
func parseIlsReconstructSeqArgs() {
	expectedNumArgs := 3
	var err error
	IlsReconstructSeqFlags := flag.NewFlatSet("IlsReconstructSeq", flag.ExitOnError)
	var precision *float32 = IlsReconstructSeqFlags.float32("precision", 0.001, "Specify the precision to use when checking for valid pDNA base (sums to 1).")

	err = IlsReconstructSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	IlsReconstructSeqFlags.Usage = func() { IlsReconstructSeqUsage(IlsReconstructSeqFlags)}

	if len(IlsReconstructSeqFlags.Args()) != expectedNumArgs {
		IlsReconstructSeqFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(IlsReconstructSeqFlags.Args()))
	}

	postProbFiles := IlsReconstructSeqFlags.Arg(0)
	reconFiles := IlsReconstructSeqFlags.Arg(1)
	outDir := IlsReconstructSeqFlags.Arg(2)

	s := IlsReconstructSeqSettings{
		PostProbFiles:	postProbFiles,
		ReconFiles: reconFiles,
		OutDir: outDir,
		Precision: precision,
	}

	ilsReconstructSeq(s)
}

func ilsReconstructSeq(s IlsReconstructSeqSettings) {
	readPostProbs, err := os.Open(s.PostProbFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", s.PostProbFiles)
	}

	fileScanner := bufio.NewScanner(readPostProbs)
	fileScanner.Split(bufio.ScanLines)
	var postProbFileLines []string

	for fileScanner.Scan() {
		postProbFileLines = append(postProbFileLines, fileScanner.Text())
	}

	readPostProbs.Close()
	postProbs := [][]float64

	for _, filepath := range postProbfileLines {
		postProbs = append(postProbs, wig.Read(filepath))
	}


	readRecons, err := os.Open(s.ReconFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", s.ReconFiles)
	}

	fileScanner := bufio.NewScanner(readRecons)
	fileScanner.Split(bufio.ScanLines)
	var ReconFileLines []string

	for fileScanner.Scan() {
		ReconFileLines = append(ReconFileLines, fileScanner.Text())
	}

	readRecons.Close()
	recons := []pFasta.PFasta

	for _, filepath := range ReconFileLines {
		recons = append(recons, pFasta.Read(filepath))
	}

	out := reconstruct.IlsReconstructSeq(postProbs, recons, s.Precision)

	pFasta.Write(s.OutDir, out)
}