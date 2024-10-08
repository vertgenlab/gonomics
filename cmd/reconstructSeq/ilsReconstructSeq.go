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

type Settings struct {
	PostProbFiles string
	ReconFiles string
	OutFile string
}

func ilsReconstructSeq(s Settings) {
	readPostProbs, err := os.Open(PostProbFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", PostProbFiles)
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


	readRecons, err := os.Open(ReconFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", ReconFiles)
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

	out := reconstruct.IlsReconstructSeq(postProbs, recons)

	pFasta.Write(Outfile, out)
}

func usage() {
	fmt.Print(
		"ilsReconstructSeq accounts for incomplete lineage sorting in ancestral sequence reconstruction, by weighting the input sequence reconstructions by the posterior probability at that position of the associated topology (out of four). Memory scales by pDNA size, so if working with large genomes, make sure to split by chromosome." + 
			"This program takes as input a txt file containing the addresses of n wig files, and a txt file containing the addresses of n corresponding pFastas."
			"This program returns a pfasta file containing a sequence representing the ----------------" + 
			"Usage:\n" + 
			"ilsReconstructSeq <posterior-probabilities>.txt <recons>.txt"
	)
}