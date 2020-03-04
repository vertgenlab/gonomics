package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

/*
	The first loop checks each position in the specified range, with the inner loop checking each
	species to determine if any species contain a gap or an N at that position or if all of the
	species have the same base at that position. If there are no gaps or Ns and at least one species
	has an alternative base position, that column of bases is added to the output fasta.

	Next, the output fasta, once constructed, is converted into a binary matrix for input for PCA.
*/

func prcompFormat(infile string, outfile string) {
	records := fasta.Read(infile)
	fasta.AllToUpper(records)
	var subFa = make([]fasta.Fasta, len(records))
	var allMatch, allValid bool
	var currentBase dna.Base
	//fmt.Printf("Len Seq: %d\n", len(records[0].Seq))
	//fmt.Printf("Len Rec: %d\n", len(records))

	//add species name to subFa
	for i := 0; i < len(records); i++ {
		subFa[i].Name = records[i].Name
	}

	for i := 0; i < len(records[0].Seq); i++ {
		//check for gaps, identity, or Ns
		currentBase = records[0].Seq[i]
		allMatch = true
		allValid = true
		if !(currentBase == dna.N || currentBase == dna.Gap) {
			for j := 1; j < len(records); j++ {
				//fmt.Printf("i = %d, j = %d\n", i, j)
				if records[j].Seq[i] != currentBase && records[j].Seq[i] != dna.Gap && records[j].Seq[i] != dna.N {
					fmt.Printf("i = %d, j = %d\n", i, j)
					allMatch = false
				}
				if records[j].Seq[i] == dna.Gap || records[j].Seq[i] == dna.N {
					allValid = false
				}
			}
			//now we add to output fasta
			if !allMatch && allValid {
				for k := 0; k < len(records); k++ {
					subFa[k].Seq = append(subFa[k].Seq, records[k].Seq[i])
				}
			}
		}
	}

	//expand bases to binary values while writing to the outfile
	file := fileio.EasyCreate(outfile)
	defer file.Close()
	var outline string

	//print header
	var headerstring string = "Sample"
	for n := 0; n < (4 * len(subFa[0].Seq)); n++ {
		headerstring = headerstring + fmt.Sprintf("\tVar_%d", n)
	}

	_, err := fmt.Fprintf(file, "%s\n", headerstring)
	common.ExitIfError(err)

	for i := 0; i < len(subFa); i++ {
		outline = subFa[i].Name
		fmt.Println(len(subFa[i].Seq))
		for j := 0; j < len(subFa[0].Seq); j++ {
			if subFa[i].Seq[j] == dna.A {
				outline = outline + "\t1\t0\t0\t0"
			} else if subFa[i].Seq[j] == dna.C {
				outline = outline + "\t0\t1\t0\t0"
			} else if subFa[i].Seq[j] == dna.G {
				outline = outline + "\t0\t0\t1\t0"
			} else if subFa[i].Seq[j] == dna.T {
				outline = outline + "\t0\t0\t0\t1"
			} else {
				fmt.Printf("Base: %d\n", subFa[i].Seq[j])
				common.ExitIfError(fmt.Errorf("Critical Failure!!\n"))
			}
		}
		outline = outline + "\n"
		_, err := fmt.Fprintf(file, "%s", outline)
		common.ExitIfError(err)
	}
}

func usage() {
	fmt.Print(
		"prcompFormat - Generates a binary input matrix for PCA in R.\n" +
			"Usage:\n" +
			"prcompFormat infile.fa outfile.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	prcompFormat(infile, outfile)
}
