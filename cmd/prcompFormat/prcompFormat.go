// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

/*
	The first loop checks each position in the specified range, with the inner loop checking each
	species to determine if any species contain a gap or an N at that position or if all of the
	species have the same base at that position. If there are no gaps or Ns and at least one species
	has an alternative base position, that column of bases is added to the output fasta.

	Next, the output fasta, once constructed, is converted into a binary matrix for input for PCA.
*/

func prcompFormat(infile string, outfile string, verbose bool) {
	records := fasta.Read(infile)
	fasta.AllToUpper(records)
	var subFa = make([]fasta.Fasta, len(records))
	var allMatch, allValid bool
	var currentBase dna.Base
	var err error

	//add species name to subFa
	for i := range records {
		subFa[i].Name = records[i].Name
	}

	for i := 0; i < len(records[0].Seq); i++ {
		//check for gaps, identity, or Ns
		currentBase = records[0].Seq[i]
		allMatch = true
		allValid = true
		if !(currentBase == dna.N || currentBase == dna.Gap) {
			for j := 1; j < len(records); j++ {
				if records[j].Seq[i] != currentBase && records[j].Seq[i] != dna.Gap && records[j].Seq[i] != dna.N {
					if verbose {
						_, err = fmt.Printf("i = %d, j = %d\n", i, j)
						exception.PanicOnErr(err)
					}
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
	var outline string

	//print header
	var headerString string = "Sample"
	for n := 0; n < (4 * len(subFa[0].Seq)); n++ {
		headerString = headerString + fmt.Sprintf("\tVar_%d", n)
	}

	_, err = fmt.Fprintf(file, "%s\n", headerString)
	exception.PanicOnErr(err)

	for i := 0; i < len(subFa); i++ {
		outline = subFa[i].Name
		if verbose {
			_, err = fmt.Println(len(subFa[i].Seq))
			exception.PanicOnErr(err)
		}
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
				log.Fatalf("Base: %d\n", subFa[i].Seq[j])
			}
		}
		outline = outline + "\n"
		_, err = fmt.Fprintf(file, "%s", outline)
		exception.PanicOnErr(err)
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"prcompFormat - Generates a binary input matrix for PCA.\n" +
			"Please see https://github.com/vertgenlab/vglDocumentation for a detailed explanation of\n" +
			"how to turn the output text file from this program into a PCA plot using R.\n" +
			"Usage:\n" +
			"prcompFormat infile.fa outfile.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var verbose *bool = flag.Bool("verbose", false, "enable debug prints.")

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

	prcompFormat(infile, outfile, *verbose)
}
