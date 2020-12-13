package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/interval"
	"log"
)

const Verbose int = 0

func lift(chainFile string, inFile string, outFile string, faFile string) {
	//first task, make tree from chainFile
	chainChan, _ := chain.GoReadToChan(chainFile)
	var chainIntervals []interval.Interval
	for val := range chainChan {
		chainIntervals = append(chainIntervals, val)
	}
	tree := interval.BuildTree(chainIntervals)
	out := fileio.EasyCreate(outFile)
	var records []*fasta.Fasta
	var currVcf *vcf.Vcf

	if faFile != "" {
		records = fasta.Read(faFile)
	}

	//second task, read in intervals, find chain, and convert to new interval
	inChan := interval.GoReadToLiftChan(inFile)
	var overlap []interval.Interval
	for i := range inChan {
		overlap = interval.Query(tree, i, "any")
		if len(overlap) > 1 {
			log.Fatalf("Multiple overlaps???")
		}
		interval.PrettyPrint(i)
		i.UpdateLift(interval.LiftIntervalWithChain(overlap[0].(*chain.Chain), i))
		interval.PrettyPrint(i)
		//special check for lifting over VCF files
		if faFile != "" {
			//faFile will be given if we are lifting over VCF data.
			currVcf = i.(*vcf.Vcf)

			//first question: does the "Ref" match the destination fa at this position.
			if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos - 1), dna.StringToBases(currVcf.Ref)) {
				//second question: does the "Alt" also match. Can occur in corner cases such as Ref=A, Alt=AAA. Currently we don't invert but write a verbose log print.
				if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos - 1), dna.StringToBases(currVcf.Alt)) && Verbose > 0 {
					log.Printf("For VCF on %s at position %d, Alt and Ref both match the fasta. Ref: %s. Alt: %s.", currVcf.Chr, currVcf.Pos, currVcf.Ref, currVcf.Alt)
				}
			//the third case handles when the alt matches but not the ref, in which case we invert the VCF.
			} else if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos - 1), dna.StringToBases(currVcf.Alt)) {
				vcf.InvertVcf(currVcf)
			} else {
				log.Fatalf("Neither the Ref nor the Alt allele matched the bases in the corresponding destination fasta location.")
			}
		}
		i.WriteToFileHandle(out)
	}
}

func usage() {
	fmt.Print(
		"lift - Moves an interval interface compatable file format between assemblies.\n" +
			"Usage:\n" +
			"lift lift.chain inFile outFile\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var faFile *string = flag.String("faFile", "", "Specify a fasta file for Lifting VCF format files. This fasta should correspond to the destination assembly.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	chainFile := flag.Arg(0)
	inFile := flag.Arg(1)
	outFile := flag.Arg(2)
	lift(chainFile, inFile, outFile, *faFile)
}
