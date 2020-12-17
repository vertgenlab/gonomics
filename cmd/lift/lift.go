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
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

const Verbose int = 0

func lift(chainFile string, inFile string, outFile string, faFile string, unMapped string, minMatch float64) {
	//TODO: minMatch must be between 0 and 1
	//first task, make tree from chainFile
	chainChan, _ := chain.GoReadToChan(chainFile)
	var chainIntervals []interval.Interval
	for val := range chainChan {
		chainIntervals = append(chainIntervals, val)
	}
	tree := interval.BuildTree(chainIntervals)
	out := fileio.EasyCreate(outFile)
	defer out.Close()

	un := fileio.EasyCreate(unMapped)
	defer un.Close()

	var records []*fasta.Fasta
	var currVcf *vcf.Vcf
	var a, b float64

	if faFile != "" {
		if !vcf.IsVcfFile(inFile) {
			log.Fatalf("Fasta files should be provided when lifting VCF files. inFile is not a valid VCF file and must end in either .vcf or .vcf.gz.")
		}
		records = fasta.Read(faFile)
	}

	//second task, read in intervals, find chain, and convert to new interval
	inChan := interval.GoReadToLiftChan(inFile)
	var overlap []interval.Interval
	for i := range inChan {
		overlap = interval.Query(tree, i, "any")//TODO: verify proper t/q contains
		if len(overlap) > 1 {
			fmt.Fprintf(un, "Record below maps multiple chains:\n")
			i.WriteToFileHandle(un)
		} else if len(overlap) == 0 {
			fmt.Fprintf(un, "Record below has no ortholog in new assembly:\n")
			i.WriteToFileHandle(un)
		} else if !minMatchPass(overlap[0].(*chain.Chain), i, minMatch) {
			a, b = interval.MatchProportion(overlap[0].(*chain.Chain), i)
			fmt.Fprintf(un, "Record below fails minMatch with a proportion of %f:\n", numbers.MinFloat64(a, b))
			fmt.Fprintf(un, "Here's the corresponding chain: Score: %d.\n", overlap[0].(*chain.Chain).Score)
			i.WriteToFileHandle(un)
		} else {
			//DEBUG: interval.PrettyPrint(i)
			i.UpdateLift(interval.LiftIntervalWithChain(overlap[0].(*chain.Chain), i))
			//DEBUG: interval.PrettyPrint(i)
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
					fmt.Fprintf(un, "Record below was lifted, but the ref and alt alleles are inverted:\n")
					i.WriteToFileHandle(un)
					vcf.InvertVcf(currVcf)
				} else {
					log.Fatalf("Neither the Ref nor the Alt allele matched the bases in the corresponding destination fasta location.")
				}
			}
			i.WriteToFileHandle(out)
		}
	}
}

//minMatchPass returns true if the interval/chain has over a 95 percent base match, false otherwise.
func minMatchPass(overlap *chain.Chain, i interval.Interval, minMatch float64) bool {
	a, b := interval.MatchProportion(overlap, i)
	if a < minMatch {
		return false
	}
	if b < minMatch {
		return false
	}
	return true
}

func usage() {
	fmt.Print(
		"lift - Moves an interval interface compatable file format between assemblies.\n" +
			"Usage:\n" +
			"lift lift.chain inFile outFile unMapped\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var faFile *string = flag.String("faFile", "", "Specify a fasta file for Lifting VCF format files. This fasta should correspond to the destination assembly.")
	var minMatch *float64 = flag.Float64("minMatch", 0.95, "Specify the minimum proportion of matching bases required for a successful lift operation.")

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
	unMapped := flag.Arg(3)
	lift(chainFile, inFile, outFile, *faFile, unMapped, *minMatch)
}
