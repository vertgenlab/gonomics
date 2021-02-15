package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"unicode/utf8"
)

const Verbose int = 0

func lift(chainFile string, inFile string, outFile string, faFile string, unMapped string, minMatch float64) {
	if minMatch < 0.0 || minMatch > 1.0 {
		log.Fatalf("minMatch must be between 0 and 1. User input: %f.\n", minMatch)
	}
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
			log.Fatalf("Fasta file is provided but lift file is not a VCF file.")
		}
		records = fasta.Read(faFile)
	}

	//TODO: General GoReadToChan header returns will allow us to avoid opening the file twice.
	if vcf.IsVcfFile(inFile) {
		tmpOpen := fileio.EasyOpen(inFile)
		header := vcf.ReadHeader(tmpOpen)
		vcf.NewWriteHeader(out, header)
		tmpOpen.Close()
	}

	//second task, read in intervals, find chain, and convert to new interval
	inChan := interval.GoReadToLiftChan(inFile)
	var overlap []interval.Interval
	for i := range inChan {
		overlap = interval.Query(tree, i, "any") //TODO: verify proper t/q contains
		if len(overlap) > 1 {
			fmt.Fprintf(un, "Record below maps to multiple chains:\n")
			i.WriteToFileHandle(un)
		} else if len(overlap) == 0 {
			fmt.Fprintf(un, "Record below has no ortholog in new assembly:\n")
			i.WriteToFileHandle(un)
		} else if !minMatchPass(overlap[0].(*chain.Chain), i, minMatch) {
			a, b = interval.MatchProportion(overlap[0].(*chain.Chain), i)
			fmt.Fprintf(un, "Record below fails minMatch with a proportion of %f. Here's the corresponding chain: %d.\n", numbers.MinFloat64(a, b), overlap[0].(*chain.Chain).Score)
			i.WriteToFileHandle(un)
		} else {
			//DEBUG: interval.PrettyPrint(i)
			i.UpdateLift(interval.LiftIntervalWithChain(overlap[0].(*chain.Chain), i)) //now i is 1-based and we can assert VCF.
			//DEBUG: interval.PrettyPrint(i)
			//special check for lifting over VCF files
			if faFile != "" {
				//faFile will be given if we are lifting over VCF data.
				currVcf = i.(*vcf.Vcf)
				if utf8.RuneCountInString(currVcf.Ref) > 1 || utf8.RuneCountInString(currVcf.Alt[0]) > 1 {
					log.Fatalf("VCF liftOver is currently not supported for INDEL records. Please filter the input VCF for substitutions and try again.") //Currently we're causing INDEL records to fatal.
				}
				if len(currVcf.Alt) > 1 {
					log.Fatalf("VCF liftOver is currently only supported for biallelic variants. Variant at %s %v is polyallelic.", currVcf.Chr, currVcf.Pos)
				}
				//first question: does the "Ref" match the destination fa at this position.
				if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Ref)) {
					//second question: does the "Alt" also match. Can occur in corner cases such as Ref=A, Alt=AAA. Currently we don't invert but write a verbose log print.
					if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) && Verbose > 0 {
						fmt.Fprintf(un, "For VCF on %s at position %d, Alt and Ref both match the fasta. Ref: %s. Alt: %s.", currVcf.Chr, currVcf.Pos, currVcf.Ref, currVcf.Alt)
					}
					//the third case handles when the alt matches but not the ref, in which case we invert the VCF.
				} else if fasta.QuerySeq(records, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) {
					fmt.Fprintf(un, "Record below was lifted, but the ref and alt alleles are inverted:\n")
					//DEBUG:log.Printf("currVcf Pos -1: %d. records base: %s.", currVcf.Pos-1, dna.BaseToString(records[0].Seq[int(currVcf.Pos-1)]))
					i.WriteToFileHandle(un)
					vcf.InvertVcf(currVcf)
				} else {
					fmt.Fprintf(un, "For the following record, neither the Ref nor the Alt allele matched the bases in the corresponding destination fasta location.\n")
					i.WriteToFileHandle(un)
				}
			}
			i.WriteToFileHandle(out)
		}
	}
}

//minMatchPass returns true if the interval/chain has over a certain percent base match (minMatch argument is the proportion, default 0.95), false otherwise.
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
			"Warning: For Vcf lift, the original headers are retained in the output without modification. Use output header information at your own risk.\n" +
			"Please note: Vcf lift is not compatable with Unix piping.\n" +
			"Please note: Vcf lift with fa crossreferencing is currently only supported for biallelic and substitution variants.\n" +
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
