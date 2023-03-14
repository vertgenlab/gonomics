// Command Group: "General Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	bedpe2 "github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"strings"
	"unicode/utf8"
)

type fileWriter interface {
	WriteToFileHandle(io.Writer)
}

func liftCoordinates(chainFile string, inFile string, outFile string, faFile string, unMapped string, minMatch float64, swapAlleleAB bool, verbose int, bedpe bool) {
	var err error
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
	un := fileio.EasyCreate(unMapped)

	var ref *fasta.Seeker
	var currVcf vcf.Vcf
	var a, b float64

	if faFile != "" {
		if !vcf.IsVcfFile(inFile) {
			log.Fatalf("Fasta file is provided but lift file is not a VCF file.")
		}
		ref = fasta.NewSeeker(faFile, "")
	}

	if vcf.IsVcfFile(inFile) {
		tmpOpen := fileio.EasyOpen(inFile)
		header := vcf.ReadHeader(tmpOpen)
		vcf.NewWriteHeader(out, header)
		err = tmpOpen.Close()
		exception.PanicOnErr(err)
	}
	// bedpe error message map
	mp := map[int]string{
		0:  "Record below has no ortholog in new assembly.",
		1:  "Record below mapped correctly",
		-1: "Record below falls below min-match threshold.",
	}

	//second task, read in intervals, find chain, and convert to new interval
	var A, B string
	var found1, found2 bool
	var halfA, halfB bed.Bed
	inChan := lift.GoReadToChan(inFile)
	var overlap []interval.Interval
	for i := range inChan {
		if bedpe {
			j := <-inChan
			iOut, mappedI := liftAndCheck(i, tree, faFile, minMatch)
			jOut, mappedJ := liftAndCheck(j, tree, faFile, minMatch)
			halfA = bed.Bed{iOut.GetChrom(), iOut.GetChromStart(), iOut.GetChromEnd()}
			halfB = bed.Bed{jOut.GetChrom(), jOut.GetChromStart(), jOut.GetChromEnd()}
			bedpeEntry := bedpe2.BedPe{halfA, halfB}

			if mappedI == 1 && mappedJ == 1 {
				bedpeEntry.(fileWriter).WriteToFileHandle(out)
			} else {
				A, found1 = mp[mappedI]
				B, found2 = mp[mappedJ]
				if found1 && found2 {
					_, err = fmt.Fprintf(un, "A: %s  B: %s\n", A, B)
					exception.PanicOnErr(err)
					bedpeEntry.(fileWriter).WriteToFileHandle(out)
				} else if found1 && !found2 {
					_, err = fmt.Fprintf(un, "A: %s  B: Record below maps to multiple chains\n", A)
					exception.PanicOnErr(err)
					bedpeEntry.(fileWriter).WriteToFileHandle(out)
				} else if !found1 && found2 {
					_, err = fmt.Fprintf(un, "A: Record below maps to multiple chains  B: %s\n", B)
					exception.PanicOnErr(err)
					bedpeEntry.(fileWriter).WriteToFileHandle(out)
				}

			}

			overlap = interval.Query(tree, i, "any")

			if len(overlap) > 1 {
				_, err = fmt.Fprintf(un, "Record below maps to multiple chains:\n")
				exception.PanicOnErr(err)
				i.WriteToFileHandle(un)
			} else if len(overlap) == 0 {
				_, err = fmt.Fprintf(un, "Record below has no ortholog in new assembly:\n")
				exception.PanicOnErr(err)
				i.WriteToFileHandle(un)
			} else if !minMatchPass(overlap[0].(chain.Chain), i, minMatch) {
				a, b = lift.MatchProportion(overlap[0].(chain.Chain), i)
				_, err = fmt.Fprintf(un, "Record below fails minMatch with a proportion of %f. Here's the corresponding chain: %d.\n", numbers.Min(a, b), overlap[0].(chain.Chain).Score)
				exception.PanicOnErr(err)
				i.WriteToFileHandle(un)
			} else {
				i = i.UpdateCoord(lift.LiftCoordinatesWithChain(overlap[0].(chain.Chain), i)).(lift.Lift) //now i is 1-based and we can assert VCF.
				//special check for lifting over VCF files
				if faFile != "" {
					//faFile will be given if we are lifting over VCF data.
					currVcf = i.(vcf.Vcf)
					if utf8.RuneCountInString(currVcf.Ref) > 1 || utf8.RuneCountInString(currVcf.Alt[0]) > 1 {
						_, err = fmt.Fprintf(un, "The following record did not lift as VCF lift is not currently supported for INDEL records.\n")
						exception.PanicOnErr(err)
						i.WriteToFileHandle(un)
						//log.Fatalf("VCF liftOver is currently not supported for INDEL records. Please filter the input VCF for substitutions and try again.") //Currently we're causing INDEL records to fatal.
					} else if len(currVcf.Alt) > 1 {
						_, err = fmt.Fprintf(un, "The following record did not lift as VCF lift is not currently supported for multiallelic sites.\n")
						exception.PanicOnErr(err)
						i.WriteToFileHandle(un)
					} else if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Ref)) { //first question: does the "Ref" match the destination fa at this position.
						//second question: does the "Alt" also match. Can occur in corner cases such as Ref=A, Alt=AAA. Currently we don't invert but write a verbose log print.
						if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) && verbose > 0 {
							_, err = fmt.Fprintf(un, "For VCF on %s at position %d, Alt and Ref both match the fasta. Ref: %s. Alt: %s.", currVcf.Chr, currVcf.Pos, currVcf.Ref, currVcf.Alt)
							exception.PanicOnErr(err)
						}
						i.(fileWriter).WriteToFileHandle(out)
						//the third case handles when the alt matches but not the ref, in which case we invert the VCF.
					} else if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) {
						_, err = fmt.Fprintf(un, "Record below was lifted, but the ref and alt alleles are inverted:\n")
						exception.PanicOnErr(err)
						i.WriteToFileHandle(un)
						currVcf = vcf.InvertVcf(currVcf)
						if swapAlleleAB {
							swapInfoAlleles(&currVcf)
						}
						i = &currVcf
						i.(fileWriter).WriteToFileHandle(out)
					} else {
						_, err = fmt.Fprintf(un, "For the following record, neither the Ref nor the Alt allele matched the bases in the corresponding destination fasta location.\n")
						exception.PanicOnErr(err)
						i.WriteToFileHandle(un)
					}
				} else {
					i.(fileWriter).WriteToFileHandle(out)
				}
			}
		}
		err = out.Close()
		exception.PanicOnErr(err)
		err = un.Close()
		exception.PanicOnErr(err)
	}
}

// QuerySeq takes in a fasta seeker and a position (name and index) and returns true if a query
// sequence of bases matches the fasta at this position.
// Note: for QuerySeq, RefPosToAlnPos is probably not required if you are using an assembly fasta
// as the reference, but if you are querying from alignment Fasta, you'll want to get the alnIndex
// before calling this function
func QuerySeq(ref *fasta.Seeker, chr string, index int, query []dna.Base) bool {
	fetchSeq, err := fasta.SeekByName(ref, chr, index, index+len(query))
	exception.PanicOnErr(err)
	return dna.CompareSeqsIgnoreCaseAndGaps(query, fetchSeq) == 0
}

// minMatchPass returns true if the interval/chain has over a certain percent base match (minMatch argument is the proportion, default 0.95), false otherwise.
func minMatchPass(overlap chain.Chain, i interval.Interval, minMatch float64) bool {
	a, b := lift.MatchProportion(overlap, i)
	if a < minMatch {
		return false
	}
	if b < minMatch {
		return false
	}
	return true
}

// swapInfoAlelles switches the values of ALLELE_A and ALLELE_B in the info field
func swapInfoAlleles(v *vcf.Vcf) {
	var foundA, foundB bool
	newInfo := []byte(v.Info)
	alleleAidx := strings.Index(v.Info, "ALLELE_A=")
	if alleleAidx == -1 {
		foundA = true
	}
	alleleAidx += len("ALLELE_A=")

	alleleBidx := strings.Index(v.Info, "ALLELE_B=")
	if alleleBidx == -1 {
		foundB = true
	}
	alleleBidx += len("ALLELE_B=")

	if foundA != foundB {
		log.Printf("WARNING: Found ALLELE_A or ALLELE_B in the following record, but not both. Record may be malformed\n%s\n", v)
		return
	}

	newInfo[alleleAidx], newInfo[alleleBidx] = newInfo[alleleBidx], newInfo[alleleAidx]
	v.Info = string(newInfo)
}

func liftAndCheck(in lift.Lift, inChain map[string]*interval.IntervalNode, faFile string, matchPerc float64) (lift.Lift, int) {
	var overlap []interval.Interval
	var a, b float64
	var currVcf vcf.Vcf

	overlap = interval.Query(inChain, in, "any")
	mapped := len(overlap)

	if len(overlap) == 1 {
		if !minMatchPass(overlap[0].(chain.Chain), in, matchPerc) {
			mapped = -1
			a, b = lift.MatchProportion(overlap[0].(chain.Chain), in)
		} else {
			in = in.UpdateCoord(lift.LiftCoordinatesWithChain(overlap[0].(chain.Chain), in)).(lift.Lift)
			//special check for lifting over VCF files
			if faFile != "" {
				//faFile will be given if we are lifting over VCF data.
				currVcf = in.(vcf.Vcf)
			}
			//special check for lifting over VCF files
			if faFile != "" {
				//faFile will be given if we are lifting over VCF data.
				currVcf = in.(vcf.Vcf)
				if utf8.RuneCountInString(currVcf.Ref) > 1 || utf8.RuneCountInString(currVcf.Alt[0]) > 1 {
					_, err = fmt.Fprintf(un, "The following record did not lift as VCF lift is not currently supported for INDEL records.\n")
					exception.PanicOnErr(err)
					i.WriteToFileHandle(un)
					//log.Fatalf("VCF liftOver is currently not supported for INDEL records. Please filter the input VCF for substitutions and try again.") //Currently we're causing INDEL records to fatal.
				} else if len(currVcf.Alt) > 1 {
					_, err = fmt.Fprintf(un, "The following record did not lift as VCF lift is not currently supported for multiallelic sites.\n")
					exception.PanicOnErr(err)
					i.WriteToFileHandle(un)
				} else if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Ref)) { //first question: does the "Ref" match the destination fa at this position.
					//second question: does the "Alt" also match. Can occur in corner cases such as Ref=A, Alt=AAA. Currently we don't invert but write a verbose log print.
					if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) && verbose > 0 {
						_, err = fmt.Fprintf(un, "For VCF on %s at position %d, Alt and Ref both match the fasta. Ref: %s. Alt: %s.", currVcf.Chr, currVcf.Pos, currVcf.Ref, currVcf.Alt)
						exception.PanicOnErr(err)
					}
					i.(fileWriter).WriteToFileHandle(out)
					//the third case handles when the alt matches but not the ref, in which case we invert the VCF.
				} else if QuerySeq(ref, currVcf.Chr, int(currVcf.Pos-1), dna.StringToBases(currVcf.Alt[0])) {
					_, err = fmt.Fprintf(un, "Record below was lifted, but the ref and alt alleles are inverted:\n")
					exception.PanicOnErr(err)
					i.WriteToFileHandle(un)
					currVcf = vcf.InvertVcf(currVcf)
					if swapAlleleAB {
						swapInfoAlleles(&currVcf)
					}
					i = &currVcf
					i.(fileWriter).WriteToFileHandle(out)
				} else {
					_, err = fmt.Fprintf(un, "For the following record, neither the Ref nor the Alt allele matched the bases in the corresponding destination fasta location.\n")
					exception.PanicOnErr(err)
					i.WriteToFileHandle(un)
				}
			} else {
				i.(fileWriter).WriteToFileHandle(out)
			}
		}
	}
	return in, mapped
}

func usage() {
	fmt.Print(
		"liftCoordinates - Lifts a compatible file format between assembly coordinates.\n" +
			"Current support for bed and vcf format files.\n" +
			"Usage:\n" +
			"liftCoordinates lift.chain inFile outFile unMapped\n" +
			"Warning: For Vcf lift, the original headers are retained in the output without modification. Use output header information at your own risk.\n" +
			"Please note: Vcf lift is not compatible with Unix piping.\n" +
			"Please note: Vcf lift with fa cross-referencing is currently only supported for biallelic and substitution variants.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var faFile *string = flag.String("faFile", "", "Specify a fasta file for Lifting VCF format files. This fasta should correspond to the destination assembly.")
	var minMatch *float64 = flag.Float64("minMatch", 0.95, "Specify the minimum proportion of matching bases required for a successful lift operation.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to enable debug prints.")
	var swapAlleleAB *bool = flag.Bool("swapAlleleAB", false, "Swap 'Allele_A' and 'Allele_B' designations in the INFO field if ref/alt are inverted during lift.")

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
	liftCoordinates(chainFile, inFile, outFile, *faFile, unMapped, *minMatch, *swapAlleleAB, *verbose)
}
