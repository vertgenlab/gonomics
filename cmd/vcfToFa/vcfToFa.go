// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfToFa(inVcfFilename string, inFaFilename string, outFaFilename string, useAlt bool, multiFaMode bool, multiFaChromName string) {
	var seqsOrdered []fasta.Fasta
	var seqsLookup fasta.FastaMap
	var vcfRecords <-chan vcf.Vcf
	var firstTime bool = true
	var prevPos int              //position of the previous variant. Used to enforce sorted order in multiFaMode.
	var multiFaAppendName string //name of the appended record in multiFaMode.
	var currAlnPos int
	var outIndex int
	var prevRefPos, prevAlnPos int = 0, 0

	// We hang on to the ordered version of the fasta sequences so that we could
	// output the edited sequences (those in the map) in the same order that
	// we originally read them, but we do not have that currently
	seqsOrdered = fasta.Read(inFaFilename)
	if multiFaMode && useAlt {
		var refSeqCopy []dna.Base
		if len(seqsOrdered) == 0 {
			log.Fatalf("Input fasta must contain at least one record in multiFaMode.")
		}
		refSeqCopy = make([]dna.Base, len(seqsOrdered[0].Seq))
		copy(refSeqCopy, seqsOrdered[0].Seq)
		multiFaAppendName = seqsOrdered[0].Name + "alt"
		seqsOrdered = append(seqsOrdered, fasta.Fasta{Name: multiFaAppendName, Seq: refSeqCopy})
		outIndex = len(seqsOrdered) - 1 //index of the output sequence.
	}
	seqsLookup = fasta.ToMap(seqsOrdered)

	// we could do some sort of header check here to make sure the fa matches the vcf
	vcfRecords, _ = vcf.GoReadToChan(inVcfFilename)

	for v := range vcfRecords {
		//first a block of code to enforce sorted order in multiFaMode.
		if firstTime && multiFaMode && v.Chr == multiFaChromName { //only check if in multiFaMode and we have the right chrom
			firstTime = false
			prevPos = v.Pos
		} else if multiFaMode {
			if v.Pos <= prevPos && v.Chr == multiFaChromName {
				log.Fatalf("Input VCF variants must be in sorted order in multiFaMode.")
			}
		}

		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}

		if multiFaMode {
			if v.Chr == multiFaChromName && useAlt { //only consider variants with the correct chrom name.
				currAlnPos = fasta.RefPosToAlnPosCounter(seqsOrdered[outIndex], v.Pos-1, prevRefPos, prevAlnPos)
				if seqsLookup[multiFaAppendName][currAlnPos] != dna.StringToBase(v.Ref) {
					log.Fatalf("Error: base in fasta didn't match ref base from VCF record.\n"+
						"CurrAlnPos: %v. VarPos(1base): %v. SequenceFound: %v. v.Ref: %v.\n", currAlnPos, v.Pos, dna.BaseToString(seqsLookup[multiFaAppendName][currAlnPos]), v.Ref)
				}
				if useAlt {
					seqsLookup[multiFaAppendName][currAlnPos] = dna.StringToBase(v.Alt[0])
				}

				//update the ref/aln pos pair
				prevRefPos = v.Pos - 1 //move to 0 base from 1-based vcf
				prevAlnPos = currAlnPos
			}
		} else {
			// check to make sure the fa matches the vcf before we edit
			if seqsLookup[v.Chr][v.Pos-1] != dna.StringToBase(v.Ref) {
				log.Fatal("Error: base in fasta didn't match ref base from VCF record\n")
			}

			if useAlt {
				seqsLookup[v.Chr][v.Pos-1] = dna.StringToBase(v.Alt[0])
			}
		}
	}

	// Even though we edited sequences in the map, we did not do any insertions, so we
	// know that the location is memory did not change.  If in the future this program
	// does insertions to the starting sequences, we will need to update the slice to
	// point to the sequences in the map
	fasta.Write(outFaFilename, seqsOrdered)
}

func usage() {
	fmt.Print(
		"vcfToFa - Use the variant data in the vcf to edit an input fasta of the reference.\n\n" +
			"Usage:\n" +
			"  vcfToFa [options] input.vcf input.fa output.fa\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var useAlt *bool = flag.Bool("useAlt", false, "Retains only biallelic variants in the output file.")
	var multiFaMode *bool = flag.Bool("multiFaMode", false, "Use a multiFa reference chromosome. Variants will be placed in reference"+
		"positions with respect to gaps. Input vcf must be sorted when using multiFaMode. Appends a new sequence with the updated sequence from the reference sequence,"+
		"defined as the first sequence in the multiFa. Currently only supports useAlt.")
	var multiFaChromName *string = flag.String("multiFaChromName", "", "Identifies the chrom of the multiFa file. Only variants matching this chrom name will be used.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inVcfFilename := flag.Arg(0)
	inFaFilename := flag.Arg(1)
	outFaFilename := flag.Arg(2)

	if !(*useAlt) {
		log.Fatal("Error: you must use at least one option or the fasta file will have no edits\n")
	}

	vcfToFa(inVcfFilename, inFaFilename, outFaFilename, *useAlt, *multiFaMode, *multiFaChromName)
}
