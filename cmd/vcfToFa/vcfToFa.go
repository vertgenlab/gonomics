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

// addAnotherCopyOfRef takes an alignment and a name.  The function will copy the reference (first sequence
//  in the alignment and place this copy at the end of the alignment with the name that is given.
func addAnotherCopyOfRef(aln []fasta.Fasta, nameOfCopy string) []fasta.Fasta {
	var refSeqCopy []dna.Base
	if len(aln) == 0 {
		log.Fatalf("Input fasta must contain at least one record in multiFaMode.")
	}
	refSeqCopy = make([]dna.Base, len(aln[0].Seq))
	copy(refSeqCopy, aln[0].Seq)
	aln = append(aln, fasta.Fasta{Name: nameOfCopy, Seq: refSeqCopy})
	return aln
}

// updateSeq takes an alignment, seqsOrdered, and a vcf record, v, and then loops through all the sample sequences to update
//  their sequences in the alignment if they have any alt alleles.  currAlnPos is the index in the alignment that corresponds
//  to the reference position in the vcf record.  The header of the vcf file is needed to loop through the samples.  We also
//  need to know sampleSeqIdxOffset, which is the index of the first sequences in the alignment that corresponds to a sample,
//  instead of an original sequence in the alignment.  Haploid is true if each sample provides only one new sequence for the alignment
func updateSeq(seqsOrdered []fasta.Fasta, currAlnPos int, header vcf.Header, v vcf.Vcf, sampleSeqIdxOffset int, haploid bool) {
	var alleleNumber int16
	var sampleIdx, sampleSeqIdx int
	for _, sampleIdx = range header.Samples {
		if len(v.Samples[sampleIdx].Alleles) != 0 {
			alleleNumber = v.Samples[sampleIdx].Alleles[0]
			if alleleNumber > 0 {
				if haploid {
					sampleSeqIdx = sampleSeqIdxOffset + sampleIdx
				} else {
					sampleSeqIdx = sampleSeqIdxOffset + 2*sampleIdx
				}
				seqsOrdered[sampleSeqIdx].Seq[currAlnPos] = dna.StringToBase(v.Alt[alleleNumber-1])
			}
			if !haploid {
				alleleNumber = v.Samples[sampleIdx].Alleles[1]
				if alleleNumber > 0 {
					sampleSeqIdx = sampleSeqIdxOffset + 2*sampleIdx + 1
					seqsOrdered[sampleSeqIdx].Seq[currAlnPos] = dna.StringToBase(v.Alt[alleleNumber-1])
				}
			}
		}
	}
}

func vcfToMultiFa(inVcfFilename string, inFaFilename string, outFaFilename string, chromName string, useAlt bool, useSamples bool, haploid bool) {
	var seqsOrdered []fasta.Fasta
	var vcfRecords <-chan vcf.Vcf
	var vcfHeader vcf.Header
	var firstTime bool = true
	var prevPos int //position of the previous variant. Used to enforce sorted order in multiFaMode.
	var currAlnPos, outIndex, prevRefPos, prevAlnPos, altSeqIdx, sampleSeqIdxOffset int
	var sampleName string

	if chromName == "" {
		log.Fatalf("Must specify a chrom name when using multiFa mode.")
	}

	vcfRecords, vcfHeader = vcf.GoReadToChan(inVcfFilename)

	seqsOrdered = fasta.Read(inFaFilename)

	if useAlt {
		seqsOrdered = addAnotherCopyOfRef(seqsOrdered, seqsOrdered[0].Name+"alt")
		altSeqIdx = len(seqsOrdered) - 1 //index of the output sequence.
	}
	if useSamples {
		sampleSeqIdxOffset = len(seqsOrdered)
		for sampleName = range vcfHeader.Samples {
			seqsOrdered = addAnotherCopyOfRef(seqsOrdered, sampleName)
			if !haploid {
				seqsOrdered = addAnotherCopyOfRef(seqsOrdered, sampleName+"_secondAllele")
			}
		}
	}

	for v := range vcfRecords {
		//first a block of code to enforce sorted order.
		if firstTime && v.Chr == chromName { //only check if we have the right chrom
			firstTime = false
			prevPos = v.Pos
		} else if v.Pos <= prevPos && v.Chr == chromName {
			log.Fatalf("Input VCF variants must be in sorted order in multiFaMode.")
		}

		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}
		if v.Chr == chromName { //only consider variants with the correct chrom name.
			currAlnPos = fasta.RefPosToAlnPosCounter(seqsOrdered[outIndex], v.Pos-1, prevRefPos, prevAlnPos)
			if seqsOrdered[0].Seq[currAlnPos] != dna.StringToBase(v.Ref) {
				log.Fatalf("Error: base in fasta didn't match ref base from VCF record.\n"+
					"CurrAlnPos: %v. VarPos(1base): %v. SequenceFound: %v. v.Ref: %v.\n", currAlnPos, v.Pos, dna.BaseToString(seqsOrdered[0].Seq[currAlnPos]), v.Ref)
			}
			if useAlt {
				seqsOrdered[altSeqIdx].Seq[currAlnPos] = dna.StringToBase(v.Alt[0])
			}
			if useSamples {
				updateSeq(seqsOrdered, currAlnPos, vcfHeader, v, sampleSeqIdxOffset, haploid)
			}
			//update the ref/aln pos pair
			prevRefPos = v.Pos - 1 //move to 0 base from 1-based vcf
			prevAlnPos = currAlnPos
		}
	}
	fasta.Write(outFaFilename, seqsOrdered)
}

func vcfToFa(inVcfFilename string, inFaFilename string, outFaFilename string, useAlt bool) {
	var seqsOrdered []fasta.Fasta
	var seqsLookup fasta.FastaMap
	var vcfRecords <-chan vcf.Vcf

	// We hang on to the ordered version of the fasta sequences so that we could
	// output the edited sequences (those in the map) in the same order that
	// we originally read them, but we do not have that currently
	seqsOrdered = fasta.Read(inFaFilename)
	seqsLookup = fasta.ToMap(seqsOrdered)

	// we could do some sort of header check here to make sure the fa matches the vcf
	vcfRecords, _ = vcf.GoReadToChan(inVcfFilename)
	for v := range vcfRecords {
		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}

		// check to make sure the fa matches the vcf before we edit
		if seqsLookup[v.Chr][v.Pos-1] != dna.StringToBase(v.Ref) {
			log.Fatal("Error: base in fasta didn't match ref base from VCF record\n")
		}

		if useAlt {
			seqsLookup[v.Chr][v.Pos-1] = dna.StringToBase(v.Alt[0])
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
	var useAlt *bool = flag.Bool("useAlt", false, "Creates a new sequence using the alternative base at every site (only uses biallelic SNPs).")
	var useSamples *bool = flag.Bool("useSamples", false, "Creates new sequences using the genotype of that individual (only uses biallic SNPs).")
	var multiFaMode *bool = flag.Bool("multiFaMode", false, "Use a multiFa reference chromosome. Variants will be placed in reference"+
		"positions with respect to gaps. Input vcf must be sorted when using multiFaMode. Appends a new sequence with the updated sequence from the reference sequence,"+
		"defined as the first sequence in the multiFa. Currently only supports useAlt.")
	var multiFaChromName *string = flag.String("multiFaChromName", "", "Identifies the chrom of the multiFa file. Only variants matching this chrom name will be used.")
	var haploid *bool = flag.Bool("haploid", false, "Only make one sequence per individual in the vcf, instead of two.  The first genotype will be used.")

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

	if !(*useAlt || *useSamples) {
		log.Fatal("Error: you must use at least one option or the fasta file will have no edits\n")
	}

	if *multiFaMode {
		vcfToMultiFa(inVcfFilename, inFaFilename, outFaFilename, *multiFaChromName, *useAlt, *useSamples, *haploid)
	} else {
		vcfToFa(inVcfFilename, inFaFilename, outFaFilename, *useAlt)
	}
}
