package pFasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/vcf"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"math/rand"
)

// checks if input pFasta has a sequence with chrom as name and returns its index
func checkIfChromInPfasta(input []PFasta, chrom string) int {
	chromInInput := false
	var answer int

	for inputIdx, inputpFa := range input {
		if inputpFa.Name == chrom {
			chromInInput = true
			answer = inputIdx
		}
	}

	if !chromInInput {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}

	return answer
}

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input []PFasta, start int, end int, outputName string, chrom string, takeCoords bool) PFasta {

	chromIdx := checkIfChromInPfasta(input, chrom)

	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input[chromIdx].Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var outName string
	if takeCoords {
		outName = fmt.Sprintf("%s:%v-%v", chrom, start, end)
	} else if len(outputName) > 0 {
		outName = outputName
	} else {
		outName = chrom
	}

	var answer = PFasta{Name: outName, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input[chromIdx].Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a pFa that has a list of subsequences of the input pFa
// defined by the regions in the bed region
// takeCoords specifies if name fields in output should be original names in region or identified by ChromStart and ChromEnd
func ExtractBed(input []PFasta, region []bed.Bed, takeCoords bool) []PFasta {
	answer := make([]PFasta, 0)
	for _, reg := range region {
		answer = append(answer, Extract(input, reg.ChromStart, reg.ChromEnd, "", reg.Chrom, takeCoords))
	}
	return answer
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution
func Sample(input []PFasta, chrom string) fasta.Fasta {
	chromIdx := checkIfChromInPfasta(input, chrom)

	var answer = fasta.Fasta{Name: input[chromIdx].Name, Seq: make([]dna.Base, len(input[chromIdx].Seq))}
	var currRand float32
	for inputIdx := range input[chromIdx].Seq {
		currRand = rand.Float32()
		if currRand < input[chromIdx].Seq[inputIdx].A {
			answer.Seq[inputIdx] = dna.A
		} else if currRand < (input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.C
		} else if currRand < (input[chromIdx].Seq[inputIdx].G + input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.G
		} else {
			answer.Seq[inputIdx] = dna.T
		}
	}

	return answer
}

// faToPfa returns a pFasta representation of the given Fasta sequence
func faToPfa(input fasta.Fasta, start int, end int) PFasta {
	answer := PFasta{Name: input.Name, Seq: make([]pDna.Float32Base, end-start)}
	if end == -1 {
		end = len(input.Seq)
	}

	fasta.ToUpper(input)

	for idx, base := range input.Seq[start:end] {
		if base == dna.A {
			answer.Seq[idx] = pDna.Float32Base{A: 1, C: 0, G: 0, T: 0}
		} else if base == dna.C {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 1, G: 0, T: 0}
		} else if base == dna.G {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 1, T: 0}
		} else if base == dna.T {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 1}
		} else if base == dna.N {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0.25, T: 0.25}
		} else if base == dna.Gap {
			log.Fatalf("Must specify a sequence without gaps.")
		}
	}

	return answer
}

// faToPfa returns a pFasta representation of the given Fasta sequence
func MultiFaToPfa(inputFaFilename string, start int, end int, chrom string) PFasta {
	inputFa := fasta.Read(inputFaFilename)
	chromInInput := false
	var answer PFasta
	for _, seq := range inputFa {
		if seq.Name == chrom {
			chromInInput = true
			answer = faToPfa(seq, start, end)
		}
	}

	if !chromInInput {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}

	return answer
}

// vcfToPfa returns a pFasta representation of the given VCF sequence, only accepts single sequence Fasta
func VcfToPfa(inVcfFilename string, inputFaFilename string) PFasta {
	var vcfRecords <-chan vcf.Vcf

	inputFa := fasta.Read(inputFaFilename)
	answer := faToPfa(inputFa[0], 0, -1)

	vcfRecords, _ = vcf.GoReadToChan(inVcfFilename)

	for v := range vcfRecords {
		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}

		if inputFa[0].Seq[v.Pos-1] != dna.StringToBase(v.Ref) {
			log.Fatal("Error: base in fasta didn't match ref base from VCF record\n")
		}

		answer.Seq[v.Pos-1] = vcfSampleToPdnaBase(v.Samples, v.Ref, v.Alt[0])
	}

	return answer
}

// vcfSampleToPdnaBase calculates the distribution of samples at a position
func vcfSampleToPdnaBase(samples []vcf.Sample, ref string, alt string) pDna.Float32Base {
	// can i just assume that everything only has 2 alleles and multiple len(samples)*2
	totalSamples := 0
	refCount := 0
	altCount := 1
	for _, s := range samples {
		for _, p := range s.Alleles {
			if p == 0 {
				refCount += 1
			} else if p == 1 {
				altCount += 1
			} else {
				fmt.Print("Second allele.")
			}
			totalSamples += 1
		}
	}

	var answer pDna.Float32Base
	if ref == "A" {
		answer.A = float32(refCount) / float32(totalSamples)
	} else if alt == "A" {
		answer.A = float32(altCount) / float32(totalSamples)
	}

	if ref == "C" {
		answer.C = float32(refCount) / float32(totalSamples)
	} else if alt == "C" {
		answer.C = float32(altCount) / float32(totalSamples)
	}

	if ref == "G" {
		answer.G = float32(refCount) / float32(totalSamples)
	} else if alt == "G" {
		answer.G = float32(altCount) / float32(totalSamples)
	}

	if ref == "T" {
		answer.T = float32(refCount) / float32(totalSamples)
	} else if alt == "T" {
		answer.T = float32(altCount) / float32(totalSamples)
	}

	return answer
}
