package pFasta

import (
	"log"
	"math/rand"
	"reflect"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
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

// faToPfa returns a pFasta representation of the given Fasta sequence, start inclusive, end exclusive
func faToPfa(input fasta.Fasta, start int, end int) PFasta {
	if end == -1 {
		end = len(input.Seq)
	} else if end > len(input.Seq) {
		log.Fatalf("Requested end argument (%v) out of range.", end)
	}

	answer := PFasta{Name: input.Name, Seq: make([]pDna.Float32Base, end-start)}

	fasta.ToUpper(input)
	var base dna.Base
	var idx int
	for idx, base = range input.Seq[start:end] {
		if base == dna.A {
			answer.Seq[idx] = pDna.Float32Base{A: 1, C: 0, G: 0, T: 0}
		} else if base == dna.C {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 1, G: 0, T: 0}
		} else if base == dna.G {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 1, T: 0}
		} else if base == dna.T {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 1}
		} else if base == dna.N {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}
		} else if base == dna.Gap {
			log.Fatalf("Must specify a sequence without gaps.")
		}
	}

	return answer
}

// MultiFaToPfa returns a pFasta representation of the given Fasta sequence
func MultiFaToPfa(inputFaFilename string, start int, end int, chrom string) PFasta {
	inputFa := fasta.Read(inputFaFilename)
	chromInInput := false
	var answer PFasta

	if len(inputFa) == 1 {
		if chrom == "" || inputFa[0].Name == chrom {
			chromInInput = true
			answer = faToPfa(inputFa[0], start, end)
		}
	} else {
		if chrom == "" {
			log.Fatalf("Error: expecting a Chrom argument for multifasta input.")
		}

		for _, seq := range inputFa {
			if seq.Name == chrom {
				chromInInput = true
				answer = faToPfa(seq, start, end)
				break
			}
		}
	}

	if chromInInput == false {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}

	return answer
}

// vcfToPfa returns a pFasta representation of the given VCF sequence, only accepts single sequence Fasta
func VcfToPfa(inVcfFilename string, inputFaFilename string, start int, end int) PFasta {
	// relax to not-biallelic
	var vcfRecords <-chan vcf.Vcf

	inputFa := fasta.Read(inputFaFilename)
	answer := faToPfa(inputFa[0], start, end)

	vcfRecords, _ = vcf.GoReadToChan(inVcfFilename)

	for v := range vcfRecords {
		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}

		if inputFa[0].Seq[v.Pos-1] != dna.StringToBase(v.Ref) {
			log.Fatal("Error: base in fasta didn't match ref base from VCF record\n")
		}

		answer.Seq[v.Pos-1] = vcfSampleToPdnaBase(v.Samples, v.Ref, v.Alt)
	}

	return answer
}

type alleleCounts struct {
	A int
	C int
	G int
	T int
}

func getFieldPointer(counts interface{}, fieldName string) *int {
	v := reflect.ValueOf(counts).Elem()
	field := v.FieldByName(fieldName)

	if field.IsValid() && field.CanAddr() {
		fieldPtr := field.Addr().Interface().(*int)
		return fieldPtr
	} else {
		fmt.Printf("Invalid field: %s\n", fieldName)
		return nil
	}
}

// vcfSampleToPdnaBase calculates the distribution of samples at a position
func vcfSampleToPdnaBase(samples []vcf.Sample, ref string, alts []string) pDna.Float32Base {
	// can i just assume that everything only has 2 alleles and multiple len(samples)*2
	// try to map it? idx in alts list = what base
	// struct {A: 0, C: 0, G: 0, T: 0}
	totalSamples := 2 * len(samples)

	var counts alleleCounts
	var mapping []*int

	// map ref and alt alleles to counts
	fieldPointer := getFieldPointer(&counts, ref)
	if fieldPointer != nil {
		mapping = append(mapping, fieldPointer)
	}

	for _, alt := range alts {
		fieldPointer = getFieldPointer(&counts, alt)
		if fieldPointer != nil {
			mapping = append(mapping, fieldPointer)
		}
	}

	for _, s := range samples {
		for _, p := range s.Alleles {
			if p == 0 {
				*mapping[0] += 1
			} else if p == 1 {
				*mapping[1] += 1
			} else if p == 2 {
				*mapping[2] += 1
			} else if p == 3 {
				*mapping[3] += 1
			} else {
				log.Fatalf("Invalid allele value > 4.")
			}
			totalSamples += 1
		}
	}

	var answer pDna.Float32Base
	answer.A = float32(counts.A) / float32(totalSamples)
	answer.C = float32(counts.C) / float32(totalSamples)
	answer.G = float32(counts.G) / float32(totalSamples)
	answer.T = float32(counts.T) / float32(totalSamples)

	return answer
}
