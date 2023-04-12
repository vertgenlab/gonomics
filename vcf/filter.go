package vcf

import (
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
)

// IsHeterozygous returns true if more than 1 allele is present in the sample
func IsHeterozygous(s Sample) bool {
	if len(s.Alleles) == 0 {
		return false
	}
	alleleBase := s.Alleles[0]
	for i := 1; i < len(s.Alleles); i++ {
		if s.Alleles[i] != alleleBase {
			return true
		}
	}
	return false
}

// IsHomozygous returns true if only 1 allele is present in the sample.
// Note that IsHomozygous also returns true for hemizygous samples.
func IsHomozygous(s Sample) bool {
	if len(s.Alleles) == 0 {
		return false
	}
	alleleBase := s.Alleles[0]
	for i := 1; i < len(s.Alleles); i++ {
		if s.Alleles[i] != alleleBase {
			return false
		}
	}
	return true
}

//IsBiallelic returns true if a vcf record has 1 alt variant, false otherwise.
func IsBiallelic(v Vcf) bool {
	return len(v.Alt) == 1
}

//IsSubstitution returns true if all of the alt fields of a vcf records are of length 1, false otherwise.
func IsSubstitution(v Vcf) bool {
	if len(v.Ref) != 1 {
		return false
	}
	for _, alt := range v.Alt {
		if len(alt) != 1 {
			return false
		}
	}
	return true
}

//IsSegregating returns true if a Vcf record is a segregating site, true if the samples of the record contain at least two allelic states (ex. not all 0 or all 1).
func IsSegregating(v Vcf) bool {
	if len(v.Samples) == 0 {
		return false //special case, no samples
	}

	var firstEncounteredBase int16
	var firstEncountered bool = false

	for _, sample := range v.Samples {
		if len(sample.Alleles) == 0 {
			continue
		}
		if !firstEncountered {
			firstEncounteredBase = sample.Alleles[0]
			firstEncountered = true
		}

		for i := range sample.Alleles {
			if sample.Alleles[i] != firstEncounteredBase {
				return true
			}
		}
	}
	return false
}

//IsPolarizable returns true if a variant can be "polarized" in a derived allele frequency spectrum, false otherwise.
func IsPolarizable(v Vcf) bool {
	if !HasAncestor(v) {
		return false
	}
	var aa string = dna.BasesToString(QueryAncestor(v))
	if len(aa) > 1 || aa == "-" || aa == "N" {
		return false
	}
	if aa != v.Ref && aa != v.Alt[0] { //if ancestral allele is equal to neither the alt or ref allele.
		return false
	}
	return true
}

//IsRefWeakAltStrong returns true if an input biallelic substitution variant has a weak Ref allele and a strong Alt allele, false otherwise.
func IsRefWeakAltStrong(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) {
		return false
	}
	if v.Ref == "A" || v.Ref == "T" {
		if v.Alt[0] == "C" || v.Alt[0] == "G" {
			return true
		}
	}
	return false
}

//IsStrongToWeak returns true if an input biallelic substitution variant has a strong Ref allele and a weak Alt allele, false otherwise.
func IsRefStrongAltWeak(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) {
		return false
	}
	if v.Ref == "C" || v.Ref == "G" {
		if v.Alt[0] == "A" || v.Alt[0] == "T" {
			return true
		}
	}
	return false
}

//IsNotRefStrongAltWeak returns true if an input biallelic substitution variant is not a strong to weak variant, false otherwise.
func IsNotRefStrongAltWeak(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) { //ensures the answer is false if we do not match this initial exclusion criteria.
		return false
	}
	return !IsRefStrongAltWeak(v)
}

//IsNotRefWeakAltStrong returns true if an input biallelic substitution variant does not have a weak Ref allele and a strong Alt allele, false otherwise.
func IsNotRefWeakAltStrong(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) { //ensures the answer is false if we do not match this initial exclusion criteria.
		return false
	}
	return !IsRefWeakAltStrong(v)
}

//IsWeakToStrongOrStrongToWeak returns true if an input biallelic substitution variant is a strong to weak variant or a weak to strong variant, false otherwise.
func IsWeakToStrongOrStrongToWeak(v Vcf) bool {
	return IsRefStrongAltWeak(v) || IsRefWeakAltStrong(v)
}

//IsNotWeakToStrongOrStrongToWeak returns true if a variant is neither a weak to strong variant nor a strong to weak variant, false otherwise.
func IsNotWeakToStrongOrStrongToWeak(v Vcf) bool {
	return IsNotRefWeakAltStrong(v) && IsNotRefStrongAltWeak(v)
}

//SampleVcf takes a set of Vcf records and returns a random subset of variants to an output VCF file. Can also retain a random subset of alleles from gVCF data (diploid, does not break allele pairs)
func SampleVcf(records []Vcf, header Header, numVariants int, numSamples int) ([]Vcf, Header) {
	var sampleList []string
	if len(header.Text) > 0 {
		sampleList = HeaderGetSampleList(header)
	}
	if numVariants > len(records) {
		log.Fatalf("The Number of requested sampled variants is greater than the Number of variants in the input file.")
	}
	//Shuffle the vcf records, our subset will be composed to the first entries in the shuffled order.
	rand.Shuffle(len(records), func(i, j int) { records[i], records[j] = records[j], records[i] })
	records = records[:numVariants] //keep only as many results as specified
	if numSamples > 0 {
		if numSamples > len(records[0].Samples) {
			log.Fatalf("More samples were requested than were present in the input VCF file.")
		}
		var sequentialSlice []int = getSampleKeepList(len(records[0].Samples), numSamples)

		if len(header.Text) > 0 {
			var outHeaderSampleList []string = make([]string, 0)
			for _, i := range sequentialSlice {
				outHeaderSampleList = append(outHeaderSampleList, sampleList[i])
			}
			header = HeaderUpdateSampleList(header, outHeaderSampleList)
		}

		var outSamples []Sample

		for i := range records {
			outSamples = make([]Sample, 0, len(sequentialSlice))
			for _, j := range sequentialSlice {
				outSamples = append(outSamples, records[i].Samples[j])
			}
			records[i].Samples = outSamples
		}
	}
	return records, header
}

func getSampleKeepList(n int, numSamples int) []int {
	var sequentialSlice []int = getSequentialSlice(n)
	rand.Shuffle(len(sequentialSlice), func(i, j int) { sequentialSlice[i], sequentialSlice[j] = sequentialSlice[j], sequentialSlice[i] })
	sequentialSlice = sequentialSlice[:numSamples] //now we have a list of samples to keep from each variant.
	return sequentialSlice
}

//returns a slice where the value is the index. Answer is of length n. ex (4) returns [0 1 2 3]
func getSequentialSlice(n int) []int {
	var answer []int = make([]int, n)
	for i := 0; i < n; i++ {
		answer[i] = i
	}
	return answer
}
