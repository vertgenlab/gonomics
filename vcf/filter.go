package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"math/rand"
	"strings"
)

//Filter returns true if a Vcf passes a set of filter criteria, false otherwise. Special empty strings "" for alt and ref automatically pass.
func Filter(v Vcf, chrom string, minPos int, maxPos int, ref string, alt []string, minQual float64, biAllelicOnly bool, substitutionsOnly bool, segregatingSitesOnly bool, removeNoAncestor bool, onlyPolarizableAncestors bool) bool {
	if !FilterRange(v, minPos, maxPos) {
		return false
	}
	if !FilterChrom(v, chrom) {
		return false
	}
	if len(ref) > 0 && !FilterRef(v, ref) {
		return false
	}
	if len(alt) > 0 && len(alt[0]) > 0 && !FilterAlt(v, alt) {
		return false
	}
	if !FilterQual(v, minQual) {
		return false
	}
	if biAllelicOnly && !IsBiallelic(v) {
		return false
	}
	if substitutionsOnly && !IsSubstitution(v) {
		return false
	}
	if segregatingSitesOnly && !IsSegregating(v) {
		return false
	}
	if removeNoAncestor && !HasAncestor(v) {
		return false
	}
	if onlyPolarizableAncestors && !IsPolarizable(v) {
		return false
	}
	return true
}

//FilterQual returns true if a Vcf quality is above an input minQual score, false otherwise.
func FilterQual(v Vcf, minQual float64) bool {
	if v.Qual < minQual {
		return false
	}
	return true
}

//FilterAlt returns true if the Alt field of a Vcf entry matches a desired input Alt field, false otherwise. Order sensitive.
func FilterAlt(v Vcf, alt []string) bool {
	if len(alt) > 0 && CompareAlt(v.Alt, alt) != 0 {
		return false
	}
	return true
}

//FilterRef returns true if the Ref field of a Vcf record matches an input string, false otherwise.
func FilterRef(v Vcf, ref string) bool {
	if strings.Compare(v.Ref, ref) != 0 {
		return false
	}
	return true
}

//FilterRange returns true if a Vcf position lies between an input minimum and maximum position (inclusive for both min and max), false otherwise.
func FilterRange(v Vcf, minPos int, maxPos int) bool {
	if v.Pos < minPos || v.Pos > maxPos {
		return false
	}
	return true
}

//FilterChrom returns true if the Chrom field of a Vcf record matches an input string, false otherwise.
func FilterChrom(v Vcf, chrom string) bool {
	if chrom != "" && v.Chr != chrom {
		return false
	}
	return true
}

//FilterNs removes all records from a slice of Vcfs that contain Ns.
func FilterNs(vcfs []Vcf) []Vcf {
	var answer []Vcf
	var noN bool
	for i := 0; i < len(vcfs); i++ {
		noN = true
		if dna.CountBase(dna.StringToBases(vcfs[i].Ref), dna.N) > 0 {
			noN = false
		}
		for j := 0; j < len(vcfs[i].Alt); j++ {
			if dna.CountBase(dna.StringToBases(vcfs[i].Alt[j]), dna.N) > 0 {
				noN = false
			}
		}
		if noN {
			answer = append(answer, vcfs[i])
		}
	}
	return answer
}

func ASFilter(v Vcf, parentOne int16, parentTwo int16, F1 int16) bool {
	if IsHomozygous(v.Samples[parentOne]) && IsHomozygous(v.Samples[parentTwo]) && IsHeterozygous(v.Samples[F1]) && v.Samples[parentOne].AlleleOne != v.Samples[parentTwo].AlleleOne {
		return true
	} else {
		return false
	}
}

func mergeSimilarVcf(a Vcf, b Vcf) Vcf {
	mergeRecord := Vcf{Chr: a.Chr, Pos: a.Pos, Id: a.Id, Ref: "", Qual: a.Qual, Filter: "Merged:SNP:INDEL", Info: a.Info, Format: a.Format, Samples: a.Samples}
	if len(a.Ref) < len(b.Ref) {
		mergeRecord.Ref += b.Ref
	} else {
		mergeRecord.Ref += a.Ref
	}
	if len(a.Alt) < len(b.Alt) {
		mergeRecord.Alt = b.Alt
	} else {
		mergeRecord.Alt = a.Alt
	}
	return mergeRecord
}

func LowQual(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return true
	} else {
		return false
	}
}

func IsHeterozygous(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return false
	}
	if genome.AlleleOne != genome.AlleleTwo {
		return true
	}
	if genome.AlleleOne == genome.AlleleTwo {
		return false
	}
	return false
}

func IsHomozygous(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return false
	}
	if genome.AlleleOne == genome.AlleleTwo {
		return true
	}
	if genome.AlleleOne != genome.AlleleTwo {
		return false
	}
	return false
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
	var firstEncountered int16 = v.Samples[0].AlleleOne

	for _, sample := range v.Samples {
		if sample.AlleleOne != firstEncountered || sample.AlleleTwo != firstEncountered {
			return true
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

//IsWeakToStrong returns true if an input biallelic substitution variant is a weak to strong variant, false otherwise.
func IsWeakToStrong(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) {
		return false
	}
	if v.Ref == "A" || v.Ref == "T" {
		if v.Alt[0] == "C" || v.Alt[0] == "G" {
			return true
		}
	}
	if v.Ref == "C" || v.Ref == "G" {
		if v.Alt[0] == "A" || v.Alt[0] == "T" {
			return true
		}
	}
	return false
}

//IsNotWeakToStrong returns true if an input biallelic substitution variant is not a weak to strong variant, false otherwise.
func IsNotWeakToStrong(v Vcf) bool {
	if !IsBiallelic(v) || !IsSubstitution(v) { //ensures the answer is false if we do not match this initial exclusion criteria.
		return false
	}
	return !IsWeakToStrong(v)
}

func getListIndex(header Header, list []string) []int16 {
	sampleHash := HeaderToMaps(header)
	var listIndex []int16 = make([]int16, len(list))
	for i := 0; i < len(listIndex); i++ {
		//look up alt allele index belonging to each string
		listIndex[i] = sampleHash.GIndex[list[i]]
	}
	return listIndex
}

func ByNames(inChan <-chan Vcf, header Header, list []string, writer *fileio.EasyWriter) {
	var listIndex []int16 = getListIndex(header, list)

	for record := range inChan {
		WriteVcf(writer, ReorderSampleColumns(record, listIndex))
	}
}

//FilterVcfPos will filter out records that appear as the same postion more than once, keeping the first one it encounters. In addition, if records contains Ns, those records will also be filtered out.
func FilterVcfPos(vcfs []Vcf) []Vcf {
	Sort(vcfs)
	var answer []Vcf
	chrVcfMap := make(map[string][]Vcf)

	var ref []dna.Base
	var alt [][]dna.Base
	var containsN bool
	for _, v := range vcfs {
		chrVcfMap[v.Chr] = append(chrVcfMap[v.Chr], v)
	}
	var i int
	var curr []Vcf
	for key := range chrVcfMap {
		curr = chrVcfMap[key]
		encountered := make(map[int]bool)
		for i = 0; i < len(curr); i++ {
			if encountered[curr[i].Pos] == true {
				//do not add
			} else {
				encountered[curr[i].Pos] = true
				containsN = false
				ref = dna.StringToBases(curr[i].Ref)
				alt = GetAltBases(curr[i].Alt)
				if dna.CountBaseInterval(ref, dna.N, 0, len(ref)) != 0 {
					containsN = true
				}
				for j := 0; j < len(alt); j++ {
					if dna.CountBaseInterval(alt[j], dna.N, 0, len(alt[j])) != 0 {
						containsN = true
					}
				}
				if !containsN {
					answer = append(answer, curr[i])
				}
			}
		}
	}
	return answer
}

//SampleVcf takes a VCF file and returns a random subset of variants to an output VCF file. Can also retain a random subset of alleles from gVCF data (diploid, does not break allele pairs)
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
	//DEBUG:fmt.Printf("lenRecords before slice: %v.\n", len(records))
	records = records[:numVariants] //keep only as many results as specified.
	//DEBUG: fmt.Printf("lenRecords after slice: %v.\n", len(records))

	if numSamples > 0 {
		if numSamples > len(records[0].Samples) {
			log.Fatalf("More samples were requested than were present in the input VCF file.")
		}
		var sequentialSlice []int = getSequentialSlice(len(records[0].Samples))
		rand.Shuffle(len(sequentialSlice), func(i, j int) { sequentialSlice[i], sequentialSlice[j] = sequentialSlice[j], sequentialSlice[i] })
		sequentialSlice = sequentialSlice[:numSamples] //now we have a list of samples to keep from each variant.

		if len(header.Text) > 0 {
			var outHeaderSampleList []string = make([]string, 0)
			for _, i := range sequentialSlice {
				outHeaderSampleList = append(outHeaderSampleList, sampleList[i])
			}

			header = HeaderUpdateSampleList(header, outHeaderSampleList)
		}

		var outSamples []GenomeSample

		for i := range records {
			outSamples = make([]GenomeSample, 0, len(sequentialSlice))
			for _, j := range sequentialSlice {
				outSamples = append(outSamples, records[i].Samples[j])
			}
			records[i].Samples = outSamples
		}
	}
	return records, header
}

//returns a slice where the value is the index. Answer is of length n. ex (4) returns [0 1 2 3]
func getSequentialSlice(n int) []int {
	var answer []int = make([]int, n)
	for i := 0; i < n; i++ {
		answer[i] = i
	}
	return answer
}
