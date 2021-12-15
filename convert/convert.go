//package convert contains functions for converting data between standard file formats. This is a high level package that avoids circular dependencies.

package convert

import (
	//DEBUG: "fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

//singleBedToFasta extracts a sub-Fasta from a reference Fasta sequence at positions specified by an input bed.
func singleBedToFasta(b bed.Bed, ref []fasta.Fasta) fasta.Fasta {
	for i := range ref {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return fasta.Fasta{}
}

//BedToFasta extracts subFastas out of a reference fasta slice comprised of the sequences of input bed regions.
func BedToFasta(b []bed.Bed, ref []fasta.Fasta) []fasta.Fasta {
	outlist := make([]fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = singleBedToFasta(b[i], ref)
	}
	return outlist
}

//SamToBed extracts the position information from a Sam entry and returns it as a bed entry.
func SamToBed(s sam.Sam) bed.Bed {
	if s.Cigar[0].Op == '*' {
		return bed.Bed{}
	} else {
		return bed.Bed{Chrom: s.RName, ChromStart: int(s.Pos - 1), ChromEnd: int(s.Pos-1) + cigar.ReferenceLength(s.Cigar), Name: s.QName, FieldsInitialized: 4}
	}
}

/* TODO: Write Sam to Bed conversion for paired reads.

func SamToBedPaired(s *sam.Sam) []*bed.Bed {
	//sort sam by QName
	//check for "properly aligned" flag
	//grab two properly paired samAln (sanme QName with strings.suffix removed), feed into helper function for bed conversion
	//add output to bedlist
} */

//SamToBedFrag converts a Sam entry into a bed based on the fragment length from which the aligned read was derived. Uses a chromInfo map to ensure fragments are called within the ends of the chromosomes.
func SamToBedFrag(s sam.Sam, fragLength int, reference map[string]chromInfo.ChromInfo) bed.Bed {
	//fatal if fragLength is shorter than sam read length
	if (fragLength <  len(s.Seq)) {
		log.Fatalf("Error: fragLength is %d, which is shorter than the sam read length %d\n", fragLength, len(s.Seq))
	}

	var answer bed.Bed

	if s.Cigar[0].Op == '*' {
		return bed.Bed{}
	} else {
		answer = bed.Bed{Chrom: s.RName, Name: s.QName, FieldsInitialized: 4}
		if sam.IsPosStrand(s) {
			answer.ChromStart = int(s.Pos - 1)
			answer.ChromEnd = numbers.Min(answer.ChromStart+fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar), reference[answer.Chrom].Size)
			answer.Strand = bed.Positive
		} else {
			answer.ChromEnd = int(s.Pos-1) + cigar.ReferenceLength(s.Cigar)
			answer.Strand = bed.Negative
			answer.ChromStart = numbers.Max(answer.ChromEnd-(fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar)), 0)
		}
		return answer
	}
}

//BedNameToWig uses bed entries from an input file to construct a Wig data structure where the Wig value is euqla to the float64-casted nameof an overlapping bed entry. Regions with no bed entries will be set to the value set by Missing (default 0 in the cmd).
func BedValuesToWig(inFile string, reference map[string]chromInfo.ChromInfo, Missing float64, method string) []wig.Wig {
	wigSlice := make([]wig.Wig, len(reference))
	var chromIndex int
	var midpoint int
	var x int
	var i int = 0
	var currentWig wig.Wig
	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig = wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1, Span: -1}
		currentWig.Values = make([]float64, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = Missing
		}
		wigSlice[i] = currentWig
		i++
	}
	bedChan := bed.GoReadToChan(inFile)
	for b := range bedChan {
		chromIndex = getWigChromIndex(b.Chrom, wigSlice)
		midpoint = bedMidpoint(b)
		if wigSlice[chromIndex].Values[midpoint] != Missing {
			log.Fatalf("Two bed entries share the same midpoint. Unable to resolve ambiguous value assignment.")
		}
		if method == "Name" {
			wigSlice[chromIndex].Values[midpoint] = common.StringToFloat64(b.Name)
		} else if method == "Score" {
			wigSlice[chromIndex].Values[midpoint] = float64(b.Score)
		} else {
			log.Fatalf("Unrecognized method.")
		}
	}
	return wigSlice
}

//BedReadsToWig returns a slice of Wig structs where the wig scores correspond to the number of input bed entries that overlap the position.
func BedReadsToWig(b []bed.Bed, reference map[string]chromInfo.ChromInfo) []wig.Wig {
	wigSlice := make([]wig.Wig, len(reference))
	var chromIndex int
	var i, x int = 0, 0
	var currentWig wig.Wig

	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig = wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1, Span: -1}
		currentWig.Values = make([]float64, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = 0
		}
		wigSlice[i] = currentWig
		i++
	}

	for j := range b {
		chromIndex = getWigChromIndex(b[j].Chrom, wigSlice)
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			wigSlice[chromIndex].Values[k]++
		}
	}
	return wigSlice
}

//bedMidpoint returns the midpoint position of an input bed entry.
func bedMidpoint(b bed.Bed) int {
	return (b.ChromEnd + b.ChromStart) / 2
}

//getWigChromIndex searches a wig slice for the wig entry with a particular name and returns the index of that entry in the slice.
func getWigChromIndex(s string, wigSlice []wig.Wig) int {
	for i := range wigSlice {
		if s == wigSlice[i].Chrom {
			return i
		}
	}
	log.Fatalf("Bed Chromosome, %s, not in reference genome.", s)
	return -1
}

//PairwiseFaToVcf takes in a pairwise multiFa alignment and writes Vcf entries for segregating sites with the first entry as the reference and the second fasta entry as the alt allele.
//This will have to be done by chromosome, as a pairwise multiFa will only have two entries, thus containing one chromosome per file.
func PairwiseFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter, substitutionsOnly bool, retainN bool) {
	var pastStart bool = false //bool check to see if we have an insertion at the start of an alignment.
	var insertion bool = false
	var deletion bool = false
	var insertionAlnPos int
	var deletionAlnPos int

	for i := range f[0].Seq { //loop through reference alignment positions
		if f[0].Seq[i] == dna.Gap { //reference is gap (insertion)
			if pastStart {
				if !insertion {
					insertionAlnPos = i - 1
				}
				insertion = true
			}
		} else if f[0].Seq[i] != f[1].Seq[i] {
			pastStart = true
			if insertion { //catches the case where an insertion, now complete, is followed directly by a snp.
				if !substitutionsOnly {
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], insertionAlnPos) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
				}
			}
			if f[1].Seq[i] == dna.Gap { //alt is gap (deletion)
				if !deletion {
					deletionAlnPos = i - 1
				}
				deletion = true
			} else if deletion { //snp immediately follows the end of a deletion
				deletion = false
				if !substitutionsOnly {
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], deletionAlnPos) + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
				}
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
					}
				} else {
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
				}
			} else { //this case is for normal substitutions
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				} else {
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
				}
			}
			insertion = false
		} else if insertion { //case where ref and alt agree now but previous bases were part of an insertion.
			pastStart = true
			insertion = false
			if !substitutionsOnly {
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], insertionAlnPos) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
			}
		} else if deletion {
			pastStart = true
			deletion = false
			if !substitutionsOnly {
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], deletionAlnPos) + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
			}
		}
	}
}
