// Package convert contains functions for converting data between standard file formats. This is a high level package that avoids circular dependencies.
package convert

import (
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedGraph"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/wig"
)

// SingleBedToFasta extracts a sub-Fasta from a reference Fasta sequence at positions specified by an input bed.
func SingleBedToFasta(b bed.Bed, ref []fasta.Fasta) fasta.Fasta {
	for i := range ref {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return fasta.Fasta{}
}

// BedToFasta extracts sub-Fastas out of a reference Fasta slice comprised of the sequences of input bed regions.
func BedToFasta(b []bed.Bed, ref []fasta.Fasta) []fasta.Fasta {
	outlist := make([]fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = SingleBedToFasta(b[i], ref)
	}
	return outlist
}

// SamToBed extracts the position information from a Sam entry and returns it as a bed entry.
func SamToBed(s sam.Sam) bed.Bed {
	if s.Cigar[0].Op == '*' {
		return bed.Bed{}
	} else {
		return bed.Bed{Chrom: s.RName, ChromStart: int(s.Pos - 1), ChromEnd: int(s.Pos-1) + cigar.ReferenceLength(s.Cigar), Name: s.QName, FieldsInitialized: 4}
	}
}

// SamToBedWithDeletions extracts the position from a Sam entry and returns it in a slice of bed. If there are deletions in the cigar, it will return multiple bed entries corresponding to each separate alignment segment.
func SamToBedWithDeletions(s sam.Sam) []bed.Bed {
	var outBeds []bed.Bed
	var startPos uint32
	var currPos uint32
	if s.Cigar[0].Op == '*' {
		return []bed.Bed{}
	} else {
		currPos = s.Pos - 1
		startPos = currPos
		for i, v := range s.Cigar {
			if v.Op == 'D' {
				outBeds = append(outBeds, bed.Bed{Chrom: s.RName, ChromStart: int(startPos), ChromEnd: int(currPos), Name: s.QName, FieldsInitialized: 4})
				startPos = currPos + uint32(v.RunLength)
				currPos = startPos
				continue
			} else if cigar.ConsumesReference(v.Op) {
				currPos += uint32(v.RunLength)
			}
			if len(s.Cigar)-1 == i {
				outBeds = append(outBeds, bed.Bed{Chrom: s.RName, ChromStart: int(startPos), ChromEnd: int(currPos), Name: s.QName, FieldsInitialized: 4})
			}
		}
	}
	return outBeds
}

// SamToBedFrag converts a Sam entry into a bed based on the fragment length from which the aligned read was derived.
// Uses a chromInfo map to ensure fragments are called within the ends of the chromosomes.
func SamToBedFrag(s sam.Sam, fragLength int, reference map[string]chromInfo.ChromInfo) bed.Bed {
	//fatal if fragLength is shorter than sam read length
	if fragLength < len(s.Seq) {
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

// BedGraphToWig uses bedGraph entries to construct a slice of Wig data structures where the Wig value is equal to the
// DataValue for the range of the bedGraph entry. Regions with no bedGraph entries will be set to the
// value set by Missing (default 0 in cmd).
func BedGraphToWig(inFile string, reference map[string]chromInfo.ChromInfo, missing float64) map[string]wig.Wig {
	wigMap := wig.MakeSkeleton(reference, missing)
	var i int
	bedChan := bedGraph.GoReadToChan(inFile)
	for b := range bedChan {
		for i = b.ChromStart; i < b.ChromEnd; i++ {
			if wigMap[b.Chrom].Values[i] != missing {
				log.Fatalf("Error in BedGraphToWig. Multiple bed entries map to the same position.")
			}
			wigMap[b.Chrom].Values[i] = b.DataValue
		}
	}
	return wigMap
}

// BedValuesToWig uses bed entries from an input file to construct a Wig data structure where the Wig value is
// equal to the float64-casted name of an overlapping bed entry. Regions with no bed entries will be set to the
// value set by Missing (default 0 in the cmd).
// useRange sets the wig value to the bed value across the range of the bed region, not just at the midpoint, as is default.
func BedValuesToWig(inFile string, reference map[string]chromInfo.ChromInfo, Missing float64, method string, useRange bool, annotationField int) map[string]wig.Wig {
	wigMap := wig.MakeSkeleton(reference, Missing)
	var midpoint, i int
	bedChan := bed.GoReadToChan(inFile)
	for b := range bedChan {
		midpoint = bedMidpoint(b)
		if !useRange && wigMap[b.Chrom].Values[midpoint] != Missing {
			log.Fatalf("Two bed entries share the same midpoint. Unable to resolve ambiguous value assignment.")
		}

		if useRange {
			for i = b.ChromStart; i < b.ChromEnd; i++ {
				if wigMap[b.Chrom].Values[i] != Missing {
					log.Fatalf("Error: overlapping bed elements detected in bed file. Run bedMerge and rerun.")
				}
				if method == "Name" {
					wigMap[b.Chrom].Values[i] = parse.StringToFloat64(b.Name)
				} else if method == "Score" {
					wigMap[b.Chrom].Values[i] = float64(b.Score)
				} else if method == "Annotation" {
					if annotationField >= len(b.Annotation) {
						log.Fatalf("Error: annotationField, %v, exceeds the length of the annotation slice in the following bed entry:\n%v", annotationField, bed.ToString(b, b.FieldsInitialized))
					}
					wigMap[b.Chrom].Values[i] = parse.StringToFloat64(b.Annotation[annotationField])
				} else {
					log.Fatalf("Unrecognized method.")
				}
			}
		} else {
			if method == "Name" {
				wigMap[b.Chrom].Values[midpoint] = parse.StringToFloat64(b.Name)
			} else if method == "Score" {
				wigMap[b.Chrom].Values[midpoint] = float64(b.Score)
			} else if method == "Annotation" {
				if annotationField >= len(b.Annotation) {
					log.Fatalf("Error: annotationField, %v, exceeds the length of the annotation slice in the following bed entry:\n%v", annotationField, bed.ToString(b, b.FieldsInitialized))
				}
				wigMap[b.Chrom].Values[midpoint] = parse.StringToFloat64(b.Annotation[annotationField])
			} else {
				log.Fatalf("Unrecognized method.")
			}
		}
	}
	return wigMap
}

// BedReadsToWig returns a slice of Wig structs where the wig scores correspond to the number of input bed entries that
// overlap the position.
func BedReadsToWig(b []bed.Bed, reference map[string]chromInfo.ChromInfo) map[string]wig.Wig {
	wigMap := wig.MakeSkeleton(reference, 0)
	for j := range b {
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			wigMap[b[j].Chrom].Values[k]++
		}
	}

	return wigMap
}

func BedReadUpdateWig(w map[string]wig.Wig, b bed.Bed) {
	if _, foundInMap := w[b.Chrom]; !foundInMap {
		log.Fatalf("Error: chrom in bed: %s, not found in reference genome.\n", b.Chrom)
	}
	for currPos := b.ChromStart; currPos < b.ChromEnd; currPos++ {
		w[b.Chrom].Values[currPos]++
	}

}

// bedMidpoint returns the midpoint position of an input bed entry.
func bedMidpoint(b bed.Bed) int {
	return (b.ChromEnd + b.ChromStart) / 2
}

// ThreeWayFaToVcf takes in a three-way multiFa alignment and writes Vcf entries for segregating sites with the first
// entry as the reference and the last two fasta entries as the alt alleles.
// This is done by chromosome since a multiFa contains only one chromosome per file.
// This function only checks for substitutions, not indels.
func ThreeWayFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter) {
	var currRefPos, currAlnPos int = 0, 0

	if len(f) != 3 {
		log.Fatalf("ThreeWayFaToVcf expects a fasta input with three entries.")
	}

	for i := range f[0].Seq {
		if f[0].Seq[i] == dna.Gap || f[1].Seq[i] == dna.Gap || f[2].Seq[i] == dna.Gap {
			continue
		}
		if f[0].Seq[i] != f[1].Seq[i] || f[0].Seq[i] != f[2].Seq[i] { // normal substitution
			currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
			currAlnPos = i

			var altAllele []string
			var samples []vcf.Sample

			if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] == f[2].Seq[i] { // 1/0; substitution at allele1 only
				altAllele = []string{dna.BaseToString(f[1].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 0},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] == f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] { // 0/1; substitution at allele2 only
				altAllele = []string{dna.BaseToString(f[2].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{0, 1},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] && f[1].Seq[i] == f[2].Seq[i] { // 1/1; substitution at both allele1 and allele2, same alt
				altAllele = []string{dna.BaseToString(f[1].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 1},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] && f[1].Seq[i] != f[2].Seq[i] { // 1/2; substitution at both allele1 and allele2, different alt
				altAllele = []string{dna.BaseToString(f[1].Seq[i]), dna.BaseToString(f[2].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 2},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			}
			vcf.WriteVcf(out, vcf.Vcf{
				Chr:     chr,
				Pos:     currRefPos + 1, // 1-based indexing for VCF
				Id:      ".",
				Ref:     dna.BaseToString(f[0].Seq[i]), // ref allele in string format
				Alt:     altAllele,
				Qual:    100.0,
				Filter:  "PASS",
				Info:    ".",
				Format:  []string{"GT"},
				Samples: samples,
			})
		}
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

// PairwiseFaToVcf takes in a pairwise multiFa alignment and writes Vcf entries for segregating sites with the first
// entry as the reference and the second fasta entry as the alt allele.
// This will have to be done by chromosome, as a pairwise multiFa will only have two entries, thus containing one chromosome per file.
func PairwiseFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter, substitutionsOnly bool, retainN bool) {
	var pastStart, insertion, deletion bool = false, false, false //first bool checks to see if we have an insertion at the start of an alignment.
	var insertionAlnPos, deletionAlnPos int
	var currRefPos, currAlnPos int = 0, 0 //0 based, like fasta. Add 1 to get vcf pos.
	if len(f) != 2 {
		log.Fatalf("PairwiseFaToVcf expects a fasta input with two entries.")
	}

	for i := range f[0].Seq { //loop through alignment positions
		if f[0].Seq[i] == dna.Gap { //reference is gap (insertion)
			if pastStart {
				if !insertion {
					insertionAlnPos = i - 1
				}
				insertion = true
			}
		} else if f[0].Seq[i] != f[1].Seq[i] { // sequences diff at the same position
			pastStart = true
			if insertion { //catches the case where an insertion, now complete, is followed directly by a snp.
				if !substitutionsOnly {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], insertionAlnPos, currRefPos, currAlnPos)
					currAlnPos = insertionAlnPos //update currAlnPos
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
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
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], deletionAlnPos, currRefPos, currAlnPos)
					currAlnPos = deletionAlnPos
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
				}
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
						currAlnPos = i
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
					}
				} else {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
					currAlnPos = i
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
				}
			} else { //this case is for normal substitutions
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
						currAlnPos = i
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				} else {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
					currAlnPos = i
					if i < len(f[0].Seq)-1 { //if there is a next base to look at
						if f[0].Seq[i+1] != dna.Gap && f[1].Seq[i+1] != dna.Gap { //if neither alt nor ref is a gap in the next pos.
							vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
						} else if substitutionsOnly { //we would also write the subsitution in the case where we don't care if the substitution precedes an INDEL because we are only reporting substitutions
							vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
						}
						// Otherwise we won't report this substitution because it wil be part of the INDEL reporting.
					} else { //for a substitution in the final position, we need not check for subsequent INDELs
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				}
			}
			insertion = false
		} else if insertion { //case where ref and alt agree now but previous bases were part of an insertion.
			pastStart = true
			insertion = false
			if !substitutionsOnly {
				currRefPos = fasta.AlnPosToRefPosCounter(f[0], insertionAlnPos, currRefPos, currAlnPos)
				currAlnPos = insertionAlnPos
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
			}
		} else if deletion {
			pastStart = true
			deletion = false
			if !substitutionsOnly && deletionAlnPos >= 0 {
				//TODO: we do not currently save deletions if they occur at the start of an alignment.
				currRefPos = fasta.AlnPosToRefPosCounter(f[0], deletionAlnPos, currRefPos, currAlnPos)
				currAlnPos = deletionAlnPos
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
			}
		}
	}
}
