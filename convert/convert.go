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
	"strings"
)

//singleBedToFasta extracts a sub-Fasta from a reference Fasta sequence at positions specified by an input bed.
func singleBedToFasta(b *bed.Bed, ref []fasta.Fasta) fasta.Fasta {
	for i := 0; i < len(ref); i++ {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return fasta.Fasta{}
}

//BedToFasta extracts subFastas out of a reference fasta slice comprised of the sequences of input bed regions.
func BedToFasta(b []*bed.Bed, ref []fasta.Fasta) []fasta.Fasta {
	outlist := make([]fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = singleBedToFasta(b[i], ref)
	}
	return outlist
}

//SamToBed extracts the position information from a SamAln entry and returns it as a bed entry.
func SamToBed(s *sam.SamAln) *bed.Bed {
	if s.Cigar[0].Op == '*' {
		return nil
	} else {
		return &bed.Bed{Chrom: s.RName, ChromStart: int(s.Pos - 1), ChromEnd: int(s.Pos-1) + cigar.ReferenceLength(s.Cigar), Name: s.QName}
	}
}

/* TODO: Write Sam to Bed conversion for paired reads.

func SamToBedPaired(s *sam.Sam) []*bed.Bed {
	//sort sam by QName
	//check for "properly aligned" flag
	//grab two properly paired samAln (sanme QName with strings.suffix removed), feed into helper function for bed conversion
	//add output to bedlist
} */

//SamToBedFrag converts a SamAln entry into a bed based on the fragment length from which the aligned read was derived. Uses a chromInfo map to ensure fragments are called within the ends of the chromosomes.
func SamToBedFrag(s *sam.SamAln, fragLength int, reference map[string]*chromInfo.ChromInfo) *bed.Bed {
	var answer *bed.Bed

	if s.Cigar[0].Op == '*' {
		return nil
	} else {
		answer = &bed.Bed{Chrom: s.RName, Name: s.QName}
		if sam.IsPosStrand(s) {
			answer.ChromStart = int(s.Pos - 1)
			answer.ChromEnd = numbers.Min(answer.ChromStart+fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar), reference[answer.Chrom].Size)
			answer.Strand = true
		} else {
			answer.ChromEnd = int(s.Pos-1) + cigar.ReferenceLength(s.Cigar)
			answer.Strand = false
			answer.ChromStart = numbers.Max(answer.ChromEnd-(fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar)), 0)
		}
		return answer
	}
}

//BedScoreToWig uses bed entries from an input file to construct a Wig data structure where the Wig value is equal to the score of an overlapping bed entry at the bed entry midpoint, and zero if no bed regions overlap.
func BedScoreToWig(infile string, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var line string
	var chromIndex int
	var midpoint int
	var startNum, endNum int
	var x int
	var i int = 0
	var doneReading bool = false
	var current *bed.Bed

	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	//DEBUG: log.Println("Completed wig skeleton, looping through bed.")

	//loop through bed line at a time
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		startNum = common.StringToInt(words[1])
		endNum = common.StringToInt(words[2])
		current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
		if len(words) >= 4 {
			current.Name = words[3]
		}
		if len(words) >= 5 {
			current.Score = common.StringToInt(words[4])
		}
		chromIndex = getWigChromIndex(current.Chrom, wigSlice)
		midpoint = bedMidpoint(current)
		if wigSlice[chromIndex].Values[midpoint].Value != 0 {
			log.Fatalf("Multiple scores for one position.")
		}

		wigSlice[chromIndex].Values[midpoint].Value = float64(current.Score)

	}
	return wigSlice
}

//BedScoreToWigRange uses bed entries from an input file to construct a Wig data structure where the Wig value at each position is equal to the score of an overlapping bed entry, and zero if no bed regions overlap.
func BedScoreToWigRange(infile string, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var line string
	var chromIndex int
	var midpoint int
	var startNum, endNum int
	var x int
	var i int = 0
	var doneReading bool = false
	var current *bed.Bed

	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	//DEBUG: log.Println("Completed wig skeleton, looping through bed.")

	//loop through bed line at a time
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		startNum = common.StringToInt(words[1])
		endNum = common.StringToInt(words[2])
		current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
		if len(words) >= 4 {
			current.Name = words[3]
		}
		if len(words) >= 5 {
			current.Score = common.StringToInt(words[4])
		}
		chromIndex = getWigChromIndex(current.Chrom, wigSlice)
		if wigSlice[chromIndex].Values[midpoint].Value != 0 {
			log.Fatalf("Multiple scores for one position.")
		}
		for k := current.ChromStart; k < current.ChromEnd; k++ {
			//DEBUG: fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, k: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, k, len(wigSlice[chromIndex].Values))
			wigSlice[chromIndex].Values[k+1].Value = float64(current.Score)
		}
	}
	return wigSlice
}

//BedReadsToWig returns a slice of Wig structs where the wig scores correspond to the number of input bed entries that overlap the position.
func BedReadsToWig(b []*bed.Bed, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var chromIndex int
	var i, x int = 0, 0
	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	for j := 0; j < len(b); j++ {
		chromIndex = getWigChromIndex(b[j].Chrom, wigSlice)
		//DEBUG: fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			//DEBUG: fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, k: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, k, len(wigSlice[chromIndex].Values))
			wigSlice[chromIndex].Values[k].Value++
		}
		//DEBUG: fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
	}
	return wigSlice
}

//bedMidpoint returns the midpoint position of an input bed entry.
func bedMidpoint(b *bed.Bed) int {
	return int(b.ChromEnd+b.ChromStart) / 2
}

//getWigChromIndex searches a wig slice for the wig entry with a particular name and returns the index of that entry in the slice.
func getWigChromIndex(s string, wigSlice []*wig.Wig) int {
	for i := 0; i < len(wigSlice); i++ {
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

	for i := 0; i < len(f[0].Seq); i++ { //loop through reference alignment positions
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
					vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], insertionAlnPos) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
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
					vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], deletionAlnPos) + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
				}
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
					}
				} else {
					vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
				}
			} else { //this case is for normal substitutions
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				} else {
					vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], i) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
				}
			}
			insertion = false
		} else if insertion { //case where ref and alt agree now but previous bases were part of an insertion.
			pastStart = true
			insertion = false
			if !substitutionsOnly {
				vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], insertionAlnPos) + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
			}
		} else if deletion {
			pastStart = true
			deletion = false
			if !substitutionsOnly {
				vcf.WriteVcf(out, &vcf.Vcf{Chr: chr, Pos: fasta.AlnPosToRefPos(f[0], deletionAlnPos) + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
			}
		}
	}
}
