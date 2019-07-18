package convert 

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/chromInfo"
	"log"
	"fmt"
)

func singleBedToFasta(b *bed.Bed, ref []*fasta.Fasta) *fasta.Fasta {
	for i := 0; i < len(ref); i++ {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return nil
}

func BedToFasta(b []*bed.Bed, ref []*fasta.Fasta) []*fasta.Fasta {
	outlist := make([]*fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = singleBedToFasta(b[i], ref)
	}
	return outlist
}

func SamToBed(s *sam.Sam) []*bed.Bed {
	var outlist []*bed.Bed
	var current *bed.Bed

	for i := 0; i < len(s.Aln); i++ {
		if !(s.Aln[i].Cigar[0].Op == '*') {
			current = &bed.Bed{Chrom: s.Aln[i].RName,ChromStart: s.Aln[i].Pos - 1, ChromEnd: s.Aln[i].Pos + cigar.ReferenceLength(s.Aln[i].Cigar), Name: s.Aln[i].QName}
			outlist = append(outlist, current)
		}
	}

	return outlist
}

/* TODO: Write Sam to Bed conversion for paired reads.

func SamToBedPaired(s *sam.Sam) []*bed.Bed {
	//sort sam by QName
	//check for "properly aligned" flag
	//grab two properly paired samAln (sanme QName with strings.suffix removed), feed into helper function for bed conversion
	//add output to bedlist
} */

func SamToBedFrag(s *sam.Sam, fragLength int64, reference map[string]*chromInfo.ChromInfo) []*bed.Bed {
	var outlist []*bed.Bed
	var current *bed.Bed

	for i := 0; i < len(s.Aln); i++ {		
		if !(s.Aln[i].Cigar[0].Op == '*') {
			current = &bed.Bed{Chrom: s.Aln[i].RName, Name: s.Aln[i].QName}
			if sam.IsPosStrand(s.Aln[i]) {
				current.ChromStart = s.Aln[i].Pos - 1
				current.ChromEnd = common.MinInt64(current.ChromStart + fragLength - cigar.NumInsertions(s.Aln[i].Cigar) + cigar.NumDeletions(s.Aln[i].Cigar), reference[current.Chrom].Size)
				current.Strand = true
			} else {
				current.ChromEnd = s.Aln[i].Pos - 1 + cigar.ReferenceLength(s.Aln[i].Cigar)
				current.Strand = false
				current.ChromStart = common.MaxInt64(current.ChromEnd - (fragLength - cigar.NumInsertions(s.Aln[i].Cigar) + cigar.NumDeletions(s.Aln[i].Cigar)), 0)
			}
			outlist = append(outlist, current)
		}
	}

	return outlist
}

func BedScoreToWig(b []*bed.Bed, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var chromIndex int
	var midpoint int
	var x int64 = 0
	var i int = 0

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
	fmt.Println("Completed wig skeleton, looping through bed.")

	//loop through bed
	for i := 0; i < len(b); i++ {
		chromIndex = getWigChromIndex(b[i].Chrom, wigSlice)
		midpoint = bedMidpoint(b[i])
		if wigSlice[chromIndex].Values[midpoint].Value != 0 {
			log.Fatalf("Multiple scores for one position.")
		}
		wigSlice[chromIndex].Values[midpoint].Value = float64(b[i].Score)
	}
	return wigSlice
}


func BedReadsToWig(b []*bed.Bed, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var chromIndex int
	var i int = 0
	var x int64 = 0
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
		fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			//fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, k: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, k, len(wigSlice[chromIndex].Values))
			wigSlice[chromIndex].Values[k].Value++
		}
		//fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
	}
	return wigSlice
}

func bedMidpoint(b *bed.Bed) int {
	return int(b.ChromEnd + b.ChromStart) / 2
}

func getWigChromIndex(s string, wigSlice []*wig.Wig) int {
	for i := 0; i < len(wigSlice); i++ {
		if s == wigSlice[i].Chrom {
			return i
		}
	}
	log.Fatalf("Bed Chromosome, %s, not in reference genome.", s)
	return -1
}
