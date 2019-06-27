package convert 

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/chromInfo"
	"log"
)

func SingleBedToFasta(b *bed.Bed, ref []*fasta.Fasta) *fasta.Fasta {
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
		outlist[i] = SingleBedToFasta(b[i], ref)
	}
	return outlist
}

func SamToBed(s *sam.Sam) []*bed.Bed {
	outlist := make([]*bed.Bed, len(s.Aln))
	var current *bed.Bed

	for i := 0; i < len(s.Aln); i++ {
		current.Chrom = s.Aln[i].RName
		current.ChromStart = s.Aln[i].Pos - 1
		current.ChromEnd = s.Aln[i].Pos + cigar.ReferenceLength(s.Aln[i].Cigar)
		current.Name = s.Aln[i].QName
		outlist = append(outlist, current)
	}

	return outlist
}

func BedToWig(b []*bed.Bed, reference []*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var chromIndex int

	//generate Wig skeleton from reference
	for i := 0; i < len(reference); i++ {
		currentWig := wig.Wig{StepType: "Fixed", Chrom: reference[i].Name, Start: 0, Step: 1}
		currentWig.Values = make([]*wig.WigValue, len(reference))
		wigSlice = append(wigSlice, &currentWig)
	}

	for j := 0; j < len(b); j++ {
		chromIndex = GetWigChromIndex(b[j].Chrom, wigSlice)
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			wigSlice[chromIndex].Values[k].Value++
		}
	}
	return wigSlice
}

func GetWigChromIndex(s string, wigSlice []*wig.Wig) int {
	for i := 0; i < len(wigSlice); i++ {
		if s == wigSlice[i].Chrom {
			return i
		}
	}
	log.Fatalf("Bed Chromosome, %s, not in reference genome.", s)
	return -1
}
