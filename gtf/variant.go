package gtf

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"math"
)

type Variant struct {
	vcf.Vcf
	RefId       string // e.g. NC_000023.10, LRG_199, NG_012232.1, NM_004006.2, LRG-199t1, NR_002196.1, NP_003997.1, etc.
	Gene        string
	PosStrand 	bool
	NearestCDS  *CDS
	CDNAPos     int // 1-base
	AAPos       int // 1-base
	AARef       []dna.AminoAcid
	AAAlt       []dna.AminoAcid
	VariantType string // e.g. Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away), FarSplice (3-10 away)
}

func GenesToIntervalTree(genes map[string]*Gene) map[string]*interval.IntervalNode {
	MoveAllCanonicalToZero(genes)
	intervals := make([]interval.Interval, len(genes))
	var i int = 0
	for k, g := range genes {
		intervals[i] = g
		delete(genes, k)
		i++
	}
	return interval.BuildTree(intervals)
}

func VcfToVariant(v *vcf.Vcf, tree map[string]*interval.IntervalNode, seq map[string][]dna.Base) (*Variant, error) {
	var answer *Variant
	var err error

	overlappingGenes := interval.Query(tree, v, "any")

	if len(overlappingGenes) > 1 {
		err = errors.New(fmt.Sprintf("Variant overlaps with multiple genes. Mutation will be based on first gene."))
	}

	var annotatingGene *Gene
	if len(overlappingGenes) > 0 {
		annotatingGene = overlappingGenes[0].(*Gene)
		answer = vcfToVariant(v, annotatingGene, seq)
	} else {
		answer = &Variant{Vcf: *v}
	}
	addVariantType(answer)
	return answer, err
}

func vcfToVariant(v *vcf.Vcf, gene *Gene, seq map[string][]dna.Base) *Variant {
	answer := new(Variant)
	answer.Vcf = *v
	answer.RefId = gene.Transcripts[0].TranscriptID
	answer.Gene = gene.GeneID
	answer.PosStrand = gene.Transcripts[0].Strand
	vcfCdsIntersect(v, gene, answer)
	if int(v.Pos) >= answer.NearestCDS.Start && int(v.Pos) <= answer.NearestCDS.End {
		findAAChange(answer, seq)
	}
	return answer
}

func vcfCdsIntersect(v *vcf.Vcf, gene *Gene, answer *Variant) {
	var cdsPos int
	var exon *Exon
	//TODO: this code may be able to be compressed
	if answer.PosStrand {
		for i := 0; i < len(gene.Transcripts[0].Exons); i++ {
			exon = gene.Transcripts[0].Exons[i]
			if exon.Cds != nil && int(v.Pos) > exon.Cds.End { // variant is further in gene
				cdsPos += exon.Cds.End - exon.Cds.Start + 1
				answer.NearestCDS = exon.Cds // Store most recent exon and move on // Catches variants past the last exon
			} else if exon.Cds != nil && int(v.Pos) <= exon.Cds.End { // variant is before end of this exon
				if int(v.Pos) < exon.Cds.Start { // Variant is NOT in CDS
					if exon.Cds.Prev == nil || exon.Cds.Start-int(v.Pos) < int(v.Pos)-gene.Transcripts[0].Exons[i-1].Cds.Start {
						answer.NearestCDS = exon.Cds
					} else {
						answer.NearestCDS = gene.Transcripts[0].Exons[i-1].Cds
					}
					break
				}
				cdsPos += int(v.Pos) - exon.Cds.Start + 1
				answer.CDNAPos = cdsPos
				answer.NearestCDS = exon.Cds
			}
		}
	} else {
		for i := 0; i < len(gene.Transcripts[0].Exons); i++ {
			exon = gene.Transcripts[0].Exons[len(gene.Transcripts[0].Exons) - 1 - i]
			if exon.Cds != nil && int(v.Pos) < exon.Cds.Start { // variant is further in gene
				cdsPos += exon.Cds.End - exon.Cds.Start + 1
				answer.NearestCDS = exon.Cds // Store most recent exon and move on // Catches variants past the last exon
			} else if exon.Cds != nil && int(v.Pos) >= exon.Cds.Start { // variant is before end of this exon
				if int(v.Pos) > exon.Cds.End { // Variant is NOT in CDS
					if exon.Cds.Next == nil || int(v.Pos)-exon.Cds.End < gene.Transcripts[0].Exons[len(gene.Transcripts[0].Exons) - 1 - i + 1].Cds.Start-int(v.Pos) {
						answer.NearestCDS = exon.Cds
					} else {
						answer.NearestCDS = gene.Transcripts[0].Exons[len(gene.Transcripts[0].Exons) - 1 - i + 1].Cds
					}
					break
				}
				// Variant IS in CDS
				cdsPos += exon.Cds.End - int(v.Pos) + 1
				answer.CDNAPos = cdsPos
				answer.NearestCDS = exon.Cds
			}
		}
	}

}

func findAAChange(variant *Variant, seq map[string][]dna.Base) {
	var refBases = make([]dna.Base, 0)
	var altBases = make([]dna.Base, 0)
	var seqPos int = int(variant.Pos) - 1
	var currCDS *CDS = variant.NearestCDS
	var aaPosOffset int = 0
	if variant.PosStrand {
		seqPos -= determineFrame(variant)
		for ; seqPos < int(variant.Pos-1); seqPos++ {
			if seqPos < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			} else if seqPos > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			refBases = append(refBases, seq[variant.Chr][seqPos])
			altBases = append(altBases, seq[variant.Chr][seqPos])
		}
		refBases = append(refBases, dna.StringToBases(variant.Ref)...)
		altBases = append(altBases, dna.StringToBases(variant.Alt)...)
		seqPos += len(dna.StringToBases(variant.Ref))
		var offset int
		for offset = 0; len(altBases)%3 != 0; offset++ {
			if seqPos+offset > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			altBases = append(altBases, seq[variant.Chr][seqPos+offset])
		}
		for offset = 0; len(refBases)%3 != 0; offset++ {
			if seqPos+offset > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			refBases = append(refBases, seq[variant.Chr][seqPos+offset])
		}
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)

		//TODO: WARNING THIS CODE MAY BE BUGGY
		if !isSynonymous(variant) && len(variant.AAAlt) > 1 {
			var codonToAdd []dna.Base
			var j int
			for len(variant.AAAlt) > 0 && variant.AARef[0] == variant.AAAlt[0] {
				variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
				aaPosOffset++
				if len(variant.AARef) == 0 {
					codonToAdd = nil
					for j = 0; j < 3; j++ {
						if (seqPos+offset)+j > currCDS.End-1 {
							seqPos = currCDS.Next.Start - 1
							currCDS = currCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][(seqPos-offset)+j])
					}
					dna.Complement(codonToAdd)
					variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}
		variant.AAPos = int(math.Round((float64(variant.CDNAPos) / 3) + 0.4)) + aaPosOffset // Add 0.4 so pos will always round up
	} else {
		var trimAA bool
		frame := determineFrame(variant)
		seqPos += frame
		lenOffset := (len(dna.StringToBases(variant.Ref)) - 1)

		if  int(variant.Pos - 1) + lenOffset > seqPos {
			seqPos += 3
			trimAA = true
			aaPosOffset--
		}
		for ; seqPos > (int(variant.Pos-1) + lenOffset); seqPos-- {
			if seqPos < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			} else if seqPos > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			refBases = append(refBases, seq[variant.Chr][seqPos])
			altBases = append(altBases, seq[variant.Chr][seqPos])
		}

		refBases = append(refBases, reverse(dna.StringToBases(variant.Ref))...)
		altBases = append(altBases, reverse(dna.StringToBases(variant.Alt))...)
		seqPos -= len(dna.StringToBases(variant.Ref))
		var offset int
		for offset = 0; len(altBases)%3 != 0; offset++ {
			if seqPos-offset < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			}
			altBases = append(altBases, seq[variant.Chr][seqPos-offset])
		}
		for offset = 0; len(refBases)%3 != 0; offset++ {
			if seqPos-offset < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			}
			refBases = append(refBases, seq[variant.Chr][seqPos-offset])
		}

		dna.Complement(refBases)
		dna.Complement(altBases)
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)

		if trimAA {
			variant.AARef = variant.AARef[:len(variant.AARef)-1]
		}

		if !isSynonymous(variant) && len(variant.AAAlt) > 1 {
			var codonToAdd []dna.Base
			var j int
			for len(variant.AAAlt) > 0 && variant.AARef[0] == variant.AAAlt[0] {
				variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
				aaPosOffset++
				if len(variant.AARef) == 0 {
					codonToAdd = nil
					for j = 0; j < 3; j++ {
						if (seqPos-offset)-j < currCDS.Start-1 {
							seqPos = currCDS.Prev.End - 1
							currCDS = currCDS.Prev
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][(seqPos-offset)-j])
					}
					dna.Complement(codonToAdd)
					variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}
		variant.AAPos = int(math.Round((float64(variant.CDNAPos) / 3) + 0.4)) + aaPosOffset // Add 0.4 so pos will always round up
	}
}

func reverse(s []dna.Base) []dna.Base {
	for i := 0; i < len(s) / 2; i++ {
		s[i], s[len(s)-1-i] = s[len(s)-1-i], s[i]
	}
	return s
}

func determineFrame(v *Variant) int {
	if v.PosStrand {
		return (((int(v.Pos)-v.NearestCDS.Start))%3 + ((3-v.NearestCDS.Frame)%3)) % 3
	} else {
		return (((v.NearestCDS.End-int(v.Pos)))%3 + ((3-v.NearestCDS.Frame)%3)) % 3
	}
}

func getCdsDist(v *Variant) int {
	switch {
	case int(v.Pos) >= v.NearestCDS.Start && int(v.Pos) <= v.NearestCDS.End: // Variant is in CDS
		return 0

	case int(v.Pos) < v.NearestCDS.Start: // Variant is before nearest CDS
		return v.NearestCDS.Start - int(v.Pos)

	default:
		return int(v.Pos) - v.NearestCDS.End  // Variant is after nearest CDS
	}
}

func isFrameshift(v *Variant) bool {
	refBases := dna.StringToBases(v.Ref)
	altBases := dna.StringToBases(v.Alt)

	start := int(v.Pos)
	refEnd := start + len(refBases) - 1
	altEnd := start + len(altBases) - 1

	var refBasesInCds int
	var altBasesInCds int

	var startOffset int
	if start < v.NearestCDS.Start {
		startOffset = v.NearestCDS.Start - start
	}

	if refEnd <= v.NearestCDS.End {
		refBasesInCds = len(refBases) - startOffset
	} else if refEnd > v.NearestCDS.End {
		refBasesInCds = len(refBases) - (refEnd - v.NearestCDS.End) - startOffset
	}
	if altEnd <= v.NearestCDS.End {
		altBasesInCds = len(altBases) - startOffset
	} else if altEnd > v.NearestCDS.End {
		altBasesInCds = len(altBases) - (altEnd - v.NearestCDS.End) - startOffset
	}

	shift := altBasesInCds - refBasesInCds
	return shift%3 != 0
}

func isNonsense(v *Variant) bool {
	for _, val := range v.AAAlt {
		if val == dna.Stop {
			return true
		}
	}
	return false
}

func isSynonymous(v *Variant) bool {
	var answer bool = true
	if len(v.AAAlt) != len(v.AARef) {
		return false
	} else {
		for i := 0; i < len(v.AARef); i++ {
			if v.AARef[i] != v.AAAlt[i] {
				answer = false
			}
		}
	}
	return answer
}
