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

// NOTE: All bases in fasta record must be uppercase
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

// TODO: WILL REQUIRE CHANGES FOR POS STRAND
func findAAChange(variant *Variant, seq map[string][]dna.Base) {
	ref := dna.StringToBases(variant.Ref)
	alt := dna.StringToBases(variant.Alt)
	var refBases = make([]dna.Base, 0)
	var altBases = make([]dna.Base, 0)
	var seqPos int = int(variant.Pos) - 1
	var currCDS *CDS = variant.NearestCDS
	var aaPosOffset int = 0
	if variant.PosStrand {
		//fmt.Println("NEW VARIANT")
		seqPos -= determineFrame(variant)

		var hasDuplication bool
		var duplicateOffset int
		var duplicateBasePos int
		if (len(ref) - len(alt))%3 == 0 && len(ref) > 1 {
			var j int
			//fmt.Println("MAKE ADJUSTMENT ON:", variant.Info)
			//fmt.Println("ORIG POS:", seqPos + 1)
			for duplicateBasePos, j = 1, 1; seq[variant.Chr][int(variant.Pos-1)+(len(ref)-1)+j] == ref[duplicateBasePos]; j++{
				duplicateOffset++
				duplicateBasePos++

				if duplicateBasePos >= 3 {
					hasDuplication = true
				}

				if duplicateBasePos == len(ref) {
					duplicateBasePos = 1
				}
			}

			if !hasDuplication {
				duplicateOffset = 0
				duplicateBasePos = 0
			}
			//fmt.Println(variant.CDNAPos)

			variant.CDNAPos += duplicateOffset
			variant.Pos += int64(duplicateOffset)
			seqPos = int(variant.Pos) - 1
			seqPos -= determineFrame(variant)
			//variant.Pos -= int64(duplicateOffset)
			//fmt.Println(seqPos + 1, duplicateOffset)
		}

		for ; seqPos < int(variant.Pos-1); seqPos++ {
			if seqPos < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			} else if seqPos > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
				if seqPos <= int(variant.Pos-1) {
					break
				}
			}
			refBases = append(refBases, seq[variant.Chr][seqPos])
			altBases = append(altBases, seq[variant.Chr][seqPos])
		}
		if duplicateOffset > 0 {
			//fmt.Println(duplicateBasePos)
			refBases = append(refBases, ref[duplicateBasePos-1:]...)
			if duplicateBasePos-1 > 0 {
				refBases = append(refBases, ref[1:duplicateBasePos-1]...)
				seqPos -= len(ref[1:duplicateBasePos-1])
			}
			altBases = append(altBases, alt[1:]...)
		} else {
			refBases = append(refBases, ref...)
			altBases = append(altBases, alt...)
		}

		//if duplicateOffset != 0 {
		//	fmt.Println("NEW VAR ON", variant.Info)
		//	fmt.Println(duplicateOffset)
		//	fmt.Println(dna.BasesToString(refBases))
		//	fmt.Println(dna.BasesToString(altBases))
		//}

		//fmt.Println(len(refBases), len(altBases))
		seqPos += len(ref)

		var offset int
		altCDS := currCDS
		altSeqPos := seqPos
		for ; len(altBases)%3 != 0; altSeqPos++ {
			if altSeqPos > altCDS.End-1 {
				altSeqPos = altCDS.Next.Start - 1
				altCDS = altCDS.Next
			}
			altBases = append(altBases, seq[variant.Chr][altSeqPos])
		}
		refCDS := currCDS
		refSeqPos := seqPos
		for ; len(refBases)%3 != 0; refSeqPos++ {
			if refSeqPos > refCDS.End-1 {
				refSeqPos = refCDS.Next.Start - 1
				refCDS = refCDS.Next
			}
			//if duplicateOffset != 0 {
			//	fmt.Println(seqPos + 1)
			//}
			refBases = append(refBases, seq[variant.Chr][refSeqPos])
		}
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)

		//if duplicateOffset != 0 {
		//	fmt.Println(dna.BasesToString(refBases))
		//	fmt.Println(dna.BasesToString(altBases))
		//
		//	fmt.Println("REF:", dna.PolypeptideToString(variant.AARef))
		//	fmt.Println("ALT:", dna.PolypeptideToString(variant.AAAlt))
		//}

		//fmt.Println(len(dna.StringToBases(variant.Ref)), len(dna.StringToBases(variant.Alt)))

		if (len(ref) - len(alt))%3 != 0 {
			var codonToAdd []dna.Base
			var j int
			for variant.AARef[0] == variant.AAAlt[0] {
				codonToAdd = nil
				variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
				//fmt.Println("BASE ADDED")
				aaPosOffset++
				if len(variant.AARef) == 0 {
					for j = 0 ; j < 3; j++ {
						if refSeqPos > refCDS.End-1 {
							refSeqPos = refCDS.Next.Start - 1
							refCDS = refCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][refSeqPos])
						refSeqPos++
					}
					variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
				}
				codonToAdd = nil
				if len(variant.AAAlt) == 0 {
					for j = 0 ; j < 3; j++ {
						if altSeqPos > altCDS.End-1 {
							altSeqPos = altCDS.Next.Start - 1
							altCDS = altCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][altSeqPos])
						altSeqPos++
					}
					variant.AAAlt = append(variant.AAAlt, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}

		//fmt.Println("END REF:", dna.PolypeptideToString(variant.AARef))
		//fmt.Println("END ALT:", dna.PolypeptideToString(variant.AAAlt))

		if !isSynonymous(variant) && len(variant.AARef) > 1 {
			//fmt.Println("ATTEMPTED ADDITION ON:", variant.Info)
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
						codonToAdd = append(codonToAdd, seq[variant.Chr][(seqPos+offset)+j])
					}
					variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}

		//fmt.Println("END REF:", dna.PolypeptideToString(variant.AARef))
		//fmt.Println("END ALT:", dna.PolypeptideToString(variant.AAAlt))

		variant.AAPos = int(math.Round((float64(variant.CDNAPos) / 3) + 0.4)) + aaPosOffset // Add 0.4 so pos will always round up
	} else {
		var trimAA bool
		seqPos += determineFrame(variant)
		lenOffset := (len(dna.StringToBases(variant.Ref)) - 1)

		//if !isRepetetive(variant, seq) {
			for int(variant.Pos-1)+lenOffset > seqPos {
				seqPos += 3
				trimAA = true
				aaPosOffset--
			}
		//}

		if seqPos > currCDS.End - 1 {
			seqPos = (currCDS.Next.Start - 1) + ((seqPos - int(variant.Pos)) - (currCDS.End - int(variant.Pos)))
			currCDS = currCDS.Next
		}

		for ; seqPos > (int(variant.Pos-1) + lenOffset); seqPos-- {
			if seqPos < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
				if seqPos == (int(variant.Pos - 1)+ lenOffset) {
					break
				}
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

		altCDS := currCDS
		altSeqPos := seqPos
		for ; len(altBases)%3 != 0; altSeqPos-- {
			if altSeqPos < altCDS.Start-1 {
				altSeqPos = altCDS.Prev.End - 1
				altCDS = altCDS.Prev
			}
			//if variant.Info == "c.3747del|p.Phe1250LeufsTer13" || variant.Info == "c.5604del|p.Gln1869SerfsTer42" ||
			//	variant.Info == "c.8307_8308del|p.Ala2770HisfsTer4" || variant.Info == "c.13104del|p.Val4369CysfsTer25" {
			//	fmt.Println("ALT:", altSeqPos + 1)
			//}
			altBases = append(altBases, seq[variant.Chr][altSeqPos])
		}
		refCDS := currCDS
		refSeqPos := seqPos
		for ; len(refBases)%3 != 0; refSeqPos-- {
			if refSeqPos < refCDS.Start-1 {
				refSeqPos = refCDS.Prev.End - 1
				refCDS = refCDS.Prev
			}
			//if variant.Info == "c.3747del|p.Phe1250LeufsTer13" || variant.Info == "c.5604del|p.Gln1869SerfsTer42" ||
			//	variant.Info == "c.8307_8308del|p.Ala2770HisfsTer4" || variant.Info == "c.13104del|p.Val4369CysfsTer25" {
			//	fmt.Println("REF:", refSeqPos + 1)
			//}
			refBases = append(refBases, seq[variant.Chr][refSeqPos])
		}

		//if variant.Info == "c.3747del|p.Phe1250LeufsTer13" || variant.Info == "c.5604del|p.Gln1869SerfsTer42" ||
		//	variant.Info == "c.8307_8308del|p.Ala2770HisfsTer4" || variant.Info == "c.13104del|p.Val4369CysfsTer25" {
		//	fmt.Println(dna.BasesToString(refBases))
		//	fmt.Println(dna.BasesToString(altBases))
		//}
		//fmt.Println(dna.BasesToString(refBases))
		//fmt.Println(dna.BasesToString(altBases))

		dna.Complement(refBases)
		dna.Complement(altBases)
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)

		//if variant.Info == "c.2625_2627dup|p.Thr876dup" ||
		//	variant.Info == "c.33513_33515dup|p.Glu11172dup" {
		//	fmt.Println(dna.PolypeptideToString(variant.AARef))
		//	fmt.Println(dna.PolypeptideToString(variant.AAAlt))
		//	fmt.Println(trimAA)
		//}

		//fmt.Println(len(variant.AARef), len(variant.AAAlt))

		if trimAA {
			if (len(ref) - len(alt))%3 == 0 {
				if variant.AARef[len(variant.AARef)-1] == variant.AAAlt[len(variant.AAAlt)-1] {
					variant.AAAlt = variant.AAAlt[:len(variant.AAAlt)-1]
					variant.AARef = variant.AARef[:len(variant.AARef)-1]
				} else {
					variant.AARef = variant.AARef[:len(variant.AARef)-1]
				}
			}
			//if variant.AARef[len(variant.AARef)-1] == variant.AAAlt[len(variant.AAAlt)-1] &&
			//	(len(ref) - len(alt))%3 == 0 {
			//	variant.AAAlt = variant.AAAlt[:len(variant.AAAlt)-1]
			//
			//}
			//variant.AARef = variant.AARef[:len(variant.AARef)-1]
		}

		//if variant.Info == "c.1758_1760del|p.Glu586_Thr587delinsAsp" {
		//	fmt.Println(dna.PolypeptideToString(variant.AARef))
		//	fmt.Println(dna.PolypeptideToString(variant.AAAlt))
		//	fmt.Println(variant.Vcf)
		//}

		//fmt.Println(len(variant.AARef), len(variant.AAAlt))
		var codonToAdd []dna.Base
		var j int

		//if variant.Info == "c.107219dup|p.Pro35741AlafsTer20" {
		//	fmt.Println("REF:", dna.PolypeptideToString(variant.AARef))
		//	fmt.Println("ALT:", dna.PolypeptideToString(variant.AAAlt))
		//}

		if !isSynonymous(variant) && len(variant.AAAlt) > 1 && len(variant.AARef) > 0{
			for len(variant.AARef) > 0 && len(variant.AAAlt) > 0 && variant.AARef[0] == variant.AAAlt[0] {
				if len(variant.AAAlt) == 2 && variant.AARef[0] == variant.AAAlt[1] {
					variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
					aaPosOffset++
					break
				}
				variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
				aaPosOffset++
				if len(variant.AARef) == 0 {
					codonToAdd = nil
					for j = 0; j < 3; j++ {
						if (refSeqPos)-j < currCDS.Start-1 {
							seqPos = currCDS.Prev.End - 1
							currCDS = currCDS.Prev
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][(refSeqPos)-j])
					}
					dna.Complement(codonToAdd)
					variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}

		//if variant.Info == "c.3747del|p.Phe1250LeufsTer13" || variant.Info == "c.5604del|p.Gln1869SerfsTer42" ||
		//	variant.Info == "c.8307_8308del|p.Ala2770HisfsTer4" || variant.Info == "c.13104del|p.Val4369CysfsTer25" {
		//	fmt.Println(dna.PolypeptideToString(variant.AARef))
		//	fmt.Println(dna.PolypeptideToString(variant.AAAlt))
		//}

		//if variant.Info == "c.33513_33515del|p.Glu11172del" {
		//	fmt.Println(dna.PolypeptideToString(variant.AARef))
		//	fmt.Println(dna.PolypeptideToString(variant.AAAlt))
		//}

		if (len(ref) - len(alt))%3 != 0 && len(variant.AARef) > 0 && len(variant.AAAlt) > 0 && variant.AARef[0] == variant.AAAlt[0] {
			if trimAA {
				trimAA = false
				refSeqPos += 3
			}
			codonToAdd = nil
			variant.AARef, variant.AAAlt = variant.AARef[1:], variant.AAAlt[1:]
			aaPosOffset++
			for len(codonToAdd) == 0 || len(codonToAdd)%3 != 0 {
				//if variant.Info == "c.3747del|p.Phe1250LeufsTer13" || variant.Info == "c.5604del|p.Gln1869SerfsTer42" ||
				//	variant.Info == "c.8307_8308del|p.Ala2770HisfsTer4" || variant.Info == "c.13104del|p.Val4369CysfsTer25" {
				//	fmt.Println(refSeqPos + 1)
				//}
				codonToAdd = append(codonToAdd, seq[variant.Chr][refSeqPos])
				refSeqPos--
				if refSeqPos < refCDS.Start-1 {
					refSeqPos = refCDS.Prev.End - 1
					refCDS = refCDS.Prev
				}
			}
			dna.Complement(codonToAdd)
			variant.AARef = append(variant.AARef, dna.TranslateSeq(codonToAdd)...)
			codonToAdd = nil
			for len(codonToAdd) == 0 || len(codonToAdd)%3 != 0 {
				codonToAdd = append(codonToAdd, seq[variant.Chr][altSeqPos])
				altSeqPos--
				if altSeqPos < altCDS.Start-1 {
					altSeqPos = altCDS.Prev.End - 1
					altCDS = altCDS.Prev
				}
			}
			dna.Complement(codonToAdd)
			variant.AAAlt = append(variant.AAAlt, dna.TranslateSeq(codonToAdd)...)
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

//func isRepetetive(v *Variant, seq map[string][]dna.Base) bool {
//	var answer = true
//	ref := dna.StringToBases(v.Ref)
//	for i := 1; i < len(ref); i++ {
//		if seq[v.Chr][int(v.Pos - 1) + i] != seq[v.Chr][int(v.Pos - 1) + len(ref) + i] {
//			answer = false
//			break
//		}
//	}
//	return answer
//}