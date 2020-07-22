package gtf

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"math"
	"reflect"
	"strings"
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
		for offset = 0; len(refBases)%3 != 0; offset++ {
			if seqPos+offset > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			refBases = append(refBases, seq[variant.Chr][seqPos+offset])
		}
		for offset = 0; len(altBases)%3 != 0; offset++ {
			if seqPos+offset > currCDS.End-1 {
				seqPos = currCDS.Next.Start - 1
				currCDS = currCDS.Next
			}
			altBases = append(altBases, seq[variant.Chr][seqPos+offset])
		}
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)
		variant.AAPos = int(math.Round((float64(variant.CDNAPos) / 3) + 0.4)) // Add 0.4 so pos will always round up
	} else {
		seqPos += determineFrame(variant)
		for ; seqPos > int(variant.Pos-1); seqPos-- {
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
		for offset = 0; len(refBases)%3 != 0; offset++ {
			if seqPos-offset < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			}
			refBases = append(refBases, seq[variant.Chr][seqPos-offset])
		}
		for offset = 0; len(altBases)%3 != 0; offset++ {
			if seqPos-offset < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
			}
			altBases = append(altBases, seq[variant.Chr][seqPos-offset])
		}
		dna.Complement(refBases)
		dna.Complement(altBases)
		variant.AARef = dna.TranslateSeq(refBases)
		variant.AAAlt = dna.TranslateSeq(altBases)
		variant.AAPos = int(math.Round((float64(variant.CDNAPos) / 3) + 0.4)) // Add 0.4 so pos will always round up
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
		return ((int(v.Pos)-v.NearestCDS.Start))%3 + ((3-v.NearestCDS.Frame)%3)
	} else {
		return ((v.NearestCDS.End-int(v.Pos)))%3 + ((3-v.NearestCDS.Frame)%3)
	}
}

// Annotation format is: Genomic | VariantType | Gene | cDNA | Protein
// TODO: Not sensitive to UTR splice junctions
func VariantToAnnotation(variant *Variant, seq map[string][]dna.Base) string {
	return addGenomic(variant) + "|" + variant.VariantType + "|" + strings.Trim(variant.Gene, "\"") + "|" + addCDNA(variant) + "|" + addProtein(variant, seq)
}

func addGenomic(v *Variant) string {
	return fmt.Sprintf("g.%s:%d%s>%s", v.Chr, v.Vcf.Pos, v.Vcf.Ref, v.Vcf.Alt)
}

func addCDNA(v *Variant) string {
	var answer string = strings.Trim(v.RefId, "\"") + ":c."
	if v.VariantType == "Intronic" || v.VariantType == "Splice" || v.VariantType == "FarSplice" {
		dist, start := getNearestCdsPos(v)
		if start {
			answer += fmt.Sprintf("%d-%d", dist, getCdsDist(v))
		} else {
			answer += fmt.Sprintf("%d+%d", dist, getCdsDist(v))
		}
	} else {
		answer += fmt.Sprintf("%d", v.CDNAPos)
	}
	if v.PosStrand {
		answer += v.Ref + ">" + v.Alt
	} else {
		ref := dna.StringToBases(v.Ref)
		alt := dna.StringToBases(v.Alt)
		dna.ReverseComplement(ref)
		dna.ReverseComplement(alt)
		answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
	}

	return answer
}

func truncateOnTer(a []dna.AminoAcid) []dna.AminoAcid {
	for i := 0 ; i < len(a) ; i++ {
		if a[i] == dna.Stop {
			return a[:i+1]
		}
	}
	return a
}

func addProtein(v *Variant, seq map[string][]dna.Base) string {
	// e.g. Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away), FarSplice (3-10 away)
	if v.VariantType == "Silent" || v.VariantType == "Missense" || v.VariantType == "Nonsense" || v.VariantType == "Frameshift" {
	} else {return ""}
	var answer string = "p."
	v.AAAlt = truncateOnTer(v.AAAlt)
	answer += fmt.Sprintf("%s%d", dna.AminoAcidToString(v.AARef[0]), v.AAPos)
	if len(v.AARef) > 1 {
		answer += "_" + dna.AminoAcidToString(v.AARef[len(v.AARef)-1]) + fmt.Sprintf("%d", v.AAPos + len(v.AARef) - 1)
	}

	refLen := len(v.AARef)
	altLen := len(v.AAAlt)
	switch {
	case refLen == 1 && altLen == 1: // Neither -> nothing is added
	case refLen == 1 && altLen > 1: // Insertion -> add "ins"
		answer += "ins"
	case refLen >= 1 && altLen <= 1: // Deletion -> add "del"
		answer += "del"
	case refLen >= 1 && altLen > 1: // Delins -> add "delins"
		answer += "delins"
	}

	if len(v.AAAlt) == 1 {
		answer += dna.AminoAcidToString(v.AAAlt[0])
	} else if len(v.AAAlt) > 5 {
		answer += fmt.Sprintf("%d", len(v.AAAlt))
	} else {
		for _, val := range v.AAAlt {
			answer += dna.AminoAcidToString(val)
		}
	}

	if v.VariantType == "Frameshift" {
		terDist := distToNextTer(v, seq)
		if terDist == 0 {
			v.VariantType = "Nonsense"
			return addProtein(v, seq)
		}
		answer += fmt.Sprintf("fsTer%d", terDist)
	}
	return answer
}

func addVariantType(v *Variant) {
	cdsDist := getCdsDist(v)
	switch {
	case v.Gene == "":
		v.VariantType = "Intergenic"
	case cdsDist > 0 && cdsDist <= 2:
		v.VariantType = "Splice"
	case cdsDist > 0 && cdsDist <= 10:
		v.VariantType = "FarSplice"
	case v.AARef == nil:
		v.VariantType = "Intronic"
	case isFrameshift(v):
		v.VariantType = "Frameshift"
	case isNonsense(v):
		v.VariantType = "Nonsense"
	case !reflect.DeepEqual(v.AARef, v.AAAlt):
		v.VariantType = "Missense"
	case reflect.DeepEqual(v.AARef, v.AAAlt):
		v.VariantType = "Silent"
	default:
		v.VariantType = "Unrecognized"
	}
}

func getNearestCdsPos(v *Variant) (pos int, start bool) {
	var currCDS *CDS = v.NearestCDS
	if v.PosStrand {
		if int(v.Pos) < v.NearestCDS.Start {
			pos = 1
		} else {
			pos = v.NearestCDS.End - v.NearestCDS.Start + 1
		}

		for currCDS.Prev != nil { // Go to first CDS to begin count
			currCDS = currCDS.Prev
			pos += currCDS.End - currCDS.Start
		}
		return pos, int(v.Pos) < v.NearestCDS.Start
	} else {
		if int(v.Pos) > v.NearestCDS.End {
			pos = 1
		} else {
			pos = v.NearestCDS.End - v.NearestCDS.Start + 1
		}

		for currCDS.Next != nil { // Go to first CDS to begin count
			currCDS = currCDS.Next
			pos += currCDS.End - currCDS.Start + 1
		}
		return pos, int(v.Pos) > v.NearestCDS.End
	}
}

func distToNextTer(v *Variant, seq map[string][]dna.Base) int {
	var answer int
	currSeq := seq[v.Chr]
	frame := determineFrame(v)
	var currCodon []dna.Base
	if v.PosStrand {
		for i := frame; i > 0; i-- {
			currCodon = append(currCodon, currSeq[int(v.Pos) - 1 - i])
		}
		seqPos := int(v.Pos) + len(dna.StringToBases(v.Ref)) - 1
		altSeq := dna.StringToBases(v.Alt)
		for _, val := range altSeq {
			currCodon = append(currCodon, val)
			if len(currCodon)%3 == 0 {
				if dna.TranslateSeq(currCodon)[0] == dna.Stop {
					return answer
				}
				answer++
				currCodon = nil
			}
		}
		currCDS := v.NearestCDS
		for {
			if seqPos > currCDS.End - 1 {
				currCDS = currCDS.Next
				seqPos = currCDS.Start - 1
			}
			currCodon = append(currCodon, currSeq[seqPos])
			seqPos++
			if len(currCodon)%3 == 0 {
				if dna.TranslateSeq(currCodon)[0] == dna.Stop {
					return answer
				}
				answer++
				currCodon = nil
			}
		}
	} else {
		for i := frame; i > 0; i-- {
			currCodon = append(currCodon, currSeq[int(v.Pos) - 1 + i])
		}
		seqPos := int(v.Pos) - 1 - len(dna.StringToBases(v.Ref))
		altSeq := reverse(dna.StringToBases(v.Alt))
		for _, val := range altSeq {
			currCodon = append(currCodon, val)
			if len(currCodon)%3 == 0 {
				dna.Complement(currCodon)
				if dna.TranslateSeq(currCodon)[0] == dna.Stop {
					return answer
				}
				answer++
				currCodon = nil
			}
		}
		currCDS := v.NearestCDS
		for {
			if seqPos < currCDS.Start - 1 { //TODO: add if currCDS.Prev != nil &&
				currCDS = currCDS.Prev
				fmt.Println(currCDS)
				seqPos = currCDS.End - 1
			}
			currCodon = append(currCodon, currSeq[seqPos])
			seqPos--
			if len(currCodon)%3 == 0 {
				dna.Complement(currCodon)
				//fmt.Println(dna.TranslateSeq(currCodon)[0])
				if dna.TranslateSeq(currCodon)[0] == dna.Stop {
					return answer
				}
				answer++
				currCodon = nil
			}
		}
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
