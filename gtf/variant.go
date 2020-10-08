package gtf

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"math"
	"reflect"
)

type vcfEffectPrediction struct {
	vcf.Vcf
	RefId          string // e.g. NC_000023.10, LRG_199, NG_012232.1, NM_004006.2, LRG-199t1, NR_002196.1, NP_003997.1, etc.
	Gene           string
	PosStrand      bool
	NearestCds     *CDS
	CdnaPos        int // 1-base
	AaPos          int // 1-base
	AaRef          []dna.AminoAcid
	AaAlt          []dna.AminoAcid
	VariantType    string // e.g. Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away), FarSplice (3-10 away)
	NextTranscript *vcfEffectPrediction
}

// GenesToIntervalTree builds a fractionally cascaded 2d interval tree for efficiently identifying genes that overlap a variant
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

// VcfToVariant determines the effects of a variant on the cDNA and amino acid sequence by querying genes in the tree made by GenesToIntervalTree
// Note that if multiple genes are found to overlap a variant this function will return the variant based on the first queried gene and throw an error
// All bases in fasta record must be uppercase
func VcfToVariant(v *vcf.Vcf, tree map[string]*interval.IntervalNode, seq map[string][]dna.Base, allTranscripts bool) (*vcfEffectPrediction, error) {
	var answer *vcfEffectPrediction
	var err error

	overlappingGenes := interval.Query(tree, v, "any")

	if len(overlappingGenes) > 1 {
		err = errors.New(fmt.Sprintf("Variant overlaps with multiple genes. Mutation will be based on first gene."))
	}

	var annotatingGene *Gene
	if len(overlappingGenes) > 0 {
		annotatingGene = overlappingGenes[0].(*Gene)
		answer = vcfToVariant(v, annotatingGene, seq, allTranscripts)
	} else {
		answer = &vcfEffectPrediction{Vcf: *v}
	}
	return answer, err
}

// vcfToVariant is a helper function that annotates the Variant struct with information from the Vcf and Gtf input
func vcfToVariant(v *vcf.Vcf, gene *Gene, seq map[string][]dna.Base, allTranscripts bool) *vcfEffectPrediction {
	answer := new(vcfEffectPrediction)
	answer.Vcf = *v
	answer.RefId = gene.Transcripts[0].TranscriptID
	answer.Gene = gene.GeneID
	answer.PosStrand = gene.Transcripts[0].Strand
	vcfCdsIntersect(v, gene, answer, 0)
	if answer.NearestCds != nil { // handle genes with no CDS
		if int(v.Pos) >= answer.NearestCds.Start && int(v.Pos) <= answer.NearestCds.End {
			findAAChange(answer, seq)
		}
		addVariantType(answer)
	}

	if allTranscripts {
		var prev *vcfEffectPrediction = answer
		for i := 1; i < len(gene.Transcripts); i++ {
			additionalTranscript := new(vcfEffectPrediction)
			additionalTranscript.Vcf = *v
			additionalTranscript.RefId = gene.Transcripts[i].TranscriptID
			additionalTranscript.Gene = gene.GeneID
			additionalTranscript.PosStrand = gene.Transcripts[i].Strand
			vcfCdsIntersect(v, gene, additionalTranscript, i)
			if answer.NearestCds != nil {
				if int(v.Pos) >= additionalTranscript.NearestCds.Start && int(v.Pos) <= additionalTranscript.NearestCds.End {
					findAAChange(additionalTranscript, seq)
				}
				addVariantType(additionalTranscript)
			}
			prev.NextTranscript = additionalTranscript
			prev = additionalTranscript
		}
	}

	return answer
}

// vcfCdsIntersect annotates the Variant struct with the cDNA position of the vcf as well as the CDS nearest to the vcf
func vcfCdsIntersect(v *vcf.Vcf, gene *Gene, answer *vcfEffectPrediction, transcriptPosInSlice int) {
	var cdsPos int
	var exon *Exon
	//TODO: this code may be able to be compressed
	if answer.PosStrand {
		for i := 0; i < len(gene.Transcripts[transcriptPosInSlice].Exons); i++ {
			exon = gene.Transcripts[transcriptPosInSlice].Exons[i]
			if exon.Cds != nil && int(v.Pos) > exon.Cds.End { // variant is further in gene
				cdsPos += exon.Cds.End - exon.Cds.Start + 1
				answer.NearestCds = exon.Cds // Store most recent exon and move on // Catches variants past the last exon
			} else if exon.Cds != nil && int(v.Pos) <= exon.Cds.End { // variant is before end of this exon
				if int(v.Pos) < exon.Cds.Start { // Variant is NOT in CDS
					if exon.Cds.Prev == nil || exon.Cds.Start-int(v.Pos) < int(v.Pos)-gene.Transcripts[transcriptPosInSlice].Exons[i-1].Cds.Start {
						answer.NearestCds = exon.Cds
					} else {
						answer.NearestCds = gene.Transcripts[transcriptPosInSlice].Exons[i-1].Cds
					}
					break
				}
				cdsPos += int(v.Pos) - exon.Cds.Start + 1
				answer.CdnaPos = cdsPos
				answer.NearestCds = exon.Cds
			}
		}
	} else {
		for i := 0; i < len(gene.Transcripts[transcriptPosInSlice].Exons); i++ {
			exon = gene.Transcripts[transcriptPosInSlice].Exons[len(gene.Transcripts[transcriptPosInSlice].Exons)-1-i]
			if exon.Cds != nil && int(v.Pos) < exon.Cds.Start { // variant is further in gene
				cdsPos += exon.Cds.End - exon.Cds.Start + 1
				answer.NearestCds = exon.Cds // Store most recent exon and move on // Catches variants past the last exon
			} else if exon.Cds != nil && int(v.Pos) >= exon.Cds.Start { // variant is before end of this exon
				if int(v.Pos) > exon.Cds.End { // Variant is NOT in CDS
					if exon.Cds.Next == nil || int(v.Pos)-exon.Cds.End < gene.Transcripts[transcriptPosInSlice].Exons[len(gene.Transcripts[transcriptPosInSlice].Exons)-1-i+1].Cds.Start-int(v.Pos) {
						answer.NearestCds = exon.Cds
					} else {
						answer.NearestCds = gene.Transcripts[transcriptPosInSlice].Exons[len(gene.Transcripts[transcriptPosInSlice].Exons)-1-i+1].Cds
					}
					break
				}
				// Variant IS in CDS
				cdsPos += exon.Cds.End - int(v.Pos) + 1
				answer.CdnaPos = cdsPos
				answer.NearestCds = exon.Cds
			}
		}
	}

}

// findAAChange annotates the Variant struct with the amino acids changed by a given variant
func findAAChange(variant *vcfEffectPrediction, seq map[string][]dna.Base) {
	ref := dna.StringToBases(variant.Ref)
	alt := dna.StringToBases(variant.Alt)
	var refBases = make([]dna.Base, 0)
	var altBases = make([]dna.Base, 0)
	var seqPos int = int(variant.Pos) - 1
	var currCDS *CDS = variant.NearestCds
	var aaPosOffset int = 0
	if variant.PosStrand {
		seqPos -= determineFrame(variant)

		var hasDuplication bool
		var duplicateOffset int
		var duplicateBasePos int
		if (len(ref)-len(alt))%3 == 0 && len(ref) > 1 {
			var j int
			for duplicateBasePos, j = 1, 1; seq[variant.Chr][int(variant.Pos-1)+(len(ref)-1)+j] == ref[duplicateBasePos]; j++ {
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

			variant.CdnaPos += duplicateOffset
			variant.Pos += int64(duplicateOffset)
			seqPos = int(variant.Pos) - 1
			seqPos -= determineFrame(variant)
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
			refBases = append(refBases, ref[duplicateBasePos-1:]...)
			if duplicateBasePos-1 > 0 {
				refBases = append(refBases, ref[1:duplicateBasePos-1]...)
				seqPos -= len(ref[1 : duplicateBasePos-1])
			}
			altBases = append(altBases, alt[1:]...)
		} else {
			refBases = append(refBases, ref...)
			altBases = append(altBases, alt...)
		}

		seqPos += len(ref)

		var offset int
		altCDS := currCDS
		altSeqPos := seqPos
		for ; len(altBases)%3 != 0; altSeqPos++ {
			if altSeqPos > altCDS.End-1 && altCDS.Next != nil {
				altSeqPos = altCDS.Next.Start - 1
				altCDS = altCDS.Next
			}
			altBases = append(altBases, seq[variant.Chr][altSeqPos])
		}
		refCDS := currCDS
		refSeqPos := seqPos
		for ; len(refBases)%3 != 0; refSeqPos++ {
			if refSeqPos > refCDS.End-1 && refCDS.Next != nil {
				refSeqPos = refCDS.Next.Start - 1
				refCDS = refCDS.Next
			}
			refBases = append(refBases, seq[variant.Chr][refSeqPos])
		}
		variant.AaRef = dna.TranslateSeq(refBases)
		variant.AaAlt = dna.TranslateSeq(altBases)

		if (len(ref)-len(alt))%3 != 0 {
			var codonToAdd []dna.Base
			var j int
			for variant.AaRef[0] == variant.AaAlt[0] {
				codonToAdd = nil
				variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
				aaPosOffset++
				if len(variant.AaRef) == 0 {
					for j = 0; j < 3; j++ {
						if refSeqPos > refCDS.End-1 && refCDS.Next != nil {
							refSeqPos = refCDS.Next.Start - 1
							refCDS = refCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][refSeqPos])
						refSeqPos++
					}
					variant.AaRef = append(variant.AaRef, dna.TranslateSeq(codonToAdd)...)
				}
				codonToAdd = nil
				if len(variant.AaAlt) == 0 {
					for j = 0; j < 3; j++ {
						if altSeqPos > altCDS.End-1 && altCDS.Next != nil {
							altSeqPos = altCDS.Next.Start - 1
							altCDS = altCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][altSeqPos])
						altSeqPos++
					}
					variant.AaAlt = append(variant.AaAlt, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}

		if !isSynonymous(variant) && len(variant.AaRef) > 1 {
			var codonToAdd []dna.Base
			var j int
			for len(variant.AaAlt) > 0 && variant.AaRef[0] == variant.AaAlt[0] {
				variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
				aaPosOffset++
				if len(variant.AaRef) == 0 {
					codonToAdd = nil
					for j = 0; j < 3; j++ {
						if (seqPos+offset)+j > currCDS.End-1 {
							seqPos = currCDS.Next.Start - 1
							currCDS = currCDS.Next
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][(seqPos+offset)+j])
					}
					variant.AaRef = append(variant.AaRef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		}

		variant.AaPos = int(math.Round((float64(variant.CdnaPos)/3)+0.4)) + aaPosOffset // Add 0.4 so pos will always round up
	} else {
		var trimAA bool
		seqPos += determineFrame(variant)
		lenOffset := (len(dna.StringToBases(variant.Ref)) - 1)

		for int(variant.Pos-1)+lenOffset > seqPos {
			seqPos += 3
			trimAA = true
			aaPosOffset--
		}

		if seqPos > currCDS.End-1 {
			seqPos = (currCDS.Next.Start - 1) + ((seqPos - int(variant.Pos)) - (currCDS.End - int(variant.Pos)))
			currCDS = currCDS.Next
		}

		for ; seqPos > (int(variant.Pos-1) + lenOffset); seqPos-- {
			if seqPos < currCDS.Start-1 {
				seqPos = currCDS.Prev.End - 1
				currCDS = currCDS.Prev
				if seqPos == (int(variant.Pos-1) + lenOffset) {
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
			altBases = append(altBases, seq[variant.Chr][altSeqPos])
		}
		refCDS := currCDS
		refSeqPos := seqPos
		for ; len(refBases)%3 != 0; refSeqPos-- {
			if refSeqPos < refCDS.Start-1 {
				refSeqPos = refCDS.Prev.End - 1
				refCDS = refCDS.Prev
			}
			refBases = append(refBases, seq[variant.Chr][refSeqPos])
		}

		dna.Complement(refBases)
		dna.Complement(altBases)
		variant.AaRef = dna.TranslateSeq(refBases)
		variant.AaAlt = dna.TranslateSeq(altBases)

		if trimAA && (len(ref)-len(alt))%3 == 0 && variant.AaRef[len(variant.AaRef)-1] == variant.AaAlt[len(variant.AaAlt)-1] {
			variant.AaAlt = variant.AaAlt[:len(variant.AaAlt)-1]
			variant.AaRef = variant.AaRef[:len(variant.AaRef)-1]
		}

		var codonToAdd []dna.Base
		var j int

		if !isSynonymous(variant) && len(variant.AaAlt) > 1 && len(variant.AaRef) > 0 {
			for len(variant.AaRef) > 0 && len(variant.AaAlt) > 0 && variant.AaRef[0] == variant.AaAlt[0] {
				if len(variant.AaAlt) > 1 && variant.AaRef[0] == variant.AaAlt[len(variant.AaAlt)-1] && (len(ref)-len(alt))%3 == 0 {
					variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
					aaPosOffset++
					break
				}
				variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
				aaPosOffset++
				if len(variant.AaRef) == 0 {
					codonToAdd = nil
					for j = 0; j < 3; j++ {
						if (refSeqPos)-j < currCDS.Start-1 {
							seqPos = currCDS.Prev.End - 1
							currCDS = currCDS.Prev
						}
						codonToAdd = append(codonToAdd, seq[variant.Chr][(refSeqPos)-j])
					}
					dna.Complement(codonToAdd)
					variant.AaRef = append(variant.AaRef, dna.TranslateSeq(codonToAdd)...)
				}
			}
		} else if !isSynonymous(variant) && len(variant.AaAlt) == 1 && len(variant.AaRef) == 1 && variant.AaAlt[0] == variant.AaRef[0] && len(ref) > len(alt) {
			if trimAA {
				refSeqPos += 3
			}
			variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
			aaPosOffset++
			codonToAdd = nil
			for j = 0; j < 3; j++ {
				if (refSeqPos)-j < currCDS.Start-1 {
					seqPos = currCDS.Prev.End - 1
					currCDS = currCDS.Prev
				}
				codonToAdd = append(codonToAdd, seq[variant.Chr][(refSeqPos)-j])
			}
			dna.Complement(codonToAdd)
			variant.AaRef = append(variant.AaRef, dna.TranslateSeq(codonToAdd)...)
		}

		if (len(ref)-len(alt))%3 != 0 && len(variant.AaRef) > 0 && len(variant.AaAlt) > 0 && variant.AaRef[0] == variant.AaAlt[0] {
			if trimAA {
				trimAA = false
				refSeqPos += 3
			}
			codonToAdd = nil
			variant.AaRef, variant.AaAlt = variant.AaRef[1:], variant.AaAlt[1:]
			aaPosOffset++
			for len(codonToAdd) == 0 || len(codonToAdd)%3 != 0 {
				codonToAdd = append(codonToAdd, seq[variant.Chr][refSeqPos])
				refSeqPos--
				if refSeqPos < refCDS.Start-1 {
					refSeqPos = refCDS.Prev.End - 1
					refCDS = refCDS.Prev
				}
			}
			dna.Complement(codonToAdd)
			variant.AaRef = append(variant.AaRef, dna.TranslateSeq(codonToAdd)...)
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
			variant.AaAlt = append(variant.AaAlt, dna.TranslateSeq(codonToAdd)...)
		}

		variant.AaPos = int(math.Round((float64(variant.CdnaPos)/3)+0.4)) + aaPosOffset // Add 0.4 so pos will always round up
	}
}

// addVariantType annotates the Variant struct with the VariantType
// Valid types include: Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice, FarSplice
// Splice is defined as 1-2 bases away from intron-exon border, Farsplice is 3-10 bases away from intron-exon border
func addVariantType(v *vcfEffectPrediction) {
	cdsDist := getCdsDist(v)
	switch {
	case v.Gene == "":
		v.VariantType = "Intergenic"
	case cdsDist > 0 && cdsDist <= 2:
		v.VariantType = "Splice"
	case cdsDist > 0 && cdsDist <= 10:
		v.VariantType = "FarSplice"
	case v.AaRef == nil:
		v.VariantType = "Intronic"
	case isFrameshift(v):
		v.VariantType = "Frameshift"
	case isNonsense(v):
		v.VariantType = "Nonsense"
	case !reflect.DeepEqual(v.AaRef, v.AaAlt):
		v.VariantType = "Missense"
	case reflect.DeepEqual(v.AaRef, v.AaAlt):
		v.VariantType = "Silent"
	default:
		v.VariantType = "Unrecognized"
	}
}

// reverse reverses the order of a slice of dna.Base
// e.g. [0 1 2] -> [2 1 0]
func reverse(s []dna.Base) []dna.Base {
	for i := 0; i < len(s)/2; i++ {
		s[i], s[len(s)-1-i] = s[len(s)-1-i], s[i]
	}
	return s
}

// determineFrame will determine the position of the variant in a codon
// This is used to determine how many bases before the variant must be retrieved to get the full codon
func determineFrame(v *vcfEffectPrediction) int {
	if v.PosStrand {
		return ((int(v.Pos)-v.NearestCds.Start)%3 + ((3 - v.NearestCds.Frame) % 3)) % 3
	} else {
		return ((v.NearestCds.End-int(v.Pos))%3 + ((3 - v.NearestCds.Frame) % 3)) % 3
	}
}

// getCdsDist determines the distance of the variant from the nearest CDS
// Returns 0 if the variant is inside the CDS
func getCdsDist(v *vcfEffectPrediction) int {
	switch {
	case int(v.Pos) >= v.NearestCds.Start && int(v.Pos) <= v.NearestCds.End: // Variant is in CDS
		return 0

	case int(v.Pos) < v.NearestCds.Start: // Variant is before nearest CDS
		return v.NearestCds.Start - int(v.Pos)

	default:
		return int(v.Pos) - v.NearestCds.End // Variant is after nearest CDS
	}
}

// isFrameshift returns true if the variant shifts the reading frame
func isFrameshift(v *vcfEffectPrediction) bool {
	refBases := dna.StringToBases(v.Ref)
	altBases := dna.StringToBases(v.Alt)

	start := int(v.Pos)
	refEnd := start + len(refBases) - 1

	var refBasesInCds int
	var altBasesInCds int

	var startOffset int
	if start < v.NearestCds.Start {
		startOffset = v.NearestCds.Start - start
	}

	if refEnd <= v.NearestCds.End {
		refBasesInCds = len(refBases) - startOffset
	} else if refEnd > v.NearestCds.End {
		refBasesInCds = len(refBases) - (refEnd - v.NearestCds.End) - startOffset
	}
	altBasesInCds = len(altBases) - startOffset
	shift := altBasesInCds - refBasesInCds
	return shift%3 != 0
}

// isNonsense returns true if the variant creates a premature stop codon
func isNonsense(v *vcfEffectPrediction) bool {
	for _, val := range v.AaAlt {
		if val == dna.Stop {
			return true
		}
	}
	return false
}

// isSynonymous returns true if the variant does not change the amino acid sequence
func isSynonymous(v *vcfEffectPrediction) bool {
	var answer bool = true
	if len(v.AaAlt) != len(v.AaRef) || len(dna.StringToBases(v.Ref)) != len(dna.StringToBases(v.Alt)) {
		return false
	} else {
		for i := 0; i < len(v.AaRef); i++ {
			if v.AaRef[i] != v.AaAlt[i] {
				answer = false
			}
		}
	}
	return answer
}
