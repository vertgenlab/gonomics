package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"math"
)

// VariantToAnnotation generates an annotation which can be appended to the INFO field of a VCF
// Annotation format is: GoEP= Genomic | VariantType | Gene | cDNA | Protein
// Genomic cDNA and Protein annotations are in HGVS variant nomenclature format https://varnomen.hgvs.org/
// TODO: Not sensitive to UTR splice junctions
func VariantToAnnotation(variant *vcfEffectPrediction, seq map[string][]dna.Base) string {
	return "GoEP=" + genomicToString(variant) + "|" + variant.VariantType + "|" + variant.Gene + "|" + cDnaToString(variant, seq) + "|" + proteinToString(variant, seq)
}

// genomicToString adds the genotype portion of the annotation
func genomicToString(v *vcfEffectPrediction) string {
	return fmt.Sprintf("g.%s:%d%s>%s", v.Chr, v.Vcf.Pos, v.Vcf.Ref, v.Vcf.Alt)
}

// cDnaToString adds the cDNA portion of the annotation and makes several adjustments for standard variant nomenclature
func cDnaToString(v *vcfEffectPrediction, seq map[string][]dna.Base) string {
	if v.VariantType == "Intronic" || v.VariantType == "Splice" || v.VariantType == "FarSplice" {
		return nonCodingToString(v, seq)
	} else {
		return codingToString(v, seq)
	}
}

// nonCodingToString inputs intronic variants and generates a cDNA annotation
func nonCodingToString(v *vcfEffectPrediction, seq map[string][]dna.Base) string {
	var answer string = v.RefId + ":c."
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	cdsPos, start := getNearestCdsPos(v)
	cdsDist := getCdsDist(v)
	if len(ref) == 1 && len(alt) == 1 { // Substitution: c.10-1A>G
		if start {
			answer += fmt.Sprintf("%d-%d", cdsPos, cdsDist)
		} else {
			answer += fmt.Sprintf("%d+%d", cdsPos, cdsDist)
		}
		if v.PosStrand {
			answer += v.Ref + ">" + v.Alt
		} else {
			dna.ReverseComplement(ref)
			dna.ReverseComplement(alt)
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		}
		return answer
	} else if len(ref) == 2 && len(alt) == 1 { // Single-Base Deletion: c.10-1del
		if v.PosStrand {
			var duplicateOffset int
			for i := 1; seq[v.Chr][int(v.Pos)+i] == ref[1]; i++ {
				duplicateOffset++
			}
			if start {
				answer += fmt.Sprintf("%d-%ddel", cdsPos, (cdsDist-1)-duplicateOffset)
			} else {
				answer += fmt.Sprintf("%d+%ddel", cdsPos, cdsDist+1+duplicateOffset)
			}
		} else {
			if start {
				answer += fmt.Sprintf("%d-%ddel", cdsPos, cdsDist+1)
			} else {
				answer += fmt.Sprintf("%d+%ddel", cdsPos, cdsDist-1)
			}
		}
		return answer
	} else if len(ref) > len(alt) { // Multi-Base Deletion: c.10-1_10-2del
		if v.PosStrand {
			var duplicateOffset int
			for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(ref)-1)+j] == ref[i]; j++ {
				duplicateOffset++
				i++
				if i == len(ref) {
					i = 1
				}
			}
			if start {
				answer += fmt.Sprintf("%d-%d_%d-%ddel", cdsPos, cdsDist-1-duplicateOffset, cdsPos, cdsDist-(len(ref)-1)-duplicateOffset)
			} else {
				answer += fmt.Sprintf("%d+%d_%d+%ddel", cdsPos, cdsDist+1+duplicateOffset, cdsPos, cdsDist+(len(ref)-1)+duplicateOffset)
			}
		} else {
			if start {
				answer += fmt.Sprintf("%d-%d_%d-%ddel", cdsPos, cdsDist+len(ref)-1, cdsPos, cdsDist+1)
			} else {
				if cdsDist-len(ref)+1 <= 0 {
					answer += fmt.Sprintf("%d_%d+%ddel", cdsPos+(cdsDist-len(ref)+1), cdsPos, cdsDist-1)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddel", cdsPos, cdsDist-len(ref)+1, cdsPos, cdsDist-1)
				}
			}
		}
		return answer
	} else if isDuplication(v, seq) { // Duplications
		if len(alt) == 2 { // Single-Base Duplication: c.10-1dup
			if v.PosStrand {
				var duplicateOffset int
				for i := 1; seq[v.Chr][int(v.Pos)+i] == alt[1]; i++ {
					duplicateOffset++
				}
				if start {
					answer += fmt.Sprintf("%d-%ddup", cdsPos, cdsDist-1-duplicateOffset)
				} else {
					answer += fmt.Sprintf("%d+%ddup", cdsPos, cdsDist+1+duplicateOffset)
				}
			} else {
				if start {
					answer += fmt.Sprintf("%d-%ddup", cdsPos, cdsDist+len(alt)-1)
				} else {
					answer += fmt.Sprintf("%d+%ddup", cdsPos, cdsDist-1)
				}
			}
		} else {
			if v.PosStrand {
				var duplicateOffset int
				for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+j] == alt[i]; j++ {
					duplicateOffset++
					i++
					if i == len(alt) {
						i = 1
					}
				}
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%ddup", cdsPos, cdsDist-duplicateOffset+(len(alt)-1)-1, cdsPos, cdsDist-duplicateOffset)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddup", cdsPos, cdsDist+(duplicateOffset-(len(alt)-1))+1, cdsPos, cdsDist+duplicateOffset)
				}
			} else {
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%ddup", cdsPos, cdsDist+len(alt)-1, cdsPos, cdsDist+1)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddup", cdsPos, cdsDist-len(alt)+1, cdsPos, cdsDist-1)
				}
			}
		}
		return answer
	} else if len(alt) > len(ref) { // Multi-Base Duplication: c.10-1_10-2dup
		if v.PosStrand {
			var duplicateOffset int
			for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+j] == alt[i]; j++ {
				duplicateOffset++
				i++
				if i == len(alt) {
					i = 1
				}
			}
			basesToAdd := alt[len(ref):]
			if duplicateOffset > 0 {
				basesToAdd = append(basesToAdd[duplicateOffset:], basesToAdd[:duplicateOffset]...)
			}
			if start {
				answer += fmt.Sprintf("%d-%d_%d-%dins%s", cdsPos, cdsDist+len(ref)-1-duplicateOffset, cdsPos, cdsDist-1-duplicateOffset, dna.BasesToString(basesToAdd))
			} else {
				answer += fmt.Sprintf("%d+%d_%d+%dins%s", cdsPos, cdsDist+duplicateOffset, cdsPos, cdsDist+1+duplicateOffset, dna.BasesToString(basesToAdd))
			}
		} else {
			tmpAlt := alt[len(ref):]
			dna.ReverseComplement(tmpAlt)
			if start {
				answer += fmt.Sprintf("%d-%d_%d-%dins%s", cdsPos, cdsDist+1, cdsPos, cdsDist+len(ref)-1, dna.BasesToString(tmpAlt))
			} else {
				answer += fmt.Sprintf("%d+%d_%d+%dins%s", cdsPos, cdsDist-1, cdsPos, cdsDist-len(ref)+1, dna.BasesToString(tmpAlt))
			}
		}
		return answer
	}
	return answer
}

// codingToString inputs variants inside the CDS and generates a cDNA annotation
func codingToString(v *vcfEffectPrediction, seq map[string][]dna.Base) string {
	var answer string = v.RefId + ":c."
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	cdsPos, _ := getNearestCdsPos(v)
	if v.PosStrand { // For Pos Strand
		if len(ref) == 1 && len(alt) == 1 { // Substitution c.10A>G
			if v.CdnaPos != 0 {
				answer += fmt.Sprintf("%d", v.CdnaPos)
			}
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletions
			var duplicateOffset int
			for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(ref)-1)+j] == ref[i]; j++ {
				duplicateOffset++
				i++
				if i == len(ref) {
					i = 1
				}
			}

			if len(ref) == 2 { // Single-Base Deletion: c.10del
				answer += fmt.Sprintf("%ddel", v.CdnaPos+len(alt)+duplicateOffset)
			} else { // Multi-Base Deletion: c.10_12del
				if v.CdnaPos+(len(ref)-1)+duplicateOffset > cdsPos {
					answer += fmt.Sprintf("%d_%d+%ddel", v.CdnaPos+1+duplicateOffset, cdsPos, (v.CdnaPos+(len(ref)-1)+duplicateOffset)-cdsPos)
				} else {
					answer += fmt.Sprintf("%d_%ddel", v.CdnaPos+1+duplicateOffset, v.CdnaPos+(len(ref)-1)+duplicateOffset)
				}
			}
		} else if len(alt) > len(ref) {
			if isDuplication(v, seq) { // Duplications
				var duplicateOffset int
				for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(alt)-1)+j] == alt[i]; j++ {
					duplicateOffset++
					i++
					if i == len(alt) {
						i = 1
					}
				}
				if len(alt) == 2 { // Single-Base Duplication: c.10dup
					answer += fmt.Sprintf("%ddup", v.CdnaPos+duplicateOffset+1)
				} else { // Multi-Base Duplication: c.10_12dup
					answer += fmt.Sprintf("%d_%ddup", v.CdnaPos+duplicateOffset+1, v.CdnaPos+duplicateOffset+1+(len(alt)-2))
				}
			} else { // Insertion: c.10_11insAT
				answer += fmt.Sprintf("%d_%dins%s", v.CdnaPos, v.CdnaPos+1, dna.BasesToString(alt[1:]))
			}
		}
	} else { // For Neg Strand
		dna.ReverseComplement(ref)
		dna.ReverseComplement(alt)
		if len(ref) == 1 && len(alt) == 1 { // Substitution: c.10A>G
			if v.CdnaPos != 0 {
				answer += fmt.Sprintf("%d", v.CdnaPos)
			}
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletions
			if len(ref) == 2 { // Single-Base Deletion: c.10del
				answer += fmt.Sprintf("%ddel", v.CdnaPos-len(alt))
			} else { // Multi-Base Deletion: c.10_12del
				answer += fmt.Sprintf("%d_%ddel", v.CdnaPos-(len(ref)-1), v.CdnaPos-1)
			}
		} else if len(alt) > len(ref) {
			if isDuplication(v, seq) { // Duplications
				if len(alt) == 2 { // Single-Base Duplication: c.10dup
					answer += fmt.Sprintf("%ddup", v.CdnaPos-(len(alt)-1))
				} else { // Multi-Base Duplication: c.10_12dup
					answer += fmt.Sprintf("%d_%ddup", v.CdnaPos-(len(alt)-1), v.CdnaPos-1)
				}
			} else { // Insertion: c.10_11insAT
				answer += fmt.Sprintf("%d_%dins%s", v.CdnaPos-1, v.CdnaPos, dna.BasesToString(alt[:len(alt)-1]))
			}
		}
	}
	return answer
}

// isDuplication returns true if the variant is likely a duplication of surrounding sequence
func isDuplication(v *vcfEffectPrediction, seq map[string][]dna.Base) bool {
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	if len(ref) > len(alt) {
		return false
	}
	var seqPos int = int(v.Pos - 1)
	for i := 0; i < len(alt); i++ {
		if alt[i] != seq[v.Chr][seqPos+i] {
			return false
		}
	}
	return true
}

// truncateOnTer inputs a slice of amino acids and truncates the slice at the first stop codon
func truncateOnTer(a []dna.AminoAcid) []dna.AminoAcid {
	for i := 0; i < len(a); i++ {
		if a[i] == dna.Stop {
			return a[:i+1]
		}
	}
	return a
}

// trimSynonymous removes identical amino acids that are present in the ref and alt amino acid slices
func trimSynonymous(alpha []dna.AminoAcid, beta []dna.AminoAcid) (a, b []dna.AminoAcid) {
	if len(alpha) > 1 && len(beta) > 1 {
		for i := 0; i < common.Min(len(alpha), len(beta)); i++ {
			if alpha[i] != beta[i] {
				return alpha[i:], beta[i:]
			}
		}
	}
	return alpha, beta
}

// proteinToString adds the protein portion of the annotation and makes several adjustments for standard variant nomenclature
func proteinToString(v *vcfEffectPrediction, seq map[string][]dna.Base) string {
	if v.VariantType == "Silent" || v.VariantType == "Missense" || v.VariantType == "Nonsense" || v.VariantType == "Frameshift" {
	} else {
		return ""
	}
	var answer string = "p."
	if v.VariantType == "Missense" && len(v.AaRef) == 0 {
		if len(v.AaAlt) == 1 {
			answer += fmt.Sprintf("%s%ddup", dna.AminoAcidToString(v.AaAlt[0]), v.CdnaPos/3)
		} else {
			fmt.Println(v.AaPos)
			answer += fmt.Sprintf("%s%d_%s%ddup", dna.AminoAcidToString(v.AaAlt[0]), v.AaPos, dna.AminoAcidToString(v.AaAlt[len(v.AaAlt)-1]), (v.CdnaPos/3)+len(v.AaAlt))
		}
		return answer
	}

	if v.VariantType == "Missense" && len(v.AaAlt) > 1 && v.AaRef[0] == v.AaAlt[len(v.AaAlt)-1] {
		if len(v.AaAlt)-1 == 1 {
			answer += fmt.Sprintf("%s%ddup", dna.AminoAcidToString(v.AaAlt[0]), v.AaPos-1)
		} else {
			answer += fmt.Sprintf("%s%d_%s%ddup", dna.AminoAcidToString(v.AaAlt[0]), v.AaPos-(len(v.AaAlt)-1), dna.AminoAcidToString(v.AaAlt[len(v.AaAlt)-2]), v.AaPos-1)
		}
		return answer
	}

	if v.VariantType == "Missense" && len(v.AaRef) == 1 && len(v.AaAlt) == 0 && len(dna.StringToBases(v.Ref)) > 3 {
		answer += fmt.Sprintf("%s%ddel", dna.AminoAcidToString(v.AaRef[0]), v.CdnaPos/3)
		return answer
	}

	if v.Info == "c.38753_38784del|p.Leu12918CysfsTer2" {
		fmt.Println(dna.PolypeptideToString(v.AaRef))
		fmt.Println(dna.PolypeptideToString(v.AaAlt))
	}

	v.AaAlt = truncateOnTer(v.AaAlt)
	v.AaRef, v.AaAlt = trimSynonymous(v.AaRef, v.AaAlt)

	answer += fmt.Sprintf("%s%d", dna.AminoAcidToString(v.AaRef[0]), v.AaPos)

	if v.VariantType == "Nonsense" {
		answer += "Ter"
		return answer
	}

	if len(v.AaRef) > 1 && v.VariantType != "Frameshift" {
		answer += "_" + dna.AminoAcidToString(v.AaRef[len(v.AaRef)-1]) + fmt.Sprintf("%d", v.AaPos+len(v.AaRef)-1)
	}

	refLen := len(v.AaRef)
	altLen := len(v.AaAlt)
	switch {
	case refLen == 1 && altLen == 1: // Neither -> nothing is added
	case refLen == 1 && altLen > 1 && v.VariantType != "Frameshift": // Insertion -> add "ins"
		answer += "delins"
	case refLen >= 1 && altLen == 0 && v.VariantType != "Frameshift": // Deletion -> add "del"
		answer += "del"
	case refLen >= 1 && altLen >= 1 && v.VariantType != "Frameshift": // Delins -> add "delins"
		answer += "delins"
	}

	if len(v.AaAlt) == 1 || (len(v.AaAlt) > 1 && v.VariantType == "Frameshift") {
		answer += dna.AminoAcidToString(v.AaAlt[0])
	} else if len(v.AaAlt) > 5 {
		answer += fmt.Sprintf("%d", len(v.AaAlt))
	} else {
		for _, val := range v.AaAlt {
			answer += dna.AminoAcidToString(val)
		}
	}

	if v.VariantType == "Frameshift" {
		terDist := distToNextTer(v, seq)
		terDist = terDist - (v.AaPos - int(math.Round((float64(v.CdnaPos)/3)+0.4)))

		if terDist == 1 || (len(v.AaAlt) > 0 && v.AaAlt[0] == dna.Stop) {
			v.VariantType = "Nonsense"
			return proteinToString(v, seq)
		}
		answer += fmt.Sprintf("fsTer%d", terDist)
	}
	return answer
}

// getNearestCdsPos determines the cDNA position of the nearest CDS end and whether the reported end is the 5' or 3' end
// Used for reporting intronic variants
func getNearestCdsPos(v *vcfEffectPrediction) (pos int, start bool) {
	var currCDS *CDS = v.NearestCds
	if v.PosStrand {
		if int(v.Pos) < v.NearestCds.Start {
			pos = 1
		} else {
			pos = v.NearestCds.End - v.NearestCds.Start + 1
		}

		for currCDS.Prev != nil { // Go to first CDS to begin count
			currCDS = currCDS.Prev
			pos += currCDS.End - currCDS.Start + 1
		}
		return pos, int(v.Pos) < v.NearestCds.Start
	} else {
		if int(v.Pos) > v.NearestCds.End {
			pos = 1
		} else {
			pos = v.NearestCds.End - v.NearestCds.Start + 1
		}

		for currCDS.Next != nil { // Go to first CDS to begin count
			currCDS = currCDS.Next
			pos += currCDS.End - currCDS.Start + 1
		}
		return pos, int(v.Pos) > v.NearestCds.End
	}
}

// distToNextTer determines the distance to the nearest stop codon
// Used for annotating frameshift variants
func distToNextTer(v *vcfEffectPrediction, seq map[string][]dna.Base) int {
	var answer int = 1
	currSeq := seq[v.Chr]
	var currCodon []dna.Base
	var originalFrame int = determineFrame(v)
	if v.PosStrand {
		frame := originalFrame
		for i := frame; i > 0; i-- {
			currCodon = append(currCodon, currSeq[int(v.Pos)-1-i])
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
		currCDS := v.NearestCds
		for {
			if seqPos > currCDS.End-1 {
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
		refLen := len(dna.StringToBases(v.Ref))
		altSeq := reverse(dna.StringToBases(v.Alt))

		if (refLen-1)-originalFrame > 0 {
			answer -= 1 + (((refLen - 2) - originalFrame) / 3)
		}

		frame := ((v.NearestCds.End-(int(v.Pos)+len(dna.StringToBases(v.Ref))-1))%3 + ((3 - v.NearestCds.Frame) % 3)) % 3
		for i := frame; i > 0; i-- {
			currCodon = append(currCodon, currSeq[int(v.Pos)+i])
		}

		var seqPos int
		if len(altSeq) < len(dna.StringToBases(v.Ref)) {
			seqPos = int(v.Pos) - 2
		} else {
			seqPos = int(v.Pos) - 1 - len(dna.StringToBases(v.Ref))
		}

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
		currCDS := v.NearestCds
		for {
			if currCDS.Prev != nil && seqPos < currCDS.Start-1 {
				currCDS = currCDS.Prev
				seqPos = currCDS.End - 1
			}

			currCodon = append(currCodon, currSeq[seqPos])
			seqPos--
			if len(currCodon)%3 == 0 {
				dna.Complement(currCodon)
				if dna.TranslateSeq(currCodon)[0] == dna.Stop {
					return answer
				}
				answer++
				currCodon = nil
			}
		}
	}
}
