package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"math"
	"reflect"
	"strings"
)

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
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	if v.VariantType == "Intronic" || v.VariantType == "Splice" || v.VariantType == "FarSplice" {
		cdsPos, start := getNearestCdsPos(v)
		cdsDist := getCdsDist(v)
		if len(ref) == 1 && len(alt) == 1 { // Substitution
			if start {
				answer += fmt.Sprintf("%d-%d", cdsPos, cdsDist)
			} else {
				answer += fmt.Sprintf("%d+%d", cdsPos, cdsDist)
			}
		} else if len(ref) > len(alt) { // Deletion
			if v.PosStrand {
				if start {
					// TODO: test this on pos strand to make sure it is right
					answer += fmt.Sprintf("%d-%d_%d-%ddel",cdsPos ,cdsDist - 1 ,cdsPos, cdsDist + len(ref) - 1)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddel",cdsPos ,cdsDist - 1 ,cdsPos, cdsDist + len(ref) - 1)
				}
			} else {
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%ddel",cdsPos ,cdsDist - len(ref) + 1 ,cdsPos, cdsDist - 1)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddel",cdsPos ,cdsDist - len(ref) + 1 ,cdsPos, cdsDist - 1)
				}
			}
		} else if len(alt) > len(ref) { // Insertion
			//TODO Make addition
		}
	} else {
		answer += fmt.Sprintf("%d", v.CDNAPos)
	}

	if v.PosStrand {
		if len(ref) == 1 && len(alt) == 1 { // Substitution
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletion
			//TODO Make addition
		} else if len(alt) > len(ref) { // Insertion
			//TODO Make addition
		}

	} else {
		dna.ReverseComplement(ref)
		dna.ReverseComplement(alt)
		if len(ref) == 1 && len(alt) == 1 {
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletion
			//TODO Make addition
		} else if len(alt) > len(ref) { // Insertion
			//TODO Make addition
		}
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

func trimSynonymous(alpha []dna.AminoAcid, beta []dna.AminoAcid) (a, b []dna.AminoAcid) {
	if len(alpha) > 1 && len(beta) > 1 {
		for i := 0 ; i < common.Min(len(alpha), len(beta)); i++ {
			if alpha[i] != beta[i] {
				return alpha[i:], beta[i:]
			}
		}
	}
	return alpha, beta
}

func addProtein(v *Variant, seq map[string][]dna.Base) string {
	// e.g. Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away), FarSplice (3-10 away)
	if v.VariantType == "Silent" || v.VariantType == "Missense" || v.VariantType == "Nonsense" || v.VariantType == "Frameshift" {
	} else {return ""}
	var answer string = "p."
	v.AAAlt = truncateOnTer(v.AAAlt)
	v.AARef, v.AAAlt = trimSynonymous(v.AARef, v.AAAlt)
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
		terDist = terDist - (v.AAPos - int(math.Round((float64(v.CDNAPos) / 3) + 0.4)))
		if terDist == 1 {
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
	var answer int = 1
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
			if currCDS.Prev != nil && seqPos < currCDS.Start - 1 {
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
