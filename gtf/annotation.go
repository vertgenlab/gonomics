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
	return "GoEP=" + addGenomic(variant) + "|" + variant.VariantType + "|" + strings.Trim(variant.Gene, "\"") + "|" + addCDNA(variant, seq) + "|" + addProtein(variant, seq)
}

func addGenomic(v *Variant) string {
	return fmt.Sprintf("g.%s:%d%s>%s", v.Chr, v.Vcf.Pos, v.Vcf.Ref, v.Vcf.Alt)
}

func addCDNA(v *Variant, seq map[string][]dna.Base) string {
	var answer string = strings.Trim(v.RefId, "\"") + ":c."
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	cdsPos, start := getNearestCdsPos(v)
	if v.VariantType == "Intronic" || v.VariantType == "Splice" || v.VariantType == "FarSplice" {
		cdsDist := getCdsDist(v)
		if len(ref) == 1 && len(alt) == 1 { // Substitution
			if start {
				answer += fmt.Sprintf("%d-%d", cdsPos, cdsDist)
			} else {
				answer += fmt.Sprintf("%d+%d", cdsPos, cdsDist)
			}
		} else if len(ref) == 2 && len(alt) == 1 {
			if v.PosStrand {
				var duplicateOffset int
				for i := 1; seq[v.Chr][int(v.Pos)+i] == ref[1]; i++ {
					duplicateOffset++
				}
				if start {
					answer += fmt.Sprintf("%d-%ddel", cdsPos, (cdsDist - 1) - duplicateOffset)
				} else {
					answer += fmt.Sprintf("%d+%ddel", cdsPos, cdsDist + 1 + duplicateOffset)
				}
			} else {
				if start {
					answer += fmt.Sprintf("%d-%ddel", cdsPos, cdsDist + 1)
				} else {
					answer += fmt.Sprintf("%d+%ddel", cdsPos, cdsDist - 1)
				}
			}
			return answer
		} else if len(ref) > len(alt) { // Deletion
			if v.PosStrand {
				var duplicateOffset int
				for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(ref)-1)+j] == ref[i]; j++{
					duplicateOffset++
					i++
					if i == len(ref) {
						i = 1
					}
				}
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%ddel",cdsPos ,cdsDist - 1 - duplicateOffset, cdsPos, cdsDist - (len(ref) - 1) - duplicateOffset)
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%ddel",cdsPos ,cdsDist + 1 + duplicateOffset ,cdsPos, cdsDist + (len(ref) - 1) + duplicateOffset)
				}
			} else {
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%ddel",cdsPos ,cdsDist + len(ref) - 1 ,cdsPos, cdsDist + 1)
				} else {
					if cdsDist - len(ref) + 1 <= 0 {
						answer += fmt.Sprintf("%d_%d+%ddel",cdsPos + (cdsDist - len(ref) + 1) ,cdsPos, cdsDist - 1)
					} else {
						answer += fmt.Sprintf("%d+%d_%d+%ddel",cdsPos ,cdsDist - len(ref) + 1 ,cdsPos, cdsDist - 1)
					}
				}
			}
			return answer
		} else if isDuplication(v, seq) {
			if len(alt) == 2 {
				if v.PosStrand {
					var duplicateOffset int
					for i := 1; seq[v.Chr][int(v.Pos)+i] == alt[1]; i++ {
						duplicateOffset++
					}
					if start {
						answer += fmt.Sprintf("%d-%ddup",cdsPos ,cdsDist - 1 - duplicateOffset)
					} else {
						answer += fmt.Sprintf("%d+%ddup",cdsPos ,cdsDist + 1 + duplicateOffset)
					}
				} else {
					if start {
						answer += fmt.Sprintf("%d-%ddup",cdsPos ,cdsDist + len(alt) - 1)
					} else {
						answer += fmt.Sprintf("%d+%ddup",cdsPos ,cdsDist - 1)
					}
				}
			} else {
				if v.PosStrand {
					var duplicateOffset int
					for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+j] == alt[i]; j++{
						duplicateOffset++
						i++
						if i == len(alt) {
							i = 1
						}
					}
					if start {
						answer += fmt.Sprintf("%d-%d_%d-%ddup",cdsPos ,cdsDist - duplicateOffset + (len(alt) - 1) - 1 ,cdsPos, cdsDist - duplicateOffset)
					} else {
						answer += fmt.Sprintf("%d+%d_%d+%ddup",cdsPos ,cdsDist + (duplicateOffset - (len(alt) - 1)) + 1 ,cdsPos, cdsDist + duplicateOffset)
					}
				} else {
					if start {
						answer += fmt.Sprintf("%d-%d_%d-%ddup",cdsPos ,cdsDist + len(alt) - 1 ,cdsPos, cdsDist + 1)
					} else {
						answer += fmt.Sprintf("%d+%d_%d+%ddup",cdsPos ,cdsDist - len(alt) + 1 ,cdsPos, cdsDist - 1)
					}
				}
			}
			return answer
		} else if len(alt) > len(ref) { // Insertion
			if v.PosStrand {
				var duplicateOffset int
				for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+j] == alt[i]; j++{
					duplicateOffset++
					i++
					if i == len(alt) {
						i = 1
					}
				}
				//fmt.Println(cdsPos, cdsDist, duplicateOffset)
				basesToAdd := alt[len(ref):]
				if duplicateOffset > 0 {
					basesToAdd = append(basesToAdd[duplicateOffset:], basesToAdd[:duplicateOffset]...)
				}
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%dins%s",cdsPos ,cdsDist + len(ref) - 1 - duplicateOffset,cdsPos, cdsDist - 1 - duplicateOffset, dna.BasesToString(basesToAdd))
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%dins%s",cdsPos ,cdsDist + duplicateOffset ,cdsPos, cdsDist + 1 + duplicateOffset, dna.BasesToString(basesToAdd))
				}
			} else {
				tmpAlt := alt[len(ref):]
				dna.ReverseComplement(tmpAlt)
				if start {
					answer += fmt.Sprintf("%d-%d_%d-%dins%s",cdsPos ,cdsDist + 1 ,cdsPos, cdsDist + len(ref) - 1, dna.BasesToString(tmpAlt))
				} else {
					answer += fmt.Sprintf("%d+%d_%d+%dins%s",cdsPos ,cdsDist - 1 ,cdsPos, cdsDist - len(ref) + 1, dna.BasesToString(tmpAlt))
				}
			}
			return answer
		}
	}

	if v.PosStrand {
		if len(ref) == 1 && len(alt) == 1 { // Substitution
			if v.CDNAPos != 0 {
				answer += fmt.Sprintf("%d", v.CDNAPos)
			}
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletion
			var duplicateOffset int
			var hasDuplication bool
			for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(ref)-1)+j] == ref[i]; j++{
				duplicateOffset++
				i++
				if i >= 3 {
					hasDuplication = true
				}
				if i == len(ref) {
					i = 1
					hasDuplication = true
				}
			}

			if !hasDuplication {
				//duplicateOffset = 0
			}

			if len(ref) == 2 {
				answer += fmt.Sprintf("%ddel", v.CDNAPos + len(alt) + duplicateOffset)
			} else {
				//fmt.Println(v.CDNAPos, duplicateOffset)
				if v.CDNAPos + (len(ref) - 1) + duplicateOffset > cdsPos {
					answer += fmt.Sprintf("%d_%d+%ddel", v.CDNAPos + 1 + duplicateOffset, cdsPos, (v.CDNAPos + (len(ref) - 1) + duplicateOffset) - cdsPos)
				} else {
					answer += fmt.Sprintf("%d_%ddel", v.CDNAPos + 1 + duplicateOffset, v.CDNAPos + (len(ref) - 1) + duplicateOffset)
				}
			}
		} else if len(alt) > len(ref) { // Insertion
			if isDuplication(v, seq) {
				var duplicateOffset int
				for i, j := 1, 1; seq[v.Chr][int(v.Pos-1)+(len(alt)-1)+j] == alt[i]; j++{
					duplicateOffset++
					i++
					if i == len(alt) {
						i = 1
					}
				}
				if len(alt) == 2 {
					answer += fmt.Sprintf("%ddup", v.CDNAPos + duplicateOffset + 1)
				} else {
					answer += fmt.Sprintf("%d_%ddup", v.CDNAPos + duplicateOffset + 1, v.CDNAPos + duplicateOffset + 1 + (len(alt) - 2))
				}
			} else {
				answer += fmt.Sprintf("%d_%dins%s", v.CDNAPos, v.CDNAPos + 1, dna.BasesToString(alt[1:]))
			}
		}

	} else {
		dna.ReverseComplement(ref)
		dna.ReverseComplement(alt)
		if len(ref) == 1 && len(alt) == 1 {
			if v.CDNAPos != 0 {
				answer += fmt.Sprintf("%d", v.CDNAPos)
			}
			answer += dna.BasesToString(ref) + ">" + dna.BasesToString(alt)
		} else if len(ref) > len(alt) { // Deletion
			if len(ref) == 2 {
				answer += fmt.Sprintf("%ddel", v.CDNAPos - len(alt))
			} else {
				answer += fmt.Sprintf("%d_%ddel", v.CDNAPos - (len(ref) - 1), v.CDNAPos - 1)
			}
		} else if len(alt) > len(ref) { // Insertion
			if isDuplication(v, seq) {
				if len(alt) == 2 {
					answer += fmt.Sprintf("%ddup", v.CDNAPos - (len(alt) - 1))
				} else {
					answer += fmt.Sprintf("%d_%ddup", v.CDNAPos - (len(alt) - 1), v.CDNAPos - 1)
				}
			} else {
				answer += fmt.Sprintf("%d_%dins%s", v.CDNAPos - 1, v.CDNAPos, dna.BasesToString(alt[:len(alt)-1]))
			}
		}
	}

	return answer
}

func isDuplication(v *Variant, seq map[string][]dna.Base) bool {
	ref := dna.StringToBases(v.Ref)
	alt := dna.StringToBases(v.Alt)
	if len(ref) != 1 || len(ref) > len(alt) {
		return false
	}
	var answer bool = true
	var seqPos int = int(v.Pos - 1)
	for i := 0 ; i < len(alt); i++ {
		if alt[i] != seq[v.Chr][seqPos + i] {
			answer = false
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

	if v.VariantType == "Missense" && len(v.AAAlt) > 1 && v.AARef[0] == v.AAAlt[len(v.AAAlt)-1] {
		if len(v.AAAlt) - 1 == 1 {
			answer += fmt.Sprintf("%s%ddup", dna.AminoAcidToString(v.AAAlt[0]), v.AAPos - 1)
		} else {
			answer += fmt.Sprintf("%s%d_%s%ddup", dna.AminoAcidToString(v.AAAlt[0]), v.AAPos - (len(v.AAAlt) - 1), dna.AminoAcidToString(v.AAAlt[len(v.AAAlt)-2]), v.AAPos - 1)
		}
		return answer
	}

	if v.VariantType == "Missense" && len(v.AARef) == 1 && len(v.AAAlt) == 0 && len(dna.StringToBases(v.Ref)) > 3 {
		answer += fmt.Sprintf("%s%ddel", dna.AminoAcidToString(v.AARef[0]), v.CDNAPos/3)
		return answer
	}

	if v.VariantType == "Missense" && len(v.AARef) == 0 {
		answer += fmt.Sprintf("%s%ddup", dna.AminoAcidToString(v.AAAlt[0]), v.CDNAPos/3)
		return answer
	}

	v.AAAlt = truncateOnTer(v.AAAlt)
	v.AARef, v.AAAlt = trimSynonymous(v.AARef, v.AAAlt)

	//fmt.Println(v.Info)
	//fmt.Println(dna.PolypeptideToString(v.AARef))
	//fmt.Println(dna.PolypeptideToString(v.AAAlt))
	//fmt.Println(v.Vcf)
	//fmt.Println(v)

	answer += fmt.Sprintf("%s%d", dna.AminoAcidToString(v.AARef[0]), v.AAPos)

	if v.VariantType == "Nonsense" {
		answer += "Ter"
		return answer
	}

	if len(v.AARef) > 1 && v.VariantType != "Frameshift" {
		answer += "_" + dna.AminoAcidToString(v.AARef[len(v.AARef)-1]) + fmt.Sprintf("%d", v.AAPos + len(v.AARef) - 1)
	}

	refLen := len(v.AARef)
	altLen := len(v.AAAlt)
	switch {
	case refLen == 1 && altLen == 1: // Neither -> nothing is added
	case refLen == 1 && altLen > 1 && v.VariantType != "Frameshift": // Insertion -> add "ins"
		answer += "delins"
	case refLen >= 1 && altLen == 0 && v.VariantType != "Frameshift": // Deletion -> add "del"
		answer += "del"
	case refLen >= 1 && altLen >= 1 && v.VariantType != "Frameshift": // Delins -> add "delins"
		answer += "delins"
	}

	if len(v.AAAlt) == 1 || (len(v.AAAlt) > 1 && v.VariantType == "Frameshift") {
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

		if terDist == 1 || (len(v.AAAlt) > 0 && v.AAAlt[0] == dna.Stop) {
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
			pos += currCDS.End - currCDS.Start + 1
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
