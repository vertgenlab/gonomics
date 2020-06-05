package dna

import (
	"log"
)

type Codon struct {
	Seq []Base
}

//moving Sophie's code for AminoAcid identities and the genetic code to the DNA package.
type AminoAcid uint8

const (
	Ala  AminoAcid = 0
	Arg  AminoAcid = 1
	Asn  AminoAcid = 2
	Asp  AminoAcid = 3
	Cys  AminoAcid = 4
	Gln  AminoAcid = 5
	Glu  AminoAcid = 6
	Gly  AminoAcid = 7
	His  AminoAcid = 8
	Ile  AminoAcid = 9
	Leu  AminoAcid = 10
	Lys  AminoAcid = 11
	Met  AminoAcid = 12
	Phe  AminoAcid = 13
	Pro  AminoAcid = 14
	Ser  AminoAcid = 15
	Thr  AminoAcid = 16
	Trp  AminoAcid = 17
	Tyr  AminoAcid = 18
	Val  AminoAcid = 19
	Stop AminoAcid = 20
)

var GeneticCode = map[string]AminoAcid{
	"TGA": AminoAcid(20), "TAA": AminoAcid(20), "TAG": AminoAcid(20),
	"GTA": AminoAcid(19), "GTC": AminoAcid(19), "GTG": AminoAcid(19), "GTT": AminoAcid(19),
	"TAT": AminoAcid(18), "TAC": AminoAcid(18),
	"TGG": AminoAcid(17),
	"ACA": AminoAcid(16), "ACG": AminoAcid(16), "ACT": AminoAcid(16), "ACC": AminoAcid(16),
	"TCA": AminoAcid(15), "TCC": AminoAcid(15), "TCG": AminoAcid(15), "TCT": AminoAcid(15), "AGT": AminoAcid(15), "AGC": AminoAcid(15),
	"CCC": AminoAcid(14), "CCT": AminoAcid(14), "CCA": AminoAcid(14), "CCG": AminoAcid(14),
	"TTT": AminoAcid(13), "TTC": AminoAcid(13),
	"ATG": AminoAcid(12),
	"AAA": AminoAcid(11), "AAG": AminoAcid(11),
	"TTA": AminoAcid(10), "TTG": AminoAcid(10), "CTC": AminoAcid(10), "CTG": AminoAcid(10), "CTA": AminoAcid(10), "CTT": AminoAcid(10),
	"ATT": AminoAcid(9), "ATC": AminoAcid(9), "ATA": AminoAcid(9),
	"CAT": AminoAcid(8), "CAC": AminoAcid(8),
	"GGG": AminoAcid(7), "GGA": AminoAcid(7), "GGT": AminoAcid(7), "GGC": AminoAcid(7),
	"GAA": AminoAcid(6), "GAG": AminoAcid(6),
	"CAA": AminoAcid(5), "CAG": AminoAcid(5),
	"TGT": AminoAcid(4), "TGC": AminoAcid(4),
	"GAT": AminoAcid(3), "GAC": AminoAcid(3),
	"AAT": AminoAcid(2), "AAC": AminoAcid(2),
	"AGA": AminoAcid(1), "AGG": AminoAcid(1), "CGC": AminoAcid(1), "CGG": AminoAcid(1), "CGA": AminoAcid(1), "CGT": AminoAcid(1),
	"GCA": AminoAcid(0), "GCG": AminoAcid(0), "GCT": AminoAcid(0), "GCC": AminoAcid(0)}

func AminoAcidToShortString(a AminoAcid) string {
	switch a {
	case Ala:
		return "A"
	case Arg:
		return "R"
	case Asn:
		return "N"
	case Asp:
		return "D"
	case Cys:
		return "C"
	case Gln:
		return "Q"
	case Glu:
		return "E"
	case Gly:
		return "G"
	case His:
		return "H"
	case Ile:
		return "I"
	case Leu:
		return "L"
	case Lys:
		return "K"
	case Met:
		return "M"
	case Phe:
		return "F"
	case Pro:
		return "P"
	case Ser:
		return "S"
	case Thr:
		return "T"
	case Trp:
		return "W"
	case Tyr:
		return "Y"
	case Val:
		return "V"
	case Stop:
		return "*"
	default:
		log.Fatalf("Error: unexpected amino acid value %v\n", a)
		return "X"
	}
}

func AminoAcidToString(a AminoAcid) string {
	switch a {
	case Ala:
		return "Ala"
	case Arg:
		return "Arg"
	case Asn:
		return "Asn"
	case Asp:
		return "Asp"
	case Cys:
		return "Cys"
	case Gln:
		return "Gln"
	case Glu:
		return "Glu"
	case Gly:
		return "Gly"
	case His:
		return "His"
	case Ile:
		return "Ile"
	case Leu:
		return "Leu"
	case Lys:
		return "Lys"
	case Met:
		return "Met"
	case Phe:
		return "Phe"
	case Pro:
		return "Pro"
	case Ser:
		return "Ser"
	case Thr:
		return "Thr"
	case Trp:
		return "Trp"
	case Tyr:
		return "Tyr"
	case Val:
		return "Val"
	case Stop:
		return "Ter"
	default:
		log.Fatalf("Error: unexpected amino acid value %v\n", a)
		return "X"
	}
}

//TODO: Add in frame as an argument.
//TODO: Initialize to len 3 and use the frame argument to update the index instead of append.
func BasesToCodons(b []Base) []*Codon {
	var frame int
	var answer []*Codon
	var current []Base
	for i := 0; i < len(b); i++ {
		if frame == 2 {
			current = append(current, b[i])
			answer = append(answer, &Codon{Seq: current})
			current = make([]Base, 0) //clear current
			frame = 0
		} else {
			current = append(current, b[i])
			frame++
		}
	}
	return answer
}

func CodonsToSeq(c []*Codon) []Base {
	var answer []Base
	for i := 0; i < len(c); i++ {
		for j := 0; j < len(c[i].Seq); j++ {
			answer = append(answer, c[i].Seq[j])
		}
	}
	return answer
}

func TranslateCodon(c *Codon) AminoAcid {
	val, ok := GeneticCode[BasesToString(c.Seq)]
	if ok {
		return val
	} else {
		log.Fatalf("Codon not in map.\n")
		return val
	}
}

func NonSynonymous(c1 *Codon, c2 *Codon) bool {
	return TranslateCodon(c1) != TranslateCodon(c2)
}

func Synonymous(c1 *Codon, c2 *Codon) bool {
	return TranslateCodon(c1) == TranslateCodon(c2) && !IsEqual(c1, c2)
}

func IsEqual(c1 *Codon, c2 *Codon) bool {
	return BasesToString(c1.Seq) == BasesToString(c2.Seq)
}

func TranslateSeq(b []Base) []AminoAcid {
	var answer []AminoAcid
	codons := BasesToCodons(b)
	for i := 0; i < len(codons); i++ {
		answer = append(answer, TranslateCodon(codons[i]))
	}
	return answer
}

func PolypeptideToShortString(a []AminoAcid) string {
	var answer string
	for i := 0; i < len(a); i++ {
		answer = answer + AminoAcidToShortString(a[i])
	}
	return answer
}

func PolypeptideToString(a []AminoAcid) string {
	var answer string
	for i := 0; i < len(a); i++ {
		answer = answer + AminoAcidToString(a[i])
	}
	return answer
}

func TranslateToShortString(b []Base) string {
	AllToUpper(b)
	a := TranslateSeq(b)
	return PolypeptideToShortString(a)
}

func TranslateToString(b []Base) string {
	AllToUpper(b)
	a := TranslateSeq(b)
	return PolypeptideToString(a)
}
