package dna

import (
	"fmt"
	"log"
)

type Codon struct {
	Seq []Base
}

//moving Sophie's code for aa identities and the genetic code to the DNA package.
type aa int

const (
	ala  aa = 0
	arg  aa = 1
	asn  aa = 2
	asp  aa = 3
	cys  aa = 4
	gln  aa = 5
	glu  aa = 6
	gly  aa = 7
	his  aa = 8
	ile  aa = 9
	leu  aa = 10
	lys  aa = 11
	met  aa = 12
	phe  aa = 13
	pro  aa = 14
	ser  aa = 15
	thr  aa = 16
	trp  aa = 17
	tyr  aa = 18
	val  aa = 19
	stop aa = 20
)

var GeneticCode = map[string]aa{
	"TGA": aa(20), "TAA": aa(20), "TAG": aa(20),
	"GTA": aa(19), "GTC": aa(19), "GTG": aa(19), "GTT": aa(19),
	"TAT": aa(18), "TAC": aa(18),
	"TGG": aa(17),
	"ACA": aa(16), "ACG": aa(16), "ACT": aa(16), "ACC": aa(16),
	"TCA": aa(15), "TCC": aa(15), "TCG": aa(15), "TCT": aa(15), "AGT": aa(15), "AGC": aa(15),
	"CCC": aa(14), "CCT": aa(14), "CCA": aa(14), "CCG": aa(14),
	"TTT": aa(13), "TTC": aa(13),
	"ATG": aa(12),
	"AAA": aa(11), "AAG": aa(11),
	"TTA": aa(10), "TTG": aa(10), "CTC": aa(10), "CTG": aa(10), "CTA": aa(10), "CTT": aa(10),
	"ATT": aa(9), "ATC": aa(9), "ATA": aa(9),
	"CAT": aa(8), "CAC": aa(8),
	"GGG": aa(7), "GGA": aa(7), "GGT": aa(7), "GGC": aa(7),
	"GAA": aa(6), "GAG": aa(6),
	"CAA": aa(5), "CAG": aa(5),
	"TGT": aa(4), "TGC": aa(4),
	"GAT": aa(3), "GAC": aa(3),
	"AAT": aa(2), "AAC": aa(2),
	"AGA": aa(1), "AGG": aa(1), "CGC": aa(1), "CGG": aa(1), "CGA": aa(1), "CGT": aa(1),
	"GCA": aa(0), "GCG": aa(0), "GCT": aa(0), "GCC": aa(0)}

func AAToString(a aa) string {
	switch a {
	case ala:
		return "A"
	case arg:
		return "R"
	case asn:
		return "N"
	case asp:
		return "D"
	case cys:
		return "C"
	case gln:
		return "Q"
	case glu:
		return "E"
	case gly:
		return "G"
	case his:
		return "H"
	case ile:
		return "I"
	case leu:
		return "L"
	case lys:
		return "K"
	case met:
		return "M"
	case phe:
		return "F"
	case ser:
		return "S"
	case thr:
		return "T"
	case trp:
		return "W"
	case tyr:
		return "Y"
	case val:
		return "V"
	case stop:
		return "stop"
	default:
		log.Fatalf("Error: unexpected amino acid value %v\n", a)
		return "X"
	}
}

//TODO: Add in frame as an argument.
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

func TranslateCodon(c *Codon) aa {
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

func TranslateSeq(b []Base) []aa {
	var answer []aa
	codons := BasesToCodons(b)
	for i := 0; i < len(codons); i++ {
		fmt.Printf("%v\n", codons[i])
		answer = append(answer, TranslateCodon(codons[i]))
	}
	return answer
}

func PolypeptideToString(a []aa) string {
	var answer string
	for i := 0; i < len(a); i++ {
		answer = answer + AAToString(a[i])
	}
	return answer
}

func TranslateToString(b []Base) string {
	AllToUpper(b)
	a := TranslateSeq(b)
	for i := 0; i < len(a); i++ {
		fmt.Println(a[i])
	}
	return PolypeptideToString(a)
}
