package dna

import (
	"errors"
	"log"
	"strings"
)

// Codon is an array of three DNA bases for genetic analysis of proteins and amino acids.
type Codon [3]Base

// AminoAcid converts the twenty canonical amino acids and stop codon into bytes.
type AminoAcid byte

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

// GeneticCode is a map of codon arrays to amino acids. Used for translating coding sequences to protein sequences.
var GeneticCode = map[Codon]AminoAcid{
	{T, G, A}: Stop, {T, A, A}: Stop, {T, A, G}: Stop,
	{G, T, A}: Val, {G, T, C}: Val, {G, T, G}: Val, {G, T, T}: Val,
	{T, A, T}: Tyr, {T, A, C}: Tyr,
	{T, G, G}: Trp,
	{A, C, A}: Thr, {A, C, G}: Thr, {A, C, T}: Thr, {A, C, C}: Thr,
	{T, C, A}: Ser, {T, C, C}: Ser, {T, C, G}: Ser, {T, C, T}: Ser, {A, G, T}: Ser, {A, G, C}: Ser,
	{C, C, C}: Pro, {C, C, T}: Pro, {C, C, A}: Pro, {C, C, G}: Pro,
	{T, T, T}: Phe, {T, T, C}: Phe,
	{A, T, G}: Met,
	{A, A, A}: Lys, {A, A, G}: Lys,
	{T, T, A}: Leu, {T, T, G}: Leu, {C, T, C}: Leu, {C, T, G}: Leu, {C, T, A}: Leu, {C, T, T}: Leu,
	{A, T, T}: Ile, {A, T, C}: Ile, {A, T, A}: Ile,
	{C, A, T}: His, {C, A, C}: His,
	{G, G, G}: Gly, {G, G, A}: Gly, {G, G, T}: Gly, {G, G, C}: Gly,
	{G, A, A}: Glu, {G, A, G}: Glu,
	{C, A, A}: Gln, {C, A, G}: Gln,
	{T, G, T}: Cys, {T, G, C}: Cys,
	{G, A, T}: Asp, {G, A, C}: Asp,
	{A, A, T}: Asn, {A, A, C}: Asn,
	{A, G, A}: Arg, {A, G, G}: Arg, {C, G, C}: Arg, {C, G, G}: Arg, {C, G, A}: Arg, {C, G, T}: Arg,
	{G, C, A}: Ala, {G, C, G}: Ala, {G, C, T}: Ala, {G, C, C}: Ala,
}

var (
	ErrLenInputSeqNotDivThree = errors.New("length of input sequence is not a factor of three. remaining bases were ignored")
)

// aaToByte is an efficient lookup for the rune corresponding to a given amino acid.
// intended to remain as a private array to help the AminoAcidToShortString function.
// panics if value input is not a valid AminoAcid.
var aaToByte = []byte{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*'}

// aaToLongString is an efficient lookup for the 3 letter string corresponding to a given amino acid.
// intended to remain as a private array to help the AminoAcidToString function.
// panics if value input is not a valid AminoAcid.
var aaToLongString = []string{"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"}

// AminoAcidToShortString converts type AminoAcid into single character amino acid symbols.
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
		log.Panicf("Error: unexpected amino acid value %v\n", a)
		return "X"
	}
}

// AminoAcidToString converts type AminoAcid into three letter amino acid symbols.
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
		log.Panicf("unrecognized amino acid: %v", a)
		return ""
	}
}

// OneLetterToAminoAcid converts a one letter amino acid byte into an AminoAcid type.
func OneLetterToAminoAcid(b byte) AminoAcid {
	switch b {
	case 'A':
		return Ala
	case 'R':
		return Arg
	case 'N':
		return Asn
	case 'D':
		return Asp
	case 'C':
		return Cys
	case 'Q':
		return Gln
	case 'E':
		return Glu
	case 'G':
		return Gly
	case 'H':
		return His
	case 'I':
		return Ile
	case 'L':
		return Leu
	case 'K':
		return Lys
	case 'M':
		return Met
	case 'F':
		return Phe
	case 'P':
		return Pro
	case 'S':
		return Ser
	case 'T':
		return Thr
	case 'W':
		return Trp
	case 'Y':
		return Tyr
	case 'V':
		return Val
	case '*':
		return Stop
	default:
		log.Panicf("cannot convert unrecognized byte %v to type AminoAcid", b)
		return Stop
	}
}

// ThreeLetterToAminoAcid converts a three letter amino acid string into an AminoAcid type
func ThreeLetterToAminoAcid(s string) AminoAcid {
	switch s {
	case "Ala":
		return Ala
	case "Arg":
		return Arg
	case "Asn":
		return Asn
	case "Asp":
		return Asp
	case "Cys":
		return Cys
	case "Gln":
		return Gln
	case "Glu":
		return Glu
	case "Gly":
		return Gly
	case "His":
		return His
	case "Ile":
		return Ile
	case "Leu":
		return Leu
	case "Lys":
		return Lys
	case "Met":
		return Met
	case "Phe":
		return Phe
	case "Pro":
		return Pro
	case "Ser":
		return Ser
	case "Thr":
		return Thr
	case "Trp":
		return Trp
	case "Tyr":
		return Tyr
	case "Val":
		return Val
	case "Ter":
		return Stop
	default:
		log.Panicf("cannot convert unrecognized string %s into type AminoAcid", s)
		return Stop
	}
}

// BasesToCodons converts a slice of DNA bases into a slice of Codons.
// Input expects bases to be in-frame. If the input sequence is not a
// factor of three the function will panic.
func BasesToCodons(b []Base) []Codon {
	frame := len(b) % 3
	answer := make([]Codon, 0, len(b)/3)
	if !IsSeqOfACGT(b) {
		log.Panicf("unrecognized base in sequence: %v, all input bases must be uppercase. N's are not allowed", b)
	}
	if frame != 0 {
		log.Panicf("input sequence length is %d. input length must be a factor of 3", len(b))
	}
	for i := 0; i < len(b); i += 3 {
		answer = append(answer, Codon{b[i], b[i+1], b[i+2]})
	}
	return answer
}

// BasesToCodonsIgnoreRemainder converts a slice of DNA bases into a slice of Codons.
// Any bases remaining after all 3-base codons have been assembled will be ignored.
func BasesToCodonsIgnoreRemainder(b []Base) []Codon {
	frame := len(b) % 3
	return BasesToCodons(b[:len(b)-frame])
}

// CodonsToBases converts a slice of Codons into a slice of DNA bases.
func CodonsToBases(c []Codon) []Base {
	answer := make([]Base, 0, len(c)*3)
	for i := range c {
		answer = append(answer, c[i][0], c[i][1], c[i][2])
	}
	return answer
}

//TranslateCodon converts an individual Codon into the corresponding AminoAcid type.
func TranslateCodon(c Codon) AminoAcid {
	value, found := GeneticCode[c]
	if !found {
		log.Panicf("Error: was not able to translate the codon: %s\n", BasesToString(CodonsToBases([]Codon{c})))
	}
	return value
}

// NonSynonymous compares two Codons and returns true if they encode different AminoAcids.
func NonSynonymous(c1 Codon, c2 Codon) bool {
	return TranslateCodon(c1) != TranslateCodon(c2)
}

// Synonymous compares two codons and returns true if the codons code for the same amino acid.
func Synonymous(c1 Codon, c2 Codon) bool {
	return TranslateCodon(c1) == TranslateCodon(c2)
}

// IsEqual compares two Codons and returns true if the underlying sequences are identical.
func IsEqual(c1 Codon, c2 Codon) bool {
	return c1 == c2
}

// TranslateSeq takes a sequence of DNA bases and translates it into a slice of Amino acids.
// Input expects bases to be in-frame. If the input sequence is not a
// factor of three the function will panic.
func TranslateSeq(b []Base) []AminoAcid {
	answer := make([]AminoAcid, len(b)/3)
	codons := BasesToCodons(b)
	for i := range codons {
		answer[i] = TranslateCodon(codons[i])
	}
	return answer
}

// TranslateSeqToTer takes a sequence of DNA bases and translates it into a slice of Amino acids.
// The translation will end after the first stop codon is reached and the function will return
// the protein sequence including the trailing stop codon. Any bases beyond the stop codon,
// or remaining bases after all 3-base codons have been made will be ignored.
func TranslateSeqToTer(b []Base) []AminoAcid {
	answer := make([]AminoAcid, 0, len(b)/3)
	codons := BasesToCodonsIgnoreRemainder(b)
	for i := range codons {
		answer = append(answer, TranslateCodon(codons[i]))
		if answer[i] == Stop {
			return answer
		}
	}
	return answer
}

// PeptideToShortString converts a slice of amino acid into a string of one character amino acid symbols.
func PeptideToShortString(a []AminoAcid) string {
	var s strings.Builder
	s.Grow(len(a))
	for i := range a {
		s.WriteByte(aaToByte[a[i]])
	}
	return s.String()
}

// PeptideToString converts a slice of AminoAcids into a string of three character amino acid symbols.
func PeptideToString(a []AminoAcid) string {
	var s strings.Builder
	s.Grow(len(a) * 3)
	for i := range a {
		s.WriteString(aaToLongString[a[i]])
	}
	return s.String()
}

// TranslateToShortString converts a sequence of DNA bases into a string of one character amino acid symbols.
// Input expects bases to be in-frame. If the input sequence is not a
// factor of three the function will panic.
func TranslateToShortString(b []Base) string {
	bCopy := make([]Base, len(b))
	copy(bCopy, b)
	AllToUpper(bCopy)
	a := TranslateSeq(bCopy)
	return PeptideToShortString(a)
}

// TranslateToString converts a sequence of DNA bases into a string of three character amino acid symbols.
// Input expects bases to be in-frame. If the input sequence is not a
// factor of three the function will panic.
func TranslateToString(b []Base) string {
	bCopy := make([]Base, len(b))
	copy(bCopy, b)
	AllToUpper(bCopy)
	a := TranslateSeq(bCopy)
	answer := PeptideToString(a)
	return answer
}

// StringToAminoAcid converts a string into type amino acid.
// If singleLetter is false, the input string will be processed by the three letter code.
func StringToAminoAcid(s string, singleLetter bool) []AminoAcid {
	var answer []AminoAcid
	if singleLetter { // single letter amino acid parsing
		answer = make([]AminoAcid, len(s))
		for i := range s {
			answer[i] = OneLetterToAminoAcid(s[i])
		}
	} else { // triple letter amino acid parsing
		answer = make([]AminoAcid, len(s)/3)
		remainder := len(s) % 3
		if remainder != 0 {
			log.Panicf("length of input string is not a factor of three")
		}
		for i := 0; i < len(s)-remainder; i += 3 {
			answer[i/3] = ThreeLetterToAminoAcid(s[i : i+3])
		}
	}
	return answer
}
