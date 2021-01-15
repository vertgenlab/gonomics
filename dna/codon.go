package dna

import (
	"errors"
	"strings"
)

// Codon is an array of three DNA bases for genetic analysis of proteins and amino acids.
type Codon [3]Base

// AminoAcid converts the twenty canonical amino acids and stop codon into bytes.
type AminoAcid byte

const (
	Ala  AminoAcid = 'A'
	Arg  AminoAcid = 'R'
	Asn  AminoAcid = 'N'
	Asp  AminoAcid = 'D'
	Cys  AminoAcid = 'C'
	Gln  AminoAcid = 'Q'
	Glu  AminoAcid = 'E'
	Gly  AminoAcid = 'G'
	His  AminoAcid = 'H'
	Ile  AminoAcid = 'I'
	Leu  AminoAcid = 'L'
	Lys  AminoAcid = 'K'
	Met  AminoAcid = 'M'
	Phe  AminoAcid = 'F'
	Pro  AminoAcid = 'P'
	Ser  AminoAcid = 'S'
	Thr  AminoAcid = 'T'
	Trp  AminoAcid = 'W'
	Tyr  AminoAcid = 'Y'
	Val  AminoAcid = 'V'
	Stop AminoAcid = '*'
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
	ErrUnrecognizedAminoAcid  = errors.New("unrecognized amino acid")
	ErrCodonNotInMap          = errors.New("codon is not in dna.GeneticCode map")
	ErrLenInputSeqNotDivThree = errors.New("length of input sequence is not a factor of three. remaining bases were ignored")
)

// AminoAcidToShortString converts type AminoAcid into single character amino acid symbols.
func AminoAcidToShortString(a AminoAcid) string {
	return string(a)
}

// AminoAcidToString converts type AminoAcid into three letter amino acid symbols.
func AminoAcidToString(a AminoAcid) (string, error) {
	switch a {
	case Ala:
		return "Ala", nil
	case Arg:
		return "Arg", nil
	case Asn:
		return "Asn", nil
	case Asp:
		return "Asp", nil
	case Cys:
		return "Cys", nil
	case Gln:
		return "Gln", nil
	case Glu:
		return "Glu", nil
	case Gly:
		return "Gly", nil
	case His:
		return "His", nil
	case Ile:
		return "Ile", nil
	case Leu:
		return "Leu", nil
	case Lys:
		return "Lys", nil
	case Met:
		return "Met", nil
	case Phe:
		return "Phe", nil
	case Pro:
		return "Pro", nil
	case Ser:
		return "Ser", nil
	case Thr:
		return "Thr", nil
	case Trp:
		return "Trp", nil
	case Tyr:
		return "Tyr", nil
	case Val:
		return "Val", nil
	case Stop:
		return "Ter", nil
	default:
		return "X", ErrUnrecognizedAminoAcid
	}
}

// BasesToCodons converts a slice of DNA bases into a slice of Codons.
// Input expects bases to be in-frame. If the input sequence is not a
// factor of three, the return codons will not include the remaining
// bases and will also return an error.
func BasesToCodons(b []Base) ([]*Codon, error) {
	var err error
	frame := len(b) % 3
	answer := make([]*Codon, 0, len(b)/3)

	if frame != 0 {
		err = ErrLenInputSeqNotDivThree
	}
	for i := 0; i < len(b)-frame; i += 3 {
		answer = append(answer, &Codon{b[i], b[i+1], b[i+2]})
	}
	return answer, err
}

// CodonsToSeq converts a slice of Codons into a slice of DNA bases.
func CodonsToSeq(c []*Codon) []Base {
	answer := make([]Base, 0, len(c)*3)
	for i := range c {
		answer = append(answer, (*c[i])[0], (*c[i])[1], (*c[i])[2])
	}
	return answer
}

//TranslateCodon converts an individual Codon into the corresponding AminoAcid type.
func TranslateCodon(c *Codon) (AminoAcid, error) {
	val, ok := GeneticCode[*c]
	if ok {
		return val, nil
	} else {
		return val, ErrCodonNotInMap
	}
}

// NonSynonymous compares two Codons and returns true if they encode different AminoAcids.
func NonSynonymous(c1 *Codon, c2 *Codon) bool {
	p1, _ := TranslateCodon(c1)
	p2, _ := TranslateCodon(c2)
	return p1 != p2
}

// Synonymous compares two codons and returns true if the codons code for the same amino acid.
func Synonymous(c1 *Codon, c2 *Codon) bool {
	p1, _ := TranslateCodon(c1)
	p2, _ := TranslateCodon(c2)
	return p1 == p2
}

// IsEqual compares two Codons and returns true if the underlying sequences are identical.
func IsEqual(c1 *Codon, c2 *Codon) bool {
	return *c1 == *c2
}

// TranslateSeq takes a sequence of DNA bases and translates it into a slice of Amino acids.
func TranslateSeq(b []Base) ([]AminoAcid, error) {
	answer := make([]AminoAcid, len(b)/3)
	codons, err := BasesToCodons(b)
	for i := range codons {
		answer[i], err = TranslateCodon(codons[i])
	}
	return answer, err
}

// PolypeptideToShortString converts a slice of amino acid into a string of one character amino acid symbols.
func PolypeptideToShortString(a []AminoAcid) string {
	var s strings.Builder
	s.Grow(len(a))
	for i := range a {
		s.WriteByte(byte(a[i]))
	}
	return s.String()
}

// PolypeptideToString converts a slice of AminoAcids into a string of three character amino acid symbols.
func PolypeptideToString(a []AminoAcid) (string, error) {
	var s strings.Builder
	s.Grow(len(a) * 3)
	var str string
	var err, finalErr error
	for i := range a {
		str, err = AminoAcidToString(a[i])
		if err != nil {
			finalErr = err
		}
		s.WriteString(str)
	}
	return s.String(), finalErr
}

//TranslateToShortString converts a sequence of DNA bases into a string of one character amino acid symbols.
func TranslateToShortString(b []Base) (string, error) {
	AllToUpper(b)
	a, err := TranslateSeq(b)
	return PolypeptideToShortString(a), err
}

// TranslateToString converts a sequence of DNA bases into a string of three character amino acid symbols.
func TranslateToString(b []Base) (string, error) {
	var finalErr error
	AllToUpper(b)
	a, err := TranslateSeq(b)
	if err != nil {
		finalErr = err
	}
	answer, err := PolypeptideToString(a)
	if err != nil {
		finalErr = err
	}
	return answer, finalErr
}

// ShortStringToPolypeptide converts a string into a slice of amino acids.
func ShortStringToPolypeptide(s string) []AminoAcid {
	answer := make([]AminoAcid, len(s))
	for i := range s {
		answer[i] = AminoAcid(s[i])
	}
	return answer
}
