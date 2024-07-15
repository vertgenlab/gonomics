// Package cigar contains functions to manipulate cigar data in the SAM file format.
// More information on cigars can be found in http://samtools.github.io/hts-specs/SAMv1.pdf
package cigar

import (
	"fmt"
	"log"
	"strconv"
	"strings"
	"unicode"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

// Cigar contains information on the runLength, operation, and DNA sequence associated with a particular cigar character.
type Cigar struct {
	RunLength int
	Op        byte
}

// Defined const for byte cigar.
const (
	Match     byte = 'M'
	Insertion byte = 'I'
	Deletion  byte = 'D'
	N         byte = 'N'
	SoftClip  byte = 'S'
	HardClip  byte = 'H'
	Padded    byte = 'P'
	Equal     byte = '='
	Mismatch  byte = 'X'
	Unmapped  byte = '*'
)

var OpTable []byte = []byte{Match, Insertion, Deletion, N, SoftClip, HardClip, Padded, Equal, Mismatch}

var Uint32Table = []uint32{
	Match: 0, Insertion: 1, Deletion: 2, N: 3, SoftClip: 4, HardClip: 5, Padded: 6, Equal: 7, Mismatch: 8,
}

// NumInsertions calculates the number of inserted bases relative to a reference genome for an input Cigar slice.
func NumInsertions(input []Cigar) int {
	var count int
	if input[0].Op == Unmapped {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for i := range input {
		if !ConsumesReference(input[i].Op) && ConsumesQuery(input[i].Op) {
			count += input[i].RunLength
		}
	}
	return count
}

// NumDeletions calculates the number of deletions relative to a reference genome for an input Cigar slice.
func NumDeletions(input []Cigar) int {
	var count int
	if input[0].Op == '*' {
		log.Panic("Cannot calculate NumDeletions from unaligned reads.")
	}
	for i := range input {
		if ConsumesReference(input[i].Op) && !ConsumesQuery(input[i].Op) {
			count += input[i].RunLength
		}
	}
	return count
}

// ToString converts a slice of Cigar structs to a string for producing readable outputs for files or standard out.
func ToString(cigars []Cigar) string {
	if len(cigars) == 0 || cigars == nil {
		return "*"
	}

	var buf strings.Builder
	var err error

	for _, c := range cigars {
		if c.Op == Unmapped {
			err = buf.WriteByte(c.Op)
			exception.PanicOnErr(err)
			break
		}
		_, err = buf.WriteString(strconv.Itoa(int(c.RunLength)))
		exception.PanicOnErr(err)
		err = buf.WriteByte(c.Op)
		exception.PanicOnErr(err)
	}
	return buf.String()
}

// FromString parses an input string into a slice of Cigar structs.
func FromString(input string) []Cigar {
	var output []Cigar
	var currentNumber string
	var currentCigar Cigar
	if input == "*" || input == "**" {
		currentCigar = Cigar{RunLength: 0, Op: Unmapped}
		return append(output, currentCigar)
	}

	for _, v := range input {
		if unicode.IsDigit(v) {
			currentNumber = currentNumber + fmt.Sprintf("%c", v)
		} else if validOp(byte(v)) {
			currentCigar := Cigar{RunLength: parse.StringToInt(currentNumber), Op: byte(v)}
			output = append(output, currentCigar)
			currentNumber = ""
		} else {
			log.Panicf("Invalid character: %c", v)
		}
	}
	return output
}

// MatchLength returns the number of bases in a Cigar slice that align to the reference.
func MatchLength(c []Cigar) int {
	var ans int
	if c[0].Op == '*' {
		log.Panic("Cannot calculate MatchLength from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesReference(v.Op) && ConsumesQuery(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

// ReferenceLength calculates the number of reference positions that a Cigar slice spans.
func ReferenceLength(c []Cigar) int {
	var ans int
	if c[0].Op == Unmapped {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesReference(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

// QueryLength calculates the length of the query read from a slice of Cigar structs.
func QueryLength(c []Cigar) int {
	var ans int
	if c[0].Op == Unmapped {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesQuery(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

// validOp returns true if a particular input byte matches any of the acceptable Cigar operation characters.
func validOp(op byte) bool {
	switch op {
	case Match, Insertion, Deletion, N, SoftClip, HardClip, Padded, Equal, Mismatch:
		return true
	default:
		return false
	}
}

// ConsumesReference returns true of the byte matches an operation character that is reference consuming for Cigars.
func ConsumesReference(op byte) bool {
	switch op {
	case Match, Deletion, N, Equal, Mismatch:
		return true
	case Insertion, SoftClip, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid byte: %c", op)
	return false
}

// ConsumesQuery returns true for input runes that match query consuming characters for Cigars.
func ConsumesQuery(r byte) bool {
	switch r {
	case Match, Insertion, SoftClip, Equal, Mismatch:
		return true
	case Deletion, N, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid byte: %c", r)
	return false
}

// Uint32ToCigar will process a uint32 slice and decode each number into a byte cigar struct.
// CIGAR operation lengths are limited to 2^28-1 in the current sam/bam formats.
func Uint32ToCigar(cigar []uint32) []Cigar {
	var answer []Cigar = make([]Cigar, len(cigar))
	for i := 0; i < len(cigar); i++ {
		answer[i] = Cigar{RunLength: int(cigar[i] >> 4), Op: OpTable[cigar[i]&0xf]}
	}
	return answer
}

// CigarToUint32 will convert a slice of []ByteCigar to a slice of []uint32.
func CigarToUint32(cigar []Cigar) []uint32 {
	var answer []uint32 = make([]uint32, len(cigar))
	for i := 0; i < len(cigar); i++ {
		answer[i] = Uint32Table[cigar[i].Op] | uint32(cigar[i].RunLength)<<4
	}
	return answer
}
