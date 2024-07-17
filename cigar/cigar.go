// Package cigar contains functions to manipulate cigar data in the SAM file format.
// More information on cigars can be found in http://samtools.github.io/hts-specs/SAMv1.pdf
package cigar

import (
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

// Cigar contains information on the runLength, operation, and DNA sequence associated with a particular cigar character.
type Cigar struct {
	RunLength int
	Op        byte
}

// Cigar operation runes as defined in the SAM format specification to represent the different alignment sequences
const (
	Match     byte = 'M' // Alignment match (can be a sequence match or mismatch)
	Insertion byte = 'I' // Insertion to the reference
	Deletion  byte = 'D' // Deletion from the reference
	Ns        byte = 'N' // Skipped region from the reference
	SoftClip  byte = 'S' // Soft clipping (clipped sequences present in SEQ)
	HardClip  byte = 'H' // Hard clipping (clipped sequences NOT present in SEQ)
	Padded    byte = 'P' // Padding (silent deletion from padded reference)
	Equal     byte = '=' // Sequence match
	Mismatch  byte = 'X' // Sequence mismatch
	Unmapped  byte = '*' // Represents an unmapped read in a SAM file.
)

// OpTable is a slice containing the runes representing CIGAR operation types
var OpTable []byte = []byte{Match, Insertion, Deletion, Ns, SoftClip, HardClip, Padded, Equal, Mismatch}

// Uint32Table is a map that provides a numeric representation of each CIGAR operation byte
var Uint32Table = map[byte]uint32{
	Match:     0,
	Insertion: 1,
	Deletion:  2,
	Ns:        3,
	SoftClip:  4,
	HardClip:  5,
	Padded:    6,
	Equal:     7,
	Mismatch:  8,
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
	if input[0].Op == Unmapped {
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
	var buf strings.Builder
	var err error

	if len(cigars) == 0 || cigars == nil {
		exception.PanicOnErr(buf.WriteByte(Unmapped))
		return buf.String()
	}

	for _, c := range cigars {
		if c.Op == Unmapped {
			exception.PanicOnErr(buf.WriteByte(Unmapped))
			return buf.String()
		}
		_, err = buf.WriteString(strconv.Itoa(c.RunLength))
		exception.PanicOnErr(err)
		exception.PanicOnErr(buf.WriteByte(c.Op))
	}
	return buf.String()
}

// FromString parses an input string into a slice of Cigar structs.
func FromString(input string) []Cigar {
	if input == "*" || input == "**" {
		return []Cigar{{RunLength: 0, Op: Unmapped}}
	}
	var ans []Cigar = make([]Cigar, 0, 1)
	var lastNum int = 0
	for i := 0; i < len(input); i++ {
		if validOp(input[i]) {
			ans = append(ans, Cigar{RunLength: parse.StringToInt(input[lastNum:i]), Op: input[i]})
			lastNum = i + 1
		}
	}
	return ans
}

// MatchLength returns the number of bases in a Cigar slice that align to the reference.
func MatchLength(c []Cigar) int {
	var ans int
	if c[0].Op == Unmapped {
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
func validOp(b byte) bool {
	switch b {
	case Match, Insertion, Deletion, Ns, SoftClip, HardClip, Padded, Equal, Mismatch:
		return true
	default:
		return false
	}
}

// ConsumesReference returns true of the byte matches an operation character that is reference consuming for Cigars.
func ConsumesReference(r byte) bool {
	switch r {
	case Match, Deletion, Ns, Equal, Mismatch:
		return true
	case Insertion, SoftClip, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid byte: %c", r)
	return false
}

// ConsumesQuery returns true for input runes that match query consuming characters for Cigars.
func ConsumesQuery(r byte) bool {
	switch r {
	case Match, Insertion, SoftClip, Equal, Mismatch:
		return true
	case Deletion, Ns, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid byte: %c", r)
	return false
}
