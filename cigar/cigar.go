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
	Op        rune
}

// Cigar operation runes as defined in the SAM format specification to represent the different alignment sequences
const (
	Match     rune = 'M' // Alignment match (can be a sequence match or mismatch)
	Insertion rune = 'I' // Insertion to the reference
	Deletion  rune = 'D' // Deletion from the reference
	Ns        rune = 'N' // Skipped region from the reference
	SoftClip  rune = 'S' // Soft clipping (clipped sequences present in SEQ)
	HardClip  rune = 'H' // Hard clipping (clipped sequences NOT present in SEQ)
	Padded    rune = 'P' // Padding (silent deletion from padded reference)
	Equal     rune = '=' // Sequence match
	Mismatch  rune = 'X' // Sequence mismatch
	Unmapped  rune = '*' // Represents an unmapped read in a SAM file.
)

// OpTable is a slice containing the runes representing CIGAR operation types
var OpTable []rune = []rune{Match, Insertion, Deletion, Ns, SoftClip, HardClip, Padded, Equal, Mismatch}

// Uint32Table is a map that provides a numeric representation of each CIGAR operation rune
var Uint32Table = map[rune]uint32{
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
		_, err = buf.WriteRune(Unmapped)
		exception.PanicOnErr(err)
		return buf.String()
	}

	for _, c := range cigars {
		if c.Op == Unmapped {
			_, err = buf.WriteRune(c.Op)
			exception.PanicOnErr(err)
			return buf.String()
		}
		_, err = buf.WriteString(strconv.Itoa(c.RunLength))
		exception.PanicOnErr(err)
		_, err = buf.WriteRune(c.Op)
		exception.PanicOnErr(err)
	}
	return buf.String()
}

// FromString parses an input string into a slice of Cigar structs.
func FromString(input string) []Cigar {
	if input == "*" || input == "**" {
		return []Cigar{{RunLength: 0, Op: Unmapped}}
	}

	var output []Cigar
	var currentNumber strings.Builder
	currentNumber.Grow(4)

	for _, r := range input {
		if unicode.IsDigit(r) {
			currentNumber.WriteRune(r)
		} else if validOp(r) {

			output = append(output, Cigar{RunLength: parse.StringToInt(currentNumber.String()), Op: r})
			currentNumber.Reset()
		} else {
			exception.PanicOnErr(fmt.Errorf("invalid character: %c", r)) // Panic on invalid char
		}
	}
	// Handle case where the input ends with a number (incomplete CIGAR)
	if currentNumber.Len() > 0 {
		output = append(output, Cigar{RunLength: parse.StringToInt(currentNumber.String())}) // Assume 'M' if op is missing
	}

	return output
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

// validOp returns true if a particular input rune matches any of the acceptable Cigar operation characters.
func validOp(r rune) bool {
	switch r {
	case Match, Insertion, Deletion, Ns, SoftClip, HardClip, Padded, Equal, Mismatch:
		return true
	default:
		return false
	}
}

// ConsumesReference returns true of the rune matches an operation character that is reference consuming for Cigars.
func ConsumesReference(r rune) bool {
	switch r {
	case Match, Deletion, Ns, Equal, Mismatch:
		return true
	case Insertion, SoftClip, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid rune: %c", r)
	return false
}

// ConsumesQuery returns true for input runes that match query consuming characters for Cigars.
func ConsumesQuery(r rune) bool {
	switch r {
	case Match, Insertion, SoftClip, Equal, Mismatch:
		return true
	case Deletion, Ns, HardClip, Padded:
		return false
	}
	log.Panicf("Invalid rune: %c", r)
	return false
}
