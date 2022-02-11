// Package cigar contains functions to manipulate cigar data in the SAM file format.
//More information on cigars can be found in http://samtools.github.io/hts-specs/SAMv1.pdf

package cigar

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
	"unicode"
)

// The Cigar struct contains information on the runLength, operation, and DNA sequence associated with a particular cigar character.
type Cigar struct {
	RunLength int
	Op        rune
	Sequence  []dna.Base
}

//NumInsertions calculates the number of inserted bases relative to a reference genome for an input Cigar slice.
func NumInsertions(input []Cigar) int {
	var count int
	if input[0].Op == '*' {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for i := 0; i < len(input); i++ {
		if !ConsumesReference(input[i].Op) && ConsumesQuery(input[i].Op) {
			count = count + input[i].RunLength
		}
	}
	return count
}

//NumDeletions calculates the number of deletions relative to a reference genome for an input Cigar slice.
func NumDeletions(input []Cigar) int {
	var count int
	if input[0].Op == '*' {
		log.Panic("Cannot calculate NumDeletions from unaligned reads.")
	}
	for i := 0; i < len(input); i++ {
		if ConsumesReference(input[i].Op) && !ConsumesQuery(input[i].Op) {
			count = count + input[i].RunLength
		}
	}
	return count
}

//ToString converts a slice of Cigar structs to a string for producing readable outputs for files or standard out.
func ToString(c []Cigar) string {
	if len(c) == 0 {
		return "*"
	}
	var output string = ""
	for _, v := range c {
		if v.Op == '*' {
			output = "*"
			break
		}
		printSeq := dna.BasesToString(v.Sequence)
		output = output + fmt.Sprintf("%v%c%s", v.RunLength, v.Op, strings.ToLower(printSeq))
	}
	return output
}

//FromString parses an input string into a slice of Cigar structs.
func FromString(input string) []Cigar {
	var output []Cigar
	var currentNumber string
	var currentCigar Cigar
	if input == "*" || input == "**" {
		currentCigar = Cigar{RunLength: 0, Op: '*'}
		return append(output, currentCigar)
	}

	for _, v := range input {
		if unicode.IsDigit(v) {
			currentNumber = currentNumber + fmt.Sprintf("%c", v)
		} else if validOp(v) {
			currentCigar := Cigar{RunLength: common.StringToInt(currentNumber), Op: v}
			output = append(output, currentCigar)
			currentNumber = ""
		} else {
			log.Panicf("Invalid character: %c", v)
		}
	}
	return output
}

//MatchLength returns the number of bases in a Cigar slice that align to the reference.
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

//ReferenceLength calculates the number of reference positions that a Cigar slice spans.
func ReferenceLength(c []Cigar) int {
	var ans int
	if c[0].Op == '*' {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesReference(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

//QueryLength calculates the length of the query read from a slice of Cigar structs.
func QueryLength(c []Cigar) int {
	var ans int
	if c[0].Op == '*' {
		log.Panic("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesQuery(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

//validOp returns true if a particular input rune matches any of the acceptable Cigar operation characters.
func validOp(r rune) bool {
	switch r {
	case 'M':
		return true
	case 'I':
		return true
	case 'D':
		return true
	case 'N':
		return true
	case 'S':
		return true
	case 'H':
		return true
	case 'P':
		return true
	case '=':
		return true
	case 'X':
		return true
	}
	return false
}

//CigarConsumesReference returns true if the Cigar operation is reference consuming, false otherwise.
func CigarConsumesReference(c Cigar) bool {
	return ConsumesReference(c.Op)
}

//ConsumesReference returns true of the rune matches an operation character that is reference consuming for Cigars.
func ConsumesReference(r rune) bool {
	switch r {
	case 'M':
		return true
	case 'I':
		return false
	case 'D':
		return true
	case 'N':
		return true
	case 'S':
		return false
	case 'H':
		return false
	case 'P':
		return false
	case '=':
		return true
	case 'X':
		return true
	}
	log.Panicf("Invalid rune: %c", r)
	return false
}

//ConsumesQuery returns true for input runes that match query consuming characters for Cigars.
func ConsumesQuery(r rune) bool {
	switch r {
	case 'M':
		return true
	case 'I':
		return true
	case 'D':
		return false
	case 'N':
		return false
	case 'S':
		return true
	case 'H':
		return false
	case 'P':
		return false
	case '=':
		return true
	case 'X':
		return true
	}
	log.Panicf("Invalid rune: %c", r)
	return false
}
