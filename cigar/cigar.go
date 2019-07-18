package cigar

import(
	"log"
	"github.com/vertgenlab/gonomics/common"
	"fmt"
	"unicode"
)

type Cigar struct {
	RunLength int64
	Op        rune
}

func NumInsertions(input []*Cigar) int64 {
	var count int64
	if input[0].Op == '*' {
		log.Fatalf("Cannot calculate NumInsertions from unaligned reads.")
	}
	for i := 0; i < len(input); i++ {
		if !ConsumesReference(input[i].Op) && ConsumesQuery(input[i].Op) {
			count = count + input[i].RunLength
		}
	}
	return count
}

func NumDeletions(input []*Cigar) int64 {
	var count int64
	if input[0].Op == '*' {
		log.Fatalf("Cannot calculate NumInsertions from unaligned reads.")
	}
	for i := 0; i < len(input); i++ {
		if ConsumesReference(input[i].Op) && !ConsumesQuery(input[i].Op) {
			count = count + input[i].RunLength
		}
	}
	return count
}

func ToString(c []*Cigar) string {
	var output string = ""
	for _, v := range c {
		if v.Op == '*' {
			output = "*"
			break
		}
		output = output + fmt.Sprintf("%v%c", v.RunLength, v.Op)
	}
	return output
}

func FromString(input string) []*Cigar {
	var output []*Cigar
	var currentNumber string
	if input == "*" || input == "**" {
		currentCigar := Cigar{RunLength: 0, Op: '*'}
		return append(output, &currentCigar)
	}

	for _, v := range input {
		if unicode.IsDigit(v) {
			currentNumber = currentNumber + fmt.Sprintf("%c", v)	
		} else if RuneIsValidCharacter(v) {
			currentCigar := Cigar{RunLength: common.StringToInt64(currentNumber), Op: v}
			output = append(output, &currentCigar)
			currentNumber = ""
		} else {
			log.Fatalf("Invalid character: %c", v)
		}
	}
/*
	if RuneIsDigit(v) { //will not work, v is not seen outside loop
		log.Fatalf("Cigar ended with digit")
	} */
	return output
}

func ReferenceLength(c []*Cigar) int64 {
	var ans int64
	if c[0].Op == '*' {
		log.Fatalf("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesReference(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

func QueryLength(c []*Cigar) int64 {
	var ans int64
	if c[0].Op == '*' {
		log.Fatalf("Cannot calculate NumInsertions from unaligned reads.")
	}
	for _, v := range c {
		if ConsumesQuery(v.Op) {
			ans = ans + v.RunLength
		}
	}
	return ans
}

/*
func RuneIsDigit(r rune) bool {
	switch r {
	case '0':
		return true
	case '1':
		return true
	case '2':
		return true
	case '3':
		return true
	case '4':
		return true
	case '5':
		return true	
	case '6':
		return true
	case '7':
		return true
	case '8':
		return true
	case '9':
		return true		
	}
	return false
}
*/

func RuneIsValidCharacter(r rune) bool {
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
	log.Fatalf("Invalid rune: %c", r)
	return false
}

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
	log.Fatalf("Invalid rune: %c", r)
	return false
}
