package main

import (
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strconv"
	"strings"
)

type operator string

const (
	equal        operator = "="
	notEqual     operator = "!="
	greater      operator = ">"
	less         operator = "<"
	greaterEqual operator = ">="
	lessEqual    operator = "<="
	present      operator = ""
)

// TODO complex expression parsing
// parseExpression parses a ';' delimited string to a slice of boolean functions to test a vcf record.
func parseExpression(input string, header vcf.Header, isFormatElseInfo bool) testingFuncs {
	var answer testingFuncs
	input = strings.Trim(input, "\"") // trim quotations from the ends
	expSlice := strings.Split(input, ";")
	for _, exp := range expSlice {
		var op operator
		var tagValuePair []string
		var tag, value string
		op = searchOp(exp)
		if op != present {
			tagValuePair = strings.Split(exp, string(op))
		} else {
			tagValuePair = []string{exp}
		}
		tag = strings.Trim(tagValuePair[0], " ")
		if len(tagValuePair) == 2 {
			value = strings.Trim(tagValuePair[1], " ")
		}
		answer = append(answer, getRelationshipTest(tag, value, op, header, isFormatElseInfo))
	}
	return answer
}

func searchOp(exp string) operator {
	if strings.Contains(exp, ">=") {
		return greaterEqual
	}

	if strings.Contains(exp, "<=") {
		return lessEqual
	}

	if strings.Contains(exp, "!=") {
		return notEqual
	}

	if strings.Contains(exp, "=") {
		return equal
	}

	if strings.Contains(exp, ">") {
		return greater
	}

	if strings.Contains(exp, "<") {
		return less
	}

	return present
}

func getRelationshipTest(tag string, value string, r operator, header vcf.Header, isFormatElseInfo bool) func(vcf.Vcf) bool {
	var tagKey vcf.Key
	if isFormatElseInfo {
		tagKey = header.Format[tag].Key
	} else {
		tagKey = header.Info[tag].Key
	}

	if tagKey.Number != "1" && tagKey.Number != "0" {
		log.Printf("WARNING: expressions for tags with multiple values will be true if any value passes filter. Tag '%s' has '%s' fields.", tag, tagKey.Number)
	}

	switch tagKey.DataType {
	case vcf.Integer:
		test := r.TestInteger()
		val, err := strconv.Atoi(value)
		if err != nil {
			log.Fatalf("Error: value '%s' is not an integer as expected for tag '%s'.", value, tag)
		}
		return func(v vcf.Vcf) bool {
			var answer bool = true
			for _, recordVal := range vcf.QueryInt(v, tagKey)[0] {
				if !test(recordVal, val) {
					answer = false
				}
			}
			return answer
		}

	case vcf.Float:
		test := r.TestFloat()
		val, err := strconv.ParseFloat(value, 64)
		if err != nil {
			log.Fatalf("Error: value '%s' is not a float as expected for tag '%s'.", value, tag)
		}
		return func(v vcf.Vcf) bool {
			var answer bool = true
			for _, recordVal := range vcf.QueryFloat(v, tagKey)[0] {
				if !test(recordVal, val) {
					answer = false
				}
			}
			return answer
		}

	case vcf.Character:
		test := r.TestCharacter()
		if len(value) != 1 {
			log.Fatalf("Error: value '%s' is not a character as expected for tag '%s'.", value, tag)
		}
		return func(v vcf.Vcf) bool {
			var answer bool = true
			for _, recordVal := range vcf.QueryRune(v, tagKey)[0] {
				if !test(recordVal, rune(value[0])) {
					answer = false
				}
			}
			return answer
		}

	case vcf.String:
		test := r.TestString()
		return func(v vcf.Vcf) bool {
			var answer bool = true
			for _, recordVal := range vcf.QueryString(v, tagKey)[0] {
				if !test(recordVal, value) {
					answer = false
				}
			}
			return answer
		}

	case vcf.Flag:
		if value != "" {
			log.Fatalf("Error: expression specified a value for '%s' but '%s' is of type flag. Flags must be in the form of '%s' or '!%s'", tag, tag, tag, tag)
		}
		return func(v vcf.Vcf) bool {
			return vcf.QueryFlag(v, tagKey)
		}

	default:
		log.Panic()
		return nil
	}
}

func (r operator) TestInteger() func(a, b int) bool {
	switch r {
	case equal:
		return func(a, b int) bool { return a == b }
	case notEqual:
		return func(a, b int) bool { return a != b }
	case greater:
		return func(a, b int) bool { return a > b }
	case less:
		return func(a, b int) bool { return a < b }
	case greaterEqual:
		return func(a, b int) bool { return a >= b }
	case lessEqual:
		return func(a, b int) bool { return a <= b }
	default:
		log.Panicf("unrecognized operator '%v'", r)
		return nil
	}
}

func (r operator) TestFloat() func(a, b float64) bool {
	switch r {
	case equal:
		return func(a, b float64) bool { return a == b }
	case notEqual:
		return func(a, b float64) bool { return a != b }
	case greater:
		return func(a, b float64) bool { return a > b }
	case less:
		return func(a, b float64) bool { return a < b }
	case greaterEqual:
		return func(a, b float64) bool { return a >= b }
	case lessEqual:
		return func(a, b float64) bool { return a <= b }
	default:
		log.Panicf("unrecognized operator '%v'", r)
		return nil
	}
}

func (r operator) TestString() func(a, b string) bool {
	switch r {
	case equal:
		return func(a, b string) bool { return a == b }
	case notEqual:
		return func(a, b string) bool { return a != b }
	case greater:
		return func(a, b string) bool { return a > b }
	case less:
		return func(a, b string) bool { return a < b }
	case greaterEqual:
		return func(a, b string) bool { return a >= b }
	case lessEqual:
		return func(a, b string) bool { return a <= b }
	default:
		log.Panicf("unrecognized operator '%v'", r)
		return nil
	}
}

func (r operator) TestCharacter() func(a, b rune) bool {
	switch r {
	case equal:
		return func(a, b rune) bool { return a == b }
	case notEqual:
		return func(a, b rune) bool { return a != b }
	case greater:
		return func(a, b rune) bool { return a > b }
	case less:
		return func(a, b rune) bool { return a < b }
	case greaterEqual:
		return func(a, b rune) bool { return a >= b }
	case lessEqual:
		return func(a, b rune) bool { return a <= b }
	default:
		log.Panicf("unrecognized operator '%v'", r)
		return nil
	}
}
