package vcf

import (
	"testing"
)

var expectedFormat = map[string]interface{}{
	"FormatF": [][]int{{1}, {2}},
	"FormatJ": [][]rune{{'F', 'E'}, {'D', 'A'}},
	"FormatK": [][]string{{"the,quick,brown,fox"}, {"jumps,over,the,lazy,dog"}},
}

var expectedInfo = map[string]interface{}{
	"InfoA":      [][]int{{10}},
	"InfoB":      [][]float64{{1, 2}},
	"InfoChar":   [][]rune{{'R', 'L'}},
	"InfoFlag":   true,
	"InfoString": [][]string{{"Howdy"}},
}

func TestParseFormat(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	if !equalMap(v.parsedFormat, expectedFormat) {
		t.Errorf("error parsing format")
	}
}

func TestParseInfo(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseInfo(data[0], header)
	if !equalMap(v.parsedInfo, expectedInfo) {
		t.Errorf("error parsing info")
	}
}

func equalMap(a, b map[string]interface{}) bool {
	if len(a) != len(b) {
		return false
	}

	for key := range a {
		switch a[key].(type) {
		case [][]int:
			if !equalInt(a[key], b[key]) {
				return false
			}
		case [][]string:
			if !equalString(a[key], b[key]) {
				return false
			}
		case [][]rune:
			if !equalRune(a[key], b[key]) {
				return false
			}
		case [][]float64:
			if !equalFloat(a[key], b[key]) {
				return false
			}
		case bool:
			if !equalBool(a[key], b[key]) {
				return false
			}
		}
	}
	return true
}

func equalInt(a, b interface{}) bool {
	var ok bool
	var aInt, bInt [][]int
	aInt, ok = a.([][]int)
	if !ok {
		return false
	}
	bInt, ok = b.([][]int)
	if !ok {
		return false
	}
	if len(aInt) != len(bInt) {
		return false
	}

	for i := range aInt {
		if len(aInt[i]) != len(bInt[i]) {
			return false
		}
		for j := range aInt[i] {
			if aInt[i][j] != bInt[i][j] {
				return false
			}
		}
	}
	return true
}

func equalString(a, b interface{}) bool {
	var ok bool
	var aString, bString [][]string
	aString, ok = a.([][]string)
	if !ok {
		return false
	}
	bString, ok = b.([][]string)
	if !ok {
		return false
	}
	if len(aString) != len(bString) {
		return false
	}

	for i := range aString {
		if len(aString[i]) != len(bString[i]) {
			return false
		}
		for j := range aString[i] {
			if aString[i][j] != bString[i][j] {
				return false
			}
		}
	}
	return true
}

func equalRune(a, b interface{}) bool {
	var ok bool
	var aRune, bRune [][]rune
	aRune, ok = a.([][]rune)
	if !ok {
		return false
	}
	bRune, ok = b.([][]rune)
	if !ok {
		return false
	}
	if len(aRune) != len(bRune) {
		return false
	}

	for i := range aRune {
		if len(aRune[i]) != len(bRune[i]) {
			return false
		}
		for j := range aRune[i] {
			if aRune[i][j] != bRune[i][j] {
				return false
			}
		}
	}
	return true
}

func equalFloat(a, b interface{}) bool {
	var ok bool
	var aFloat, bFloat [][]float64
	aFloat, ok = a.([][]float64)
	if !ok {
		return false
	}
	bFloat, ok = b.([][]float64)
	if !ok {
		return false
	}
	if len(aFloat) != len(bFloat) {
		return false
	}

	for i := range aFloat {
		if len(aFloat[i]) != len(bFloat[i]) {
			return false
		}
		for j := range aFloat[i] {
			if aFloat[i][j] != bFloat[i][j] {
				return false
			}
		}
	}
	return true
}

func equalBool(a, b interface{}) bool {
	var ok bool
	var aBool, bBool bool
	aBool, ok = a.(bool)
	if !ok {
		return false
	}
	bBool, ok = b.(bool)
	if !ok {
		return false
	}
	if aBool != bBool {
		return false
	}
	return true
}

func TestQueryInt(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	v = ParseInfo(v, header)

	if !equalInt(QueryInt(v, header.Info["InfoA"].Key), expectedInfo["InfoA"]) {
		t.Errorf("troble querying read info")
	}
	if !equalInt(QueryInt(v, header.Format["FormatF"].Key), expectedFormat["FormatF"]) {
		t.Errorf("troble querying read info")
	}
}

func TestQueryRune(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	v = ParseInfo(v, header)

	if !equalRune(QueryRune(v, header.Info["InfoChar"].Key), expectedInfo["InfoChar"]) {
		t.Errorf("troble querying read info")
	}
	if !equalRune(QueryRune(v, header.Format["FormatJ"].Key), expectedFormat["FormatJ"]) {
		t.Errorf("troble querying read info")
	}
}

func TestQueryString(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	v = ParseInfo(v, header)

	if !equalString(QueryString(v, header.Info["InfoString"].Key), expectedInfo["InfoString"]) {
		t.Errorf("troble querying read info")
	}
	if !equalString(QueryString(v, header.Format["FormatK"].Key), expectedFormat["FormatK"]) {
		t.Errorf("troble querying read info")
	}
}

func TestQueryFloat(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	v = ParseInfo(v, header)

	if !equalFloat(QueryFloat(v, header.Info["InfoB"].Key), expectedInfo["InfoB"]) {
		t.Errorf("troble querying read info")
	}
}

func TestQueryFlag(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	v = ParseInfo(v, header)

	if !equalBool(QueryFlag(v, header.Info["InfoFlag"].Key), expectedInfo["InfoFlag"]) {
		t.Errorf("troble querying read info")
	}
}
