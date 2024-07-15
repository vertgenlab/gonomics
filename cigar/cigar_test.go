package cigar

import (
	"bytes"
	"log"
	"os"
	"strings"
	"testing"
)

var c1 Cigar = Cigar{RunLength: 35, Op: Match}
var c2 Cigar = Cigar{RunLength: 2, Op: Insertion}
var c3 Cigar = Cigar{RunLength: 16, Op: Deletion}
var c4 Cigar = Cigar{RunLength: 5, Op: Match}
var c5 Cigar = Cigar{RunLength: 2, Op: Mismatch}
var c6 Cigar = Cigar{RunLength: 3, Op: Equal}
var c7 Cigar = Cigar{RunLength: 3, Op: SoftClip}
var unknown []Cigar = []Cigar{{RunLength: 0, Op: Unmapped}}

func TestFromString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
	var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
	var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := FromString(cigarsString)

	for i := 0; i < len(cigarCheck); i++ {
		if !isEqual(cigars[i], cigarCheck[i]) {
			t.Errorf("Error with FromString")
		}
	}
}

func TestToString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
	var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
	var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := ToString(cigars)

	if strings.Compare(cigarCheck, cigarsString) != 0 {
		t.Errorf("Error with ToString")
	}
}

func isEqual(a Cigar, b Cigar) bool {
	return (a.RunLength == b.RunLength && a.Op == b.Op)
}

func TestMatchLength(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected int
	}{
		{[]Cigar{c1, c2, c3}, 35},
		{[]Cigar{c4, c5, c6, c7}, 10},
	}

	for _, c := range tests {
		result := MatchLength(c.input)
		if c.expected != result {
			t.Errorf("Error: MatchLength() %d != %d\n", result, c.expected)
		}
	}
}

func TestReferenceLength(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected int
	}{
		{[]Cigar{c1, c2, c3}, 51},
		{[]Cigar{c4, c5, c6, c7}, 10},
	}

	for _, c := range tests {
		result := ReferenceLength(c.input)
		if c.expected != result {
			t.Errorf("Error: Incorrect ReferenceLength() %d != %d\n", result, c.expected)
		}
	}
	var buf bytes.Buffer
	log.SetOutput(&buf)

	defer func() {
		log.SetOutput(os.Stderr) // Restore default
		r := recover()
		if r == nil {
			log.Printf("Expected ReferenceLength() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumInsertions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	ReferenceLength(unknown)
}

func TestQueryLength(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected int
	}{
		{[]Cigar{c1, c2, c3}, 37},
		{[]Cigar{c4, c5, c6, c7}, 13},
	}

	for _, c := range tests {
		result := QueryLength(c.input)
		if c.expected != result {
			t.Errorf("Error: Incorrect QueryLength() %d != %d\n", result, c.expected)
		}
	}
	var buf bytes.Buffer
	log.SetOutput(&buf)

	defer func() {
		log.SetOutput(os.Stderr) // Restore default
		r := recover()
		if r == nil {
			log.Printf("Expected QueryLength() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumInsertions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	QueryLength(unknown)
}

func TestConsumesReference(t *testing.T) {
	expected := 51
	var result int
	for _, c := range []Cigar{c1, c2, c3} {
		if ConsumesReference(c.Op) {
			result = result + c.RunLength
		}
	}
	if expected != result {
		t.Errorf("Error: Incorrect ConsumesReference() %d != %d\n", result, expected)
	}
}

func TestConsumesQuery(t *testing.T) {
	expected := 37
	var result int
	for _, c := range []Cigar{c1, c2, c3} {
		if ConsumesQuery(c.Op) {
			result = result + c.RunLength
		}
	}
	if expected != result {
		t.Errorf("Error: IncorrectConsumesQuery() %d != %d\n", result, expected)
	}
}

func BenchmarkCigarToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ToString([]Cigar{c1, c2, c3})
	}
}

func BenchmarkCigarFromString(b *testing.B) {
	var cigarsString = "35M2I16D"
	b.ReportAllocs()
	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		FromString(cigarsString)
	}
}
