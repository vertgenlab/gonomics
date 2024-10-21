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
var unknown []Cigar = []Cigar{} // len(cigars) == 0 standard for unmapped

func TestIsUnmapped(t *testing.T) {
	tests := []struct {
		cigar    []Cigar
		expected bool
	}{
		{[]Cigar{}, true},
		{[]Cigar{{1, Match}, {1, Equal}, {1, Mismatch}}, false},
		{make([]Cigar, 0), true},
		{nil, true}, // a slice that is nil will currently return true 0xff
		{FromString("*"), true},
		{FromString("0xff"), true},
		{FromString("150M"), false},
		{make([]Cigar, 1), false}, // TODO: Consider the result in which a slice is allocated memory, but contains empty values: i.e Cigar{0, nil}.
	}
	for _, test := range tests {
		result := IsUnmapped(test.cigar)
		if test.expected != result {
			t.Errorf("Error: Incorrect IsUnmapped() result. %s (%t) != %t", ToString(test.cigar), result, test.expected)
		}
	}
}

func TestNumInsertions(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected int
	}{
		{[]Cigar{c1, c2, c3}, 2},
	}

	for _, c := range tests {
		result := NumInsertions(c.input)
		if c.expected != result {
			t.Errorf("Error: Incorrect NumInsertions() %d != %d\n", result, c.expected)
		}
	}
	var buf bytes.Buffer
	log.SetOutput(&buf)

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			log.Printf("Expected NumInsertions() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumInsertions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from unmapped []Cigar
	NumInsertions(unknown)
}

func TestNumDeletions(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected int
	}{
		{[]Cigar{c1, c2, c3}, 16},
	}

	for _, c := range tests {
		result := NumDeletions(c.input)
		if c.expected != result {
			t.Errorf("Error: Incorrect NumDeletions() %d != %d\n", result, c.expected)
		}
	}
	var buf bytes.Buffer
	log.SetOutput(&buf)

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			log.Printf("Expected NumDeletions() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumDeletions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from unmapped []Cigar
	NumDeletions(unknown)
}
func TestFromString(t *testing.T) {
	tests := []struct {
		input    string
		expected []Cigar
	}{
		{"35M2I16D", []Cigar{c1, c2, c3}},
		{"*", []Cigar{}},
		{"*", make([]Cigar, 0)},
	}
	for _, test := range tests {
		result := FromString(test.input)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: Incorrect FromString() %s != %s\n", ToString(result), ToString(test.expected))
		}
	}
}

func TestToString(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected string
	}{
		{[]Cigar{c1, c2, c3}, "35M2I16D"},
		{make([]Cigar, 0), "*"},
		{[]Cigar{}, "*"},
	}
	for _, test := range tests {
		result := ToString(test.input)
		if strings.Compare(test.expected, result) != 0 {
			t.Errorf("Error: Incorrect ToString() %s != %s\n", result, test.expected)
		}
	}
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
	var buf bytes.Buffer
	log.SetOutput(&buf)

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			log.Printf("Expected MatchLength() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate MatchLength from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from unmapped []Cigar
	MatchLength(unknown)
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

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			log.Printf("Expected ReferenceLength() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumInsertions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from unmapped []Cigar
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

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			t.Logf("Expected QueryLength() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Cannot calculate NumInsertions from unaligned reads.\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from unmapped []Cigar
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
	var buf bytes.Buffer
	log.SetOutput(&buf)

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			t.Logf("Expected ConsumesReference() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Invalid byte: $\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from Op == '$'
	ConsumesReference('$')
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
	var buf bytes.Buffer
	log.SetOutput(&buf)

	// restore log.Panic(err) from unmapped cigar and check if err is caught
	defer func() {
		log.SetOutput(os.Stderr) // Restore logs
		r := recover()
		if r == nil {
			t.Logf("Expected ConsumesQuery() log.panic, but none occurred.\n")
		}
		if !strings.Contains(buf.String(), "Invalid byte: /\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	// err panic from Op == '/'
	ConsumesQuery('/')
}

func BenchmarkCigarToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ToString([]Cigar{c1, c2, c3})
	}
}

func BenchmarkCigarFromString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		FromString("35M2I16D")
	}
}
