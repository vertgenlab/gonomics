package vcf

import (
	"bytes"
	"testing"
)

func TestWriteVcf(t *testing.T) {
	// Create a Vcf struct for testing
	record := Vcf{
		Chr:     "1",
		Pos:     1234567,
		Id:      "rs123456",
		Ref:     "A",
		Alt:     []string{"G", "T"},
		Qual:    99.0,
		Filter:  "PASS",
		Info:    "DP=100",
		Format:  []string{"GT", "AD"},
		Samples: []Sample{{[]int16{0, 1}, []bool{false, false}, []string{""}}},
	}
	expected := "1\t1234567\trs123456\tA\tG,T\t99\tPASS\tDP=100\tGT:AD\t0/1\n"

	// Create a buffer to write to
	var buf bytes.Buffer

	// Write the Vcf record to the buffer
	WriteVcf(&buf, record)

	// Get the result as a string
	result := buf.String()

	// Check if the result matches the expected string
	if result != expected {
		t.Errorf("expected %q, got %q", expected, result)
	}
}

func TestIsVcfFile(t *testing.T) {
	tests := []struct {
		filename string
		expected bool
	}{
		{"test.vcf", true},
		{"test.vcf.gz", true},
		{"estAncestorSequence.fa", false},
	}

	for _, test := range tests {
		result := IsVcfFile(test.filename)
		if result != test.expected {
			t.Errorf("Error: IsVcfFile(%q) = %v; want %v", test.filename, result, test.expected)
		}
	}
}