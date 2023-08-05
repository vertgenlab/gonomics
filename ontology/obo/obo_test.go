package obo

import (
	"testing"
)

// I/O fidelity test.
func TestOboReadAndWrite(t *testing.T) {
	records1, header1 := Read("testdata/test.obo")
	Write("testdata/out.obo", records1, header1)
	records2, header2 := Read("testdata/out.obo")
	if !EqualHeader(header1, header2) {
		t.Errorf("Error: Obo package failed I/0 fidelity test for header.")
	}
	if !AllAreEqual(records1, records2) {
		t.Errorf("Error: Obo package failed I/O fidelity test for records.")
	}
}
