package obo

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
)

// I/O fidelity test.
func TestOboReadAndWrite(t *testing.T) {
	var failed bool = false
	records1, header1 := Read("testdata/test.obo", true)
	Write("testdata/out.obo", records1, header1)
	records2, header2 := Read("testdata/out.obo", true)
	if !EqualHeader(header1, header2) {
		failed = true
		t.Errorf("Error: Obo package failed I/0 fidelity test for header.")
	}
	if !AllAreEqual(records1, records2) {
		failed = true
		t.Errorf("Error: Obo package failed I/O fidelity test for records.")
	}
	if !failed {
		err := os.Remove("testdata/out.obo")
		exception.PanicOnErr(err)
	}
}
