package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var IntervalContactsTests = []struct {
	BedPeFile string
	InFile string
	OutFile string
	ExpectedFile string
}{
	{"testdata/contacts.bedpe", "testdata/input.bed", "testdata/test.out.bed", "testdata/expected.out.bed"},
	{"testdata/contacts.bedpe", "testdata/input.vcf", "testdata/test.vcf.out.bed", "testdata/expected.vcf.out.bed"},
}

func TestIntervalContacts(t *testing.T) {
	var err error
	for _, v := range IntervalContactsTests {
		intervalContacts(v.BedPeFile, v.InFile, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in intervalContacts.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
