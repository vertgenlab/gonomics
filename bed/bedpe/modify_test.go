package bedpe

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

func TestContactsToMidpoints(t *testing.T) {
	var records []BedPe
	records = Read("testdata/BedPeFileTest.bedpe")
	contactsToMidpoints(records)
	Write("testdata/temp.ContactMidpoints.bedpe", records)

	equal := fileio.AreEqual("testdata/expectedContactsMidpoints.bedpe", "testdata/temp.ContactMidpoints.bedpe")

	if !equal {
		t.Errorf("Expected bedpe file did not match test output.")
	} else {
		os.Remove("testdata/temp.ContactMidpoints.bedpe")
	}
}
