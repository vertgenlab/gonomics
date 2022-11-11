package bedpe

import (
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestContactsToMidpoints(t *testing.T) {
	records := Read("testdata/BedPeFileTest.bedpe")
	out := contactsToMidpoints(records)
	Write("testdata/temp.ContactMidpoints.bedpe", out)

	equal := fileio.AreEqual("testdata/expectedContactsMidpoints.bedpe", "testdata/temp.ContactMidpoints.bedpe")

	if !equal {
		t.Errorf("Expected bedpe file did not match test output.")
	} else {
		os.Remove("testdata/temp.ContactMidpoints.bedpe")
	}

}
