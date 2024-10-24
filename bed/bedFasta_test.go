package bed

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ToLowerTests = []struct {
	InFastaFile  string
	InBedFile    string
	OutFile      string
	ExpectedFile string
}{
	{InFastaFile: "testdata/toLower.fa",
		InBedFile:    "testdata/toLower.bed",
		OutFile:      "testdata/toLower.out.fa",
		ExpectedFile: "testdata/expected.toLower.fa"},
}

func TestToLower(t *testing.T) {
	var records []fasta.Fasta
	var regions []Bed
	var file *fileio.EasyWriter
	for _, v := range ToLowerTests {
		records = fasta.Read(v.InFastaFile)
		regions = Read(v.InBedFile)
		ToLower(records, regions)
		file = fileio.EasyCreate(v.OutFile)
		fasta.WriteToFileHandle(file, records, 50)
		err := file.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: in ToLower, OutFile did not match ExpectedFile.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
