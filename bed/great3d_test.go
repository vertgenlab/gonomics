package bed

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FillSpaceTests = []struct {
	InputFile string
	Genome    map[string]chromInfo.ChromInfo
	OutFile   string
	Expected  string
}{
	{
		InputFile: "testdata/FillSpace.Input.bed",
		Genome:    map[string]chromInfo.ChromInfo{"chr1": {Name: "chr1", Size: 600}, "chr2": {Name: "chr2", Size: 60}},
		OutFile:   "testdata/tmp.FillSpace.bed",
		Expected:  "testdata/FillSpace.Expected.bed",
	},
}

func TestFillSpaceNoHiddenValue(t *testing.T) {
	var err error
	var records, answer []Bed
	for _, v := range FillSpaceTests {
		records = Read(v.InputFile)
		answer = FillSpaceNoHiddenValue(records, v.Genome)
		Write(v.OutFile, answer)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in FillSpace. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var FillThreeDSpaceTests = []struct {
	InputFile string
	Genome    map[string]chromInfo.ChromInfo
	OutFile   string
	Expected  string
}{
	{
		InputFile: "testdata/FillSpace.Hidden.Input.bed",
		Genome:    map[string]chromInfo.ChromInfo{"chr1": {Name: "chr1", Size: 600}, "chr2": {Name: "chr2", Size: 60}},
		OutFile:   "testdata/tmp.Hidden.FillSpace.bed",
		Expected:  "testdata/FillSpace.Hidden.Expected.bed",
	},
}

func TestFillSpaceHiddenValue(t *testing.T) {
	var err error
	var records, answer []Bed
	for _, v := range FillThreeDSpaceTests {
		records = Read(v.InputFile)
		answer = FillSpaceHiddenValue(records, v.Genome)
		Write(v.OutFile, answer)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in FillSpace. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
