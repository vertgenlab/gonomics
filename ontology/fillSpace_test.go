package ontology

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestFill3dSpace(t *testing.T) {
	var genes = []bed.Bed{{Chrom: "chr1", ChromStart: 2, ChromEnd: 3, Name: "first", Score: 0}, {Chrom: "chr1", ChromStart: 13, ChromEnd: 14, Name: "second", Score: 0}, {Chrom: "chr1", ChromStart: 500, ChromEnd: 501, Name: "third", Score: 0}, {Chrom: "chr2", ChromStart: 10, ChromEnd: 40, Name: "fourth", Score: 0}}
	var contacts = []bedpe.BedPe{{bed.Bed{Chrom: "chr1", ChromStart: 80, ChromEnd: 81, Name: "", Score: 0}, bed.Bed{Chrom: "chr1", ChromStart: 300, ChromEnd: 301, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr2", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr2", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr3", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr3", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}}
	var size = []chromInfo.ChromInfo{{Name: "chr1", Size: 600}, {Name: "chr2", Size: 100}}
	var err error
	answer := Fill3dSpace(contacts, genes, chromInfo.SliceToMap(size))
	bed.Write("testdata/fill3dSpaceOut.bed", answer)

	if !bed.AllAreEqual(answer, bed.Read("testdata/expected.fill3dSpace.bed")) {
		t.Errorf("Error: output didn't match expected file.")
	} else {
		err = os.Remove("testdata/fill3dSpaceOut.bed")
		exception.PanicOnErr(err)
	}
}

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
	var records, answer []bed.Bed
	for _, v := range FillSpaceTests {
		records = bed.Read(v.InputFile)
		answer = FillSpaceNoHiddenValue(records, v.Genome)
		bed.Write(v.OutFile, answer)
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
	var records, answer []bed.Bed
	for _, v := range FillThreeDSpaceTests {
		records = bed.Read(v.InputFile)
		answer = FillSpaceHiddenValue(records, v.Genome)
		bed.Write(v.OutFile, answer)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in FillSpace. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
