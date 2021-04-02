package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"reflect"
	"sync"
	"testing"
)

var r001 = Sam{
	QName: "r001",
	Flag:  99,
	MapQ:  30,
	RName: "ref",
	Pos:   7,
	Cigar: cigar.FromString("8M2I4M1D3M"),
	RNext: "=",
	PNext: 37,
	TLen:  39,
	Seq:   dna.StringToBases("TTAGATAAAGGATACTG"),
	Qual:  "*",
	Extra: "",
}

var r002 = Sam{
	QName: "r002",
	Flag:  0,
	MapQ:  30,
	RName: "ref",
	Pos:   9,
	Cigar: cigar.FromString("3S6M1P1I4M"),
	RNext: "*",
	PNext: 0,
	TLen:  0,
	Seq:   dna.StringToBases("AAAAGATAAGGATA"),
	Qual:  "*",
	Extra: "",
}

var r003 = Sam{
	QName: "r003",
	Flag:  0,
	MapQ:  30,
	RName: "ref",
	Pos:   9,
	Cigar: cigar.FromString("5S6M"),
	RNext: "*",
	PNext: 0,
	TLen:  0,
	Seq:   dna.StringToBases("GCCTAAGCTAA"),
	Qual:  "*",
	Extra: "SA:Z:ref,29,-,6H5M,17,0;",
}

var r004 = Sam{
	QName: "r004",
	Flag:  0,
	MapQ:  30,
	RName: "ref",
	Pos:   16,
	Cigar: cigar.FromString("6M14N5M"),
	RNext: "*",
	PNext: 0,
	TLen:  0,
	Seq:   dna.StringToBases("ATAGCTTCAGC"),
	Qual:  "*",
	Extra: "",
}

var r003Supplemental = Sam{
	QName: "r003",
	Flag:  2064,
	MapQ:  17,
	RName: "ref",
	Pos:   29,
	Cigar: cigar.FromString("6H5M"),
	RNext: "*",
	PNext: 0,
	TLen:  0,
	Seq:   dna.StringToBases("TAGGC"),
	Qual:  "*",
	Extra: "SA:Z:ref,9,+,5S6M,30,1;",
}

var r001Supplemental = Sam{
	QName: "r001",
	Flag:  147,
	MapQ:  30,
	RName: "ref",
	Pos:   37,
	Cigar: cigar.FromString("9M"),
	RNext: "=",
	PNext: 7,
	TLen:  -39,
	Seq:   dna.StringToBases("CAGCGGCAT"),
	Qual:  "*",
	Extra: "NM:i:1",
}

func TestRead(t *testing.T) {
	actual, _ := Read("testdata/small.sam")

	if len(actual) != 6 {
		t.Error("problem reading sam")
	} else {
		expected := []Sam{r001, r002, r003, r004, r003Supplemental, r001Supplemental}
		for i := range actual {
			if !Equal(actual[i], expected[i]) {
				t.Error("problem reading sam")
			}
		}
	}

}

func TestReadNextIntactHeader(t *testing.T) {
	file := fileio.EasyOpen("testdata/small.sam")
	expected := []Sam{r001, r002, r003, r004, r003Supplemental, r001Supplemental}
	var i int
	for val, done := ReadNext(file); !done; val, done = ReadNext(file) {
		if !Equal(val, expected[i]) {
			t.Error("problem with sam.ReadNext not skipping header")
		}
		i++
	}
}

func TestReadToChan(t *testing.T) {
	samChan, _ := GoReadToChan("testdata/small.sam")
	expected := []Sam{r001, r002, r003, r004, r003Supplemental, r001Supplemental}
	var i int
	for val := range samChan {
		if !Equal(val, expected[i]) {
			t.Error("problem reading sam to channel")
		}
		i++
	}
}

func TestWriteFromChan(t *testing.T) {
	data, header := GoReadToChan("testdata/small.sam")
	var wg sync.WaitGroup
	wg.Add(1)
	WriteFromChan(data, "testdata/small.sam.tmp", header, &wg)
	wg.Wait()
	if !fileio.AreEqual("testdata/small.sam", "testdata/small.sam.tmp") {
		t.Error("problem with WriteFromChan")
	}
	err := os.Remove("testdata/small.sam.tmp")
	if err != nil {
		t.Errorf("Deleting temp file %s gave an error.", "testdata/small.sam.tmp")
	}
}

func TestGenerateHeader(t *testing.T) {
	_, expected := Read("testdata/small.sam")
	actual := GenerateHeader(expected.Chroms, expected.Text[2:], Coordinate, Reference)
	if len(actual.Text) != len(expected.Text) {
		t.Error("problem with generate header")
	}

	for i := range actual.Text {
		if actual.Text[i] != expected.Text[i] {
			t.Error("problem with generate header")
		}
	}

	// messy check to make sure that they were parsed equally
	// should always be true since the input text was identical
	// per the above check
	if !reflect.DeepEqual(actual, expected) {
		t.Error("generate header was parsed differently from the expected header")
	}
}

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.sam"},
}

// testing for runtime errors on real files
func TestReadReal(t *testing.T) {
	for _, test := range readWriteTests {
		Read(test.filename)
	}
}

// testing for runtime errors on real files
func TestReadAndWriteReal(t *testing.T) {
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual, header := Read(test.filename)
		Write(tempFile, actual, header)
		if !fileio.AreEqual(tempFile, test.filename) {
			t.Errorf("problem writing file")
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
