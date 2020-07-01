package axt

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/chrM_gasacu1.axt"},
}

func TestReadToChan(t *testing.T) {
	for _, test := range readWriteTests {
		file := Read(test.filename)
		testFile := fileio.EasyOpen(test.filename)
		reader := make(chan *Axt)
		defer testFile.Close()
		go ReadToChan(testFile, reader)
		var index int = 0
		for each := range reader {
			if !isEqual(each, file[index]) {
				t.Errorf("Error: Read to chan function does not equal standard read funtion\n")
			}
			index++
		}

	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Axt
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		if !AllEqual(Read(tempFile), Read("testdata/chrM_gasacu1.axt")) {
			t.Errorf("Axt files are not the same")
		}

		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}

func TestAxtSwap(t *testing.T) {
	var queryLen int64 = 50
	aTest := &Axt{
		RName:      "TargetGenome",
		RStart:     0,
		REnd:       100,
		QName:      "QueryGenome",
		QStart:     20,
		QEnd:       40,
		QStrandPos: false,
	}

	aSwap := QuerySwap(aTest, queryLen)
	log.Printf("%s", ToString(aSwap, 0))
}
