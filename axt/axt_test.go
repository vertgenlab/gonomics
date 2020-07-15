package axt

import (
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/dna"
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
	var targetLen int64 = 10
	var queryLen int64 = 10
	aTest := &Axt{
		RName:      "TargetGenome",
		RStart:     1,
		REnd:       10,
		QName:      "QueryGenome",
		QStart:     1,
		QEnd:       10,
		QStrandPos: false,
		//RSeq: dna.StringToBases("CGTCCCTCGA"),
		//QSeq: dna.StringToBases("CGTCCCTCGA"),
	}

	ans := &Axt{
		RName:      "QueryGenome",
		RStart:     1,
		REnd:       10,
		QName:      "TargetGenome",
		QStart:     1,
		QEnd:       10,
		QStrandPos: false,
		//RSeq: dna.StringToBases("CGTCCCTCGA"),
		//QSeq: dna.StringToBases("CGTCCCTCGA"),
	}

	aSwap := SwapBoth(aTest, targetLen, queryLen)
	//TODO: Add logic to test that can check reverse complemented sequences
	if !isEqual(aSwap, ans) {
		log.Printf("%s", ToString(aSwap, 0))
		log.Printf("%s", ToString(ans, 0))
		t.Errorf("Error: Swap axt did not yield the expected value...\n")
	}
}
