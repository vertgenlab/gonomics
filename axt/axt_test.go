package axt

import (
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
	var index int
	for _, test := range readWriteTests {
		file := Read(test.filename)
		reader, headerLines := GoReadToChan(test.filename)
		if len(headerLines) != 15 || headerLines[0] != "# header test start" || headerLines[len(headerLines)-1] != "#   header test end" {
			t.Errorf("Error: header was not as expected\n")
		}
		index = 0
		for each := range reader {
			if !isEqual(each, file[index]) {
				t.Errorf("Error: Read to chan function does not equal standard read function\n")
			}
			index++
		}

	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []Axt
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		if !allEqual(Read(tempFile), Read("testdata/chrM_gasacu1.axt")) {
			t.Errorf("Axt files are not the same")
		}

		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}

//TODO: Finish rev. comp sequences for testing swap
func TestAxtSwap(t *testing.T) {
	var targetLen int = 10
	var queryLen int = 10
	aTest := Axt{
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

	ans := Axt{
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

	Swap(&aTest, targetLen, queryLen)
	//TODO: Add logic to test that can check reverse complemented sequences
	if !isEqual(aTest, ans) {
		log.Printf("%s", ToString(aTest, 0))
		log.Printf("%s", ToString(ans, 0))
		t.Errorf("Error: Swap axt did not yield the expected value...\n")
	}
}
