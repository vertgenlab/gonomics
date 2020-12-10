package fastq

import (
	"os"
	"reflect"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.fastq"},
}

func TestRead(t *testing.T) {
	var tmpFilename string
	for _, test := range readWriteTests {
		tmpFilename = test.filename + ".tmp"
		fastq := Read(test.filename)
		Write(tmpFilename, fastq)
		fastqTwo := Read(tmpFilename)
		if !reflect.DeepEqual(fastq, fastqTwo) {
			t.Errorf("Error: Read,Write,Read was not equal to the Read\n")
		}
		os.Remove(tmpFilename)
	}
}
