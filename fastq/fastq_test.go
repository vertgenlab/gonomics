package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"reflect"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.fastq"},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		fastq := Read(test.filename)
		Write(test.filename+".tmp", fastq)
		fastqTwo := Read(test.filename + ".tmp")
		if !reflect.DeepEqual(fastq, fastqTwo) {
			common.ExitIfError(fmt.Errorf("Error: Read,Write,Read was not equal to the Read\n"))
		}
	}
}

func TestReadPair(t *testing.T) {
	fqPair := ReadPairs("testdata/simReads_R1.fq", "testdata/simReads_R2.fq")
	for i := 0; i < len(fqPair); i++ {
		fmt.Printf("%s\n%s\n", fqPair[i].Fwd.Name, fqPair[i].Rev.Name)
	}
}
