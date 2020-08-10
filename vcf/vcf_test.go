package vcf

import (
	"log"
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.vcf"},
}

func TestReadToChan(t *testing.T) {
	alpha := Read("testdata/test.vcf")
	var beta []*Vcf
	vcfPipe, _ := GoReadToChan("testdata/test.vcf")

	for vcfs := range vcfPipe {
		beta = append(beta, vcfs)
	}
	if AllEqual(alpha, beta) {
		log.Printf("PASS: go ReadToChan function matches standard read funtion\n")
	} else {
		t.Errorf("Error there might be a bug in the ReadToChan implementation\n")
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Vcf
	for _, test := range readWriteTests {
		tempFile := "tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		alpha := Read(tempFile)
		beta := Read("testdata/test.vcf")
		log.Printf("Looks good to me!\n")
		log.Printf("alpha=%d, beta=%d", len(alpha), len(beta))
		if !AllEqual(alpha, beta) {
			t.Errorf("Error: Read and write files do not match\n")
		}
		os.Remove(tempFile)
	}
}

func TestReadToChanTwo(t *testing.T) {
	alpha := GoReadGVcf("testdata/test.vcf")
	var savedFromAlpha []*Vcf
	for v := range alpha.Vcfs {
		savedFromAlpha = append(savedFromAlpha, v)
	}
	alpha.File.Close()
	vcfData, _ := GoReadToChan("testdata/test.vcf")
	var i int = 0
	for each := range vcfData {
		if !isEqual(each, savedFromAlpha[i]) {
			t.Errorf("Error: Read channels are not matching\n")
		}
		i++
	}
}
