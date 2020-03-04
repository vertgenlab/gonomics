package vcf

import (
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.vcf"},
}

//TODO: need to finish writing the write function
func TestWriteAndRead(t *testing.T) {
	var actual []*Vcf
	for _, test := range readWriteTests {
		tempFile := test.filename //+ ".tmp"

		actual = Read(tempFile)

		//Write(tempFile, actual)
		//alpha := ReadFile(tempFile)
		//beta := ReadFile("testdata/test.vcf")
		PrintVcf(actual)
		//if !AllEqual(alpha.Vcf, beta.Vcf) {
		//	t.Errorf("VCF files are not the same")
		//}
		//err := os.Remove(tempFile)
		//if err != nil {
		//	t.Errorf("Deleting temp file %s gave an error.", tempFile)
		//}
	}
}
