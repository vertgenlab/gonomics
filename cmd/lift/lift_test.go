package main

import (
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/bed"
	"testing"
	"os"
)

var LiftTests = []struct {
	inputFile string
	expectedOutputFile string
	chainFile string
	faFile string
}{
	{"testdata/input.bed", "testdata/expected.bed", "testdata/test.chain", ""},
	{"testdata/Pollard.HARs.hg19.trimmed.bed", "testdata/Pollard.HARs.hg38.UCSC.trimmed.bed", "testdata/hg19ToHg38.over.chain", ""},
}

func TestLift(t *testing.T) {
	for _, v := range LiftTests {
		if vcf.IsVcfFile(v.inputFile) {
			lift(v.chainFile, v.inputFile, "tmp.vcf", v.faFile)
			records := vcf.Read("tmp.vcf")
			expected := vcf.Read(v.expectedOutputFile)
			if !vcf.AllEqual(records, expected) {
				t.Errorf("Error in Lift for vcf. Input: %+v. Expected: %+v.", records, expected)
			}
			err := os.Remove("tmp.vcf")
			if err != nil {
				common.ExitIfError(err)
			}
		} else {
			lift(v.chainFile, v.inputFile, "tmp.bed", v.faFile)
			records := bed.Read("tmp.bed")
			expected := bed.Read(v.expectedOutputFile)
			if !bed.AllAreEqual(records, expected) {
				t.Errorf("Error in Lift for bed. Input: %+v. Expected: %+v.", records, expected)
			}
			//err := os.Remove("tmp.bed")
			//if err != nil {
			//	common.ExitIfError(err)
			//}
		}
	}
}