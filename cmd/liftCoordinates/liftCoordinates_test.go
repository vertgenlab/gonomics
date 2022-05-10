package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var LiftTests = []struct {
	inputFile          string
	expectedOutputFile string
	chainFile          string
	faFile             string
	verbose            int
}{
	//{"testdata/input.bed", "testdata/expected.bed", "testdata/test.chain", ""},
	//{"testdata/Pollard.HARs.hg19.trimmed.bed", "testdata/Pollard.HARs.hg38.UCSC.trimmed.bed", "testdata/hg19ToHg38.over.chain", ""},
	{"testdata/input.vcf", "testdata/expected.vcf", "testdata/test.chain", "testdata/test.fa", 0},
}

func TestLift(t *testing.T) {
	for _, v := range LiftTests {
		liftCoordinates(v.chainFile, v.inputFile, "tmp.vcf", v.faFile, "tmp.unmapped", 0.95, false, 0)

		if vcf.IsVcfFile(v.inputFile) {
			liftCoordinates(v.chainFile, v.inputFile, "tmp.vcf", v.faFile, "tmp.unmapped", 0.95, false, 0)
			records, _ := vcf.Read("tmp.vcf")
			expected, _ := vcf.Read(v.expectedOutputFile)
			if !vcf.AllEqual(records, expected) {
				fmt.Println("TEST VALUES")
				fmt.Println(records)
				fmt.Println(expected)
				t.Errorf("Error in Lift for vcf.")
			}
			err := os.Remove("tmp.vcf")
			if err != nil {
				common.ExitIfError(err)
			}
			err = os.Remove("tmp.unmapped")
			if err != nil {
				common.ExitIfError(err)
			}
		} else {
			liftCoordinates(v.chainFile, v.inputFile, "tmp.bed", v.faFile, "tmp.unmapped", 0.95, false, 0)
			records := bed.Read("tmp.bed")
			expected := bed.Read(v.expectedOutputFile)
			if !bed.AllAreEqual(records, expected) {
				t.Errorf("Error in Lift for bed.")
			}
			err := os.Remove("tmp.bed")
			if err != nil {
				common.ExitIfError(err)
			}
			err = os.Remove("tmp.unmapped")
			if err != nil {
				common.ExitIfError(err)
			}
		}
	}
}
