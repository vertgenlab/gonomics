package gtf

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var GenesToPromoterBedTests = []struct {
	gFile string
	cFile string
	upstream int
	downstream int
	outFile string
	expectedFile string
}{
	{"testdata/CFTR.test.gtf", "testdata/chr7.chrom.sizes", 500, 2000, "testdata/tmp.bed", "testdata/GenesToPromoter.expected.bed"},
}

func TestGenesToPromoterBed(t *testing.T) {
	var genes map[string]*Gene
	var answer []bed.Bed
	var err error
	var c map[string]chromInfo.ChromInfo
	for _, v := range GenesToPromoterBedTests {
		genes = Read(v.gFile)
		c = chromInfo.ReadToMap(v.cFile)
		answer = GenesToPromoterBed(genes, c, v.upstream, v.downstream)
		bed.Write(v.outFile, answer)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in GenesToPromoterBed, output did not match expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

var GenesToTssBedTests = []struct {
	gFile string
	cFile string
	outFile string
	expectedFile string
}{
	{"testdata/CFTR.test.gtf", "testdata/chr7.chrom.sizes", "testdata/tmp.Tss.bed", "testdata/GenesToTss.expected.bed"},
}

func TestGenesToTssBed(t *testing.T) {
	var genes map[string]*Gene
	var answer []bed.Bed
	var err error
	var c map[string]chromInfo.ChromInfo
	for _, v := range GenesToTssBedTests {
		genes = Read(v.gFile)
		c = chromInfo.ReadToMap(v.cFile)
		answer = GenesToTssBed(genes, c)
		bed.Write(v.outFile, answer)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in GenesToPromoterBed, output did not match expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
