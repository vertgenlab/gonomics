package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var bedCountBamTests = []struct {
	inBam     string
	inBed     string
	outFile   string
	expFile   string
	normalize bool
}{{
	"testdata/in.sort.bam", "testdata/in.bed", "testdata/out.txt", "testdata/exp.txt", false,
}, {
	"testdata/in.sort.bam", "testdata/in.bed", "testdata/out.norm.txt", "testdata/exp.norm.txt", true,
}}

func TestBedCountBam(t *testing.T) {
	var err error
	var s settings
	for _, v := range bedCountBamTests {
		s = settings{
			inBam:   v.inBam,
			inBed:   v.inBed,
			outFile: v.outFile,
			norm:    v.normalize,
		}
		bedCountBam(s)
		if !fileio.AreEqual(v.outFile, v.expFile) {
			t.Errorf("Error in bedCountBam. Expected and output files are not the same. Exp: %s, out %s", v.expFile, v.outFile)
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
