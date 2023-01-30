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
	gFile             string
	cFile             string
	upstream          int
	downstream        int
	outFile           string
	expectedFile      string
	expectedCanonFile string
}{
	{gFile: "testdata/CFTR.test.gtf",
		cFile:             "testdata/chr7.chrom.sizes",
		upstream:          500,
		downstream:        2000,
		outFile:           "testdata/tmp.bed",
		expectedFile:      "testdata/GenesToPromoter.expected.bed",
		expectedCanonFile: "testdata/GenesToPromoter.expected.bed"},
}

func TestGenesToCanonicalBeds(t *testing.T) {
	var genes map[string]*Gene
	var answer []bed.Bed
	var err error
	var c map[string]chromInfo.ChromInfo
	for _, v := range GenesToPromoterBedTests {
		genes = Read(v.gFile)
		c = chromInfo.ReadToMap(v.cFile)
		answer = GenesToCanonicalBeds(genes, c, v.upstream, v.downstream)
		bed.SortByCoord(answer)
		bed.Write(v.outFile, answer)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in GenesToCanonicalTranscriptsBed.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
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
		bed.SortByCoord(answer)
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
	gFile             string
	cFile             string
	outFile           string
	expectedFile      string
	expectedCanonFile string
}{
	{gFile: "testdata/CFTR.test.gtf",
		cFile:             "testdata/chr7.chrom.sizes",
		outFile:           "testdata/tmp.Tss.bed",
		expectedFile:      "testdata/GenesToTss.expected.bed",
		expectedCanonFile: "testdata/GenesToTss.expected.bed"},
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
		bed.SortByCoord(answer)
		bed.Write(v.outFile, answer)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in GenesToPromoterBed, output did not match expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

func TestGenesToCanonicalTranscriptsTssBed(t *testing.T) {
	var genes map[string]*Gene
	var answer []bed.Bed
	var err error
	var c map[string]chromInfo.ChromInfo
	for _, v := range GenesToTssBedTests {
		genes = Read(v.gFile)
		c = chromInfo.ReadToMap(v.cFile)
		answer = GenesToCanonicalTranscriptsTssBed(genes, c)
		bed.SortByCoord(answer)
		bed.Write(v.outFile, answer)
		if !fileio.AreEqual(v.outFile, v.expectedCanonFile) {
			t.Errorf("Error in GenesToCanonicalTranscriptsTssBed, output did not match expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
