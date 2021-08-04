package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"math"
	"os"
	"testing"
)

var bedFilterTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	MinScore     int
	MaxScore     int
	MinLength    int
	MaxLength    int
	MinStart     int
	MaxStart     int
	MinEnd       int
	MaxEnd       int
	MinNameFloat float64
	MaxNameFloat float64
	Chrom        string
}{
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/tmp.bed",
		ExpectedFile: "testdata/testOut.bed",
		MinScore:     2,
		MaxScore:     10,
		MinLength:    0,
		MaxLength:    numbers.MaxInt,
		MinStart:     5,
		MaxStart:     250,
		MinEnd:       15,
		MaxEnd:       800,
		MinNameFloat: -1 * math.MaxFloat64,
		MaxNameFloat: math.MaxFloat64,
		Chrom:        "chr10",
	},
	{InFile: "testdata/test2.bed",
		OutFile:      "testdata/tmp2.bed",
		ExpectedFile: "testdata/testOut2.bed",
		MinScore:     2,
		MaxScore:     10,
		MinLength:    0,
		MaxLength:    numbers.MaxInt,
		MinStart:     5,
		MaxStart:     250,
		MinEnd:       15,
		MaxEnd:       800,
		MinNameFloat: 2.5,
		MaxNameFloat: 50.0,
		Chrom:        "chr10",
	},
}

func TestBedFilter(t *testing.T) {
	var s Settings
	var err error
	for _, i := range bedFilterTests {
		s = Settings{
			InFile:       i.InFile,
			OutFile:      i.OutFile,
			MinScore:     i.MinScore,
			MaxScore:     i.MaxScore,
			MinLength:    i.MinLength,
			MaxLength:    i.MaxLength,
			MinStart:     i.MinStart,
			MaxStart:     i.MaxStart,
			MinEnd:       i.MinEnd,
			MaxEnd:       i.MaxEnd,
			MinNameFloat: i.MinNameFloat,
			MaxNameFloat: i.MaxNameFloat,
			Chrom:        i.Chrom,
		}
		bedFilter(s)
		if !fileio.AreEqual(i.OutFile, i.ExpectedFile) {
			t.Errorf("Error in bedFilter.")
		} else {
			err = os.Remove(i.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
