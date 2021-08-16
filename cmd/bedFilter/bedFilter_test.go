package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"os"
	"testing"
)

var BedFilterTests = []struct {
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
	Chrom        string
	SubSet       float64
	RandSeed     bool
	SetSeed      int64
}{
	{"testdata/test.bed",
		"testdata/tmp.bed",
		"testdata/expected.bed",
		0,
		1000,
		3,
		1000,
		5,
		999999,
		10,
		1000010,
		"chr1",
		1.0,
		false,
		0,
	},
	{"testdata/test.bed",
		"testdata/tmp.bed",
		"testdata/expected.SubSet.bed",
		-1 * numbers.MaxInt,
		numbers.MaxInt,
		0,
		numbers.MaxInt,
		0,
		numbers.MaxInt,
		0,
		numbers.MaxInt,
		"",
		0.5,
		false,
		0,
	},
}

func TestBedFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BedFilterTests {
		s = Settings{
			InFile:    v.InFile,
			OutFile:   v.OutFile,
			MinScore:  v.MinScore,
			MaxScore:  v.MaxScore,
			MinLength: v.MinLength,
			MaxLength: v.MaxLength,
			MinStart:  v.MinStart,
			MaxStart:  v.MaxStart,
			MinEnd:    v.MinEnd,
			MaxEnd:    v.MaxEnd,
			Chrom:     v.Chrom,
			SubSet:    v.SubSet,
			RandSeed:  v.RandSeed,
			SetSeed:   v.SetSeed,
		}
		bedFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedFilter.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
