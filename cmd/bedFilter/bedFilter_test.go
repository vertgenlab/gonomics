package main

import (
	"math"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
)

var BedFilterTests = []struct {
	InFile                string
	OutFile               string
	ExpectedFile          string
	MinScore              int
	MaxScore              int
	MinLength             int
	MaxLength             int
	MinStart              int
	MaxStart              int
	MinEnd                int
	MaxEnd                int
	MinNameFloat          float64
	MaxNameFloat          float64
	NameEquals            string
	NameNotEquals         string
	MaxAnnotationFloat    float64
	MinAnnotationFloat    float64
	AnnotationFilterField int
	Chrom                 string
	SubSet                float64
	SetSeed               int64
}{
	{InFile: "testdata/test.bed",
		OutFile:               "testdata/tmp.bed",
		ExpectedFile:          "testdata/expected.bed",
		MinScore:              0,
		MaxScore:              1000,
		MinLength:             3,
		MaxLength:             1000,
		MinStart:              5,
		MaxStart:              999999,
		MinEnd:                10,
		MaxEnd:                1000010,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		MinAnnotationFloat:    -1 * math.MaxFloat64,
		MaxAnnotationFloat:    math.MaxFloat64,
		AnnotationFilterField: 0,
		Chrom:                 "chr1",
		SubSet:                1.0,
		SetSeed:               0,
	},
	{InFile: "testdata/test.bed",
		OutFile:               "testdata/tmp.bed",
		ExpectedFile:          "testdata/expected.SubSet.bed",
		MinScore:              -1 * numbers.MaxInt,
		MaxScore:              numbers.MaxInt,
		MinLength:             0,
		MaxLength:             numbers.MaxInt,
		MinStart:              0,
		MaxStart:              numbers.MaxInt,
		MinEnd:                0,
		MaxEnd:                numbers.MaxInt,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		MinAnnotationFloat:    -1 * math.MaxFloat64,
		MaxAnnotationFloat:    math.MaxFloat64,
		AnnotationFilterField: 0,
		Chrom:                 "",
		SubSet:                0.5,
		SetSeed:               0,
	},
	{InFile: "testdata/test.annotationFilter.bed",
		OutFile:               "testdata/tmp.annotationFilter.bed",
		ExpectedFile:          "testdata/expected.annotationFilter.bed",
		MinScore:              -1 * numbers.MaxInt,
		MaxScore:              numbers.MaxInt,
		MinLength:             0,
		MaxLength:             numbers.MaxInt,
		MinStart:              0,
		MaxStart:              numbers.MaxInt,
		MinEnd:                0,
		MaxEnd:                numbers.MaxInt,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		MinAnnotationFloat:    -10,
		MaxAnnotationFloat:    10,
		AnnotationFilterField: 0,
		Chrom:                 "",
		SubSet:                1,
		SetSeed:               0,
	},
	{InFile: "testdata/test.annotationFilter.secondField.bed",
		OutFile:               "testdata/tmp.annotationFilter.secondField.bed",
		ExpectedFile:          "testdata/expected.annotationFilter.secondField.bed",
		MinScore:              -1 * numbers.MaxInt,
		MaxScore:              numbers.MaxInt,
		MinLength:             0,
		MaxLength:             numbers.MaxInt,
		MinStart:              0,
		MaxStart:              numbers.MaxInt,
		MinEnd:                0,
		MaxEnd:                numbers.MaxInt,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		MinAnnotationFloat:    -10,
		MaxAnnotationFloat:    10,
		AnnotationFilterField: 1,
		Chrom:                 "",
		SubSet:                1,
		SetSeed:               0,
	},
	{InFile: "testdata/test.nameFilter.bed",
		OutFile:               "testdata/tmp.nameFilter.bed",
		ExpectedFile:          "testdata/expected.nameFilter.bed",
		MinScore:              -1 * numbers.MaxInt,
		MaxScore:              numbers.MaxInt,
		MinLength:             0,
		MaxLength:             numbers.MaxInt,
		MinStart:              0,
		MaxStart:              numbers.MaxInt,
		MinEnd:                0,
		MaxEnd:                numbers.MaxInt,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		NameEquals:            "NameMatch",
		NameNotEquals:         "",
		MinAnnotationFloat:    -1 * numbers.MaxInt,
		MaxAnnotationFloat:    numbers.MaxInt,
		AnnotationFilterField: 1,
		Chrom:                 "",
		SubSet:                1,
		SetSeed:               0,
	},
	{InFile: "testdata/test.nameFilter.nonMatch.bed",
		OutFile:               "testdata/tmp.nameFilter.nonMatch.bed",
		ExpectedFile:          "testdata/expected.nameFilter.nonMatch.bed",
		MinScore:              -1 * numbers.MaxInt,
		MaxScore:              numbers.MaxInt,
		MinLength:             0,
		MaxLength:             numbers.MaxInt,
		MinStart:              0,
		MaxStart:              numbers.MaxInt,
		MinEnd:                0,
		MaxEnd:                numbers.MaxInt,
		MinNameFloat:          -1 * math.MaxFloat64,
		MaxNameFloat:          math.MaxFloat64,
		NameEquals:            "",
		NameNotEquals:         "NameNonMatch",
		MinAnnotationFloat:    -1 * numbers.MaxInt,
		MaxAnnotationFloat:    numbers.MaxInt,
		AnnotationFilterField: 1,
		Chrom:                 "",
		SubSet:                1,
		SetSeed:               0,
	},
}

func TestBedFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BedFilterTests {
		s = Settings{
			InFile:                v.InFile,
			OutFile:               v.OutFile,
			MinScore:              v.MinScore,
			MaxScore:              v.MaxScore,
			MinLength:             v.MinLength,
			MaxLength:             v.MaxLength,
			MinStart:              v.MinStart,
			MaxStart:              v.MaxStart,
			MinEnd:                v.MinEnd,
			MaxEnd:                v.MaxEnd,
			MinNameFloat:          v.MinNameFloat,
			MaxNameFloat:          v.MaxNameFloat,
			NameEquals:            v.NameEquals,
			NameNotEquals:         v.NameNotEquals,
			MinAnnotationFloat:    v.MinAnnotationFloat,
			MaxAnnotationFloat:    v.MaxAnnotationFloat,
			AnnotationFilterField: v.AnnotationFilterField,
			Chrom:                 v.Chrom,
			SubSet:                v.SubSet,
			SetSeed:               v.SetSeed,
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
