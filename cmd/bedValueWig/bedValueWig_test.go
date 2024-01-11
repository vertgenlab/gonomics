package main

import (
	"math"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var bedMaxWigTests = []struct {
	InputBed        string
	InputWig        string
	InputChromSizes string
	OutputFile      string
	ExpectedFile    string
	MinFlag         bool
	AverageFlag     bool
	NormFlag        bool
	NoDataValue     float64
}{
	{InputBed: "testdata/testBed.bed",
		InputWig:        "testdata/startOneStepOne.wig",
		InputChromSizes: "testdata/fake.chrom.sizes",
		OutputFile:      "testdata/testBMWOutput.bed",
		ExpectedFile:    "testdata/testBMWExpected.bed",
		MinFlag:         false,
		AverageFlag:     false,
		NormFlag:        false,
		NoDataValue:     math.MaxFloat64},
	{InputBed: "testdata/testBed.bed",
		InputWig:        "testdata/startOneStepOne.wig",
		InputChromSizes: "testdata/fake.chrom.sizes",
		OutputFile:      "testdata/testBMWOutputNorm.bed",
		ExpectedFile:    "testdata/testBMWExpectedNormFlagStep1.bed",
		MinFlag:         false,
		AverageFlag:     false,
		NormFlag:        true,
		NoDataValue:     math.MaxFloat64},
	{InputBed: "testdata/testBed.bed",
		InputWig:        "testdata/startOneStepOne.wig",
		InputChromSizes: "testdata/fake.chrom.sizes",
		OutputFile:      "testdata/testMinOutput.bed",
		ExpectedFile:    "testdata/testMinExpected.bed",
		MinFlag:         true,
		AverageFlag:     false,
		NormFlag:        false,
		NoDataValue:     math.MaxFloat64},
	{InputBed: "testdata/testBed.bed",
		InputWig:        "testdata/startOneStepOne.wig",
		InputChromSizes: "testdata/fake.chrom.sizes",
		OutputFile:      "testdata/testAverageOutput.bed",
		ExpectedFile:    "testdata/testAverageExpected.bed",
		MinFlag:         false,
		AverageFlag:     true,
		NormFlag:        false,
		NoDataValue:     math.MaxFloat64},
	{InputBed: "testdata/testBed.bed",
		InputWig:        "testdata/testNoValue.wig",
		InputChromSizes: "testdata/fake.chrom.sizes",
		OutputFile:      "testdata/testNoData.bed",
		ExpectedFile:    "testdata/testNoDataExpected.bed",
		MinFlag:         false,
		AverageFlag:     false,
		NormFlag:        false,
		NoDataValue:     -10},
}

func TestBedValueWig(t *testing.T) {
	var err error
	for _, v := range bedMaxWigTests {
		s := Settings{
			Infile:      v.InputBed,
			WigFile:     v.InputWig,
			SizesFile:   v.InputChromSizes,
			OutFile:     v.OutputFile,
			MinFlag:     v.MinFlag,
			AverageFlag: v.AverageFlag,
			NormFlag:    v.NormFlag,
			NoDataValue: v.NoDataValue,
		}
		bedValueWig(s)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in bedMaxWig, the output beds is not as expected")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}
